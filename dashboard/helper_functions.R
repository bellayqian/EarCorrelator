library(dplyr)
library(tidyr)
library(nnet)
library(matrixcalc)
library(MASS)
library(Matrix)
library(caret)
library(lme4)
library(nlme)
set.seed(1)

## Workflow 
extract_variance_components <- function(mix_model) {
  # Get variance components from model
  vc <- VarCorr(mix_model)
  var_residual <- attr(vc, "sc")^2                # Residual variance
  var_individual <- as.numeric(vc$SID.4[1, 1])    # Individual-level variance (SID)
  var_ear <- as.numeric(vc$SID.EAR[1, 1])         # Ear-level variance (SID:EAR)
  
  # Extract frequency variances for each phenotype
  frequency_order <- c("250", "500", "1K", "2K", "3K", "4K", "6K", "8K")
  
  # Initialize lists to store covariance components by phenotype
  phenotype_cov_matrices <- list()
  
  # Extract variances for each phenotype
  for (phenotype in 1:4) {
    group_name <- switch(phenotype,
                         "SID.3",  # phenotype 1 -> SID.3
                         "SID.2",  # phenotype 2 -> SID.2
                         "SID.1",  # phenotype 3 -> SID.1
                         "SID")    # phenotype 4 -> SID
    
    if (group_name %in% names(vc)) {
      # Get the relevant variance component matrix
      group_vc <- vc[[group_name]]
      
      # Create a properly ordered covariance matrix
      n_freqs <- length(frequency_order)
      ordered_cov <- matrix(0, nrow = n_freqs, ncol = n_freqs)
      rownames(ordered_cov) <- colnames(ordered_cov) <- frequency_order
      
      # Fill the covariance matrix with values in the correct order
      for (i in seq_along(frequency_order)) {
        for (j in seq_along(frequency_order)) {
          freq_i <- frequency_order[i]
          freq_j <- frequency_order[j]
          
          row_name <- paste0("subtype", phenotype, ":Frequency", freq_i)
          col_name <- paste0("subtype", phenotype, ":Frequency", freq_j)
          
          if (row_name %in% rownames(group_vc) && col_name %in% colnames(group_vc)) {
            ordered_cov[i, j] <- as.numeric(group_vc[row_name, col_name])
          } 
        }
      }
      phenotype_cov_matrices[[phenotype]] <- ordered_cov
    } else {
      warning(paste("Could not find group", group_name, "in variance components"))
      # Create empty matrix with correct dimensions
      empty_mat <- matrix(0, length(frequency_order), length(frequency_order))
      rownames(empty_mat) <- colnames(empty_mat) <- frequency_order
      phenotype_cov_matrices[[phenotype]] <- empty_mat
    }
  }
  
  return(list(
    var_individual = var_individual,
    var_ear = var_ear,
    var_residual = var_residual,
    freq_variances = phenotype_cov_matrices
  ))
}

build_cov_matrix <- function(label_idx1, label_idx2, cov_components) {
  # Extract variance components
  var_individual <- cov_components$var_individual
  var_ear <- cov_components$var_ear
  var_residual <- cov_components$var_residual
  
  # Define frequencies
  freqs <- c("250", "500", "1K", "2K", "3K", "4K", "6K", "8K")
  n_freqs <- length(freqs)
  
  # For mixed-model approach, we need 16x16 matrix to match combined ear data
  matrix_size <- n_freqs * 2
  
  # Create the full covariance matrix (16x16)
  cov_matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  
  # Create meaningful row/column names
  rownames(cov_matrix) <- colnames(cov_matrix) <- c(paste0("Left_", freqs), paste0("Right_", freqs))
  
  # Determine if same phenotype or different phenotypes
  same_phenotype <- (label_idx1 == label_idx2)
  
  # Get the covariance matrices for each phenotype
  cov_matrix1 <- cov_components$freq_variances[[label_idx1]]
  cov_matrix2 <- cov_components$freq_variances[[label_idx2]]
  
  # 1. Upper-left block (Left ear - phenotype label_idx1)
  for (i in 1:n_freqs) {
    for (j in 1:n_freqs) {
      if (i == j) {
        # Diagonal: var_individual + var_ear + freq_variance from matrix + var_residual
        cov_matrix[i, j] <- var_individual + var_ear + cov_matrix1[i, i] + var_residual
      } else {
        # Off-diagonal: var_individual + var_ear + covariance from matrix
        cov_matrix[i, j] <- var_individual + var_ear + cov_matrix1[i, j]
      }
    }
  }
  
  # 2. Lower-right block (Right ear - phenotype label_idx2)
  for (i in 1:n_freqs) {
    for (j in 1:n_freqs) {
      row_idx <- i + n_freqs
      col_idx <- j + n_freqs
      
      if (i == j) {
        # Diagonal: var_individual + var_ear + freq_variance from matrix + var_residual
        cov_matrix[row_idx, col_idx] <- var_individual + var_ear + cov_matrix2[i, i] + var_residual
      } else {
        # Off-diagonal: var_individual + var_ear + covariance from matrix
        cov_matrix[row_idx, col_idx] <- var_individual + var_ear + cov_matrix2[i, j]
      }
    }
  }
  
  # 3. Cross-ear blocks (Upper-right and Lower-left)
  if (same_phenotype) {
    # Same phenotype - use covariance structure from matrix
    for (i in 1:n_freqs) {
      for (j in 1:n_freqs) {
        if (i == j) {
          # Same frequency, same phenotype: var_individual + freq_variance
          cov_matrix[i, j + n_freqs] <- var_individual + cov_matrix1[i, i]
          cov_matrix[i + n_freqs, j] <- var_individual + cov_matrix1[i, i]
        } else {
          # Different frequencies: var_individual + covariance
          cov_matrix[i, j + n_freqs] <- var_individual + cov_matrix1[i, j]
          cov_matrix[i + n_freqs, j] <- var_individual + cov_matrix1[i, j]
        }
      }
    }
  } else {
    # Different phenotypes - only individual effects shared
    for (i in 1:n_freqs) {
      for (j in 1:n_freqs) {
        # All cross-ear elements just have var_individual when phenotypes differ
        cov_matrix[i, j + n_freqs] <- var_individual
        cov_matrix[i + n_freqs, j] <- var_individual
      }
    }
  }
  
  return(cov_matrix)
}

calculate_initial_estimates <- function(data, frequencies) {
  # Calculate mean vectors for each phenotype (theta)
  initial_theta <- data %>%
    group_by(Label) %>%
    summarise(across(all_of(frequencies), mean)) %>%
    arrange(Label) %>%
    ungroup() %>%
    dplyr::select(-Label) %>%
    as.matrix()
  
  # Calculate empirical SDs for each phenotype
  empirical_sds <- data %>%
    group_by(Label) %>%
    summarise(across(all_of(frequencies), sd)) %>%
    arrange(Label) %>%
    ungroup() %>%
    dplyr::select(-Label) %>%
    as.matrix()
  
  # Calculate marginal probabilities
  labels <- sort(unique(data$Label))
  marginal_probs <- prop.table(table(factor(data$Label, levels = labels)))
  
  # Return estimates
  list(
    initial_theta = initial_theta,
    empirical_sds = empirical_sds,
    initial_omega = marginal_probs
  )
}

dmvnorm <- function(x, mean, sigma, log = FALSE) {
  k <- length(x)
  
  if (any(is.na(x)) || any(is.na(mean)) || any(is.na(sigma))) {
    return(if(log) -Inf else 0)
  }
  
  # Add small constant to diagonal for numerical stability
  # Use a small fraction of the maximum diagonal element to preserve matrix structure
  max_diag <- max(diag(sigma), na.rm = TRUE)
  if (!is.finite(max_diag) || max_diag <= 0) max_diag <- 1e-6
  
  sigma_reg <- sigma + diag(1e-8 * max_diag, k)
  
  # Check if resulting matrix is positive definite
  eig <- tryCatch(
    eigen(sigma_reg, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(-1, k)
  )
  
  if (any(eig <= 0)) {
    # Add more regularization if needed
    sigma_reg <- sigma_reg + diag(1e-6, k)
  }
  
  # Use Cholesky decomposition for better numerical stability
  chol_result <- tryCatch(
    chol(sigma_reg),
    error = function(e) {
      # If Cholesky fails, use nearPD with conservative settings
      mat <- nearPD(sigma_reg, corr = FALSE, keepDiag = TRUE)
      chol(mat$mat)
    }
  )
  
  # Compute quadratic form using Cholesky decomposition
  diff <- x - mean
  quad_form <- sum((backsolve(chol_result, diff, transpose = TRUE))^2)
  
  # Compute log determinant from Cholesky
  log_det <- 2 * sum(log(diag(chol_result)))
  
  # Calculate log density
  logdens <- -0.5 * (k * log(2 * pi) + log_det + quad_form)
  
  if (log) return(logdens)
  return(exp(logdens))
}

evaluate_accuracy <- function(results, data, permutation_invariant = TRUE) {
  # Create evaluation dataframe
  eval_df <- data.frame(
    SID = character(),
    EAR = character(),
    True = numeric(),
    Predicted = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process results into dataframe
  for (sid in names(results$classifications)) {
    subject_data <- data[data$SID == sid, ]
    predictions <- results$classifications[[sid]]
    
    # Handle both single and paired ear cases
    if (length(predictions) == 1) {
      # Single ear case
      eval_df <- rbind(eval_df, data.frame(
        SID = sid,
        EAR = subject_data$EAR,
        True = subject_data$Label,
        Predicted = predictions
      ))
    } else {
      subject_data <- subject_data[order(subject_data$EAR), ]
      
      for (i in 1:length(predictions)) {
        eval_df <- rbind(eval_df, data.frame(
          SID = sid,
          EAR = subject_data$EAR[i],
          True = subject_data$Label[i],
          Predicted = predictions[i]
        ))
      }
    }
  }
  
  # Calculate ear-level accuracy (standard)
  ear_accuracy <- mean(eval_df$True == eval_df$Predicted, na.rm = TRUE)
  
  # For person-level accuracy
  if (permutation_invariant) {
    # Group by subject ID
    person_correct <- sapply(split(eval_df, eval_df$SID), function(subj_data) {
      n_ears <- nrow(subj_data)
      
      if (n_ears == 1) {
        # Single ear case - direct comparison
        return(all(subj_data$True == subj_data$Predicted, na.rm = TRUE))
      } else if (n_ears == 2) {
        # For paired ears, check both normal and flipped assignments
        true_labels <- subj_data$True
        pred_labels <- subj_data$Predicted
        
        # Option 1: Direct match (L→L, R→R)
        exact_match <- all(true_labels == pred_labels, na.rm = TRUE)
        
        # Option 2: Flipped match (L→R, R→L)
        flipped_match <- all(true_labels == rev(pred_labels), na.rm = TRUE)
        
        return(exact_match || flipped_match)
      } else {
        # For more than 2 ears (unusual case)
        # Check if multisets are the same (same labels regardless of order)
        return(all(sort(subj_data$True) == sort(subj_data$Predicted), na.rm = TRUE))
      }
    })
    
    person_accuracy <- mean(person_correct, na.rm = TRUE)
  } else {
    # Standard person-level accuracy calculation (no permutation invariance)
    person_results <- eval_df %>%
      group_by(SID) %>%
      summarize(
        correct = all(True == Predicted, na.rm = TRUE),
        n_ears = n()
      )
    
    person_accuracy <- mean(person_results$correct, na.rm = TRUE)
  }
  
  # Calculate counts by ear number
  ear_counts <- eval_df %>%
    group_by(SID) %>%
    summarize(n_ears = n())
  
  n_single <- sum(ear_counts$n_ears == 1)
  n_paired <- sum(ear_counts$n_ears == 2)
  
  list(
    ear_accuracy = ear_accuracy,
    person_accuracy = person_accuracy,
    permutation_invariant = permutation_invariant,
    counts = list(
      single_ear_subjects = n_single,
      paired_ear_subjects = n_paired
    ),
    detailed_results = eval_df
  )
}

format_predictions <- function(predictions_df, pred_column = "naive_pred") {
  # Split data by SID
  sids <- unique(predictions_df$SID)
  classifications_list <- list()
  probabilities_list <- list()
  
  n_classes <- 4
  
  for(sid in sids) {
    # Get data for this subject and ensure it's ordered by EAR
    subject_data <- predictions_df[predictions_df$SID == sid, ]
    subject_data <- subject_data[order(subject_data$EAR), ]
    n_ears <- nrow(subject_data)
    
    # Extract predictions and ensure they're numeric
    predictions <- as.numeric(subject_data[[pred_column]])
    
    if(n_ears == 1) {
      # Single ear case
      classifications_list[[as.character(sid)]] <- predictions
      probabilities_list[[as.character(sid)]] <- matrix(0, nrow = 1, ncol = n_classes)
      probabilities_list[[as.character(sid)]][1, predictions] <- 1.0
      
    } else {
      # Paired ear case
      classifications_list[[as.character(sid)]] <- predictions
      prob_matrix <- matrix(0, nrow = n_ears, ncol = n_classes)
      for(i in 1:n_ears) {
        prob_matrix[i, predictions[i]] <- 1.0
      }
      probabilities_list[[as.character(sid)]] <- prob_matrix
    }
  }
  
  # Add names to the classification vectors for clarity
  for(sid in names(classifications_list)) {
    subject_data <- predictions_df[predictions_df$SID == as.numeric(sid), ]
    names(classifications_list[[sid]]) <- subject_data$EAR[order(subject_data$EAR)]
  }
  
  return(list(
    classifications = classifications_list,
    probabilities = probabilities_list
  ))
}

naive_QDA <- function(test_data, theta, frequencies, train_cov_matrices = NULL, train_data = NULL) {
  n_phenotypes <- nrow(theta)
  all_classifications <- list()
  all_probabilities <- list()
  
  # If train_cov_matrices not provided, calculate from training data
  if (is.null(train_cov_matrices)) {
    if (is.null(train_data)) {
      stop("Either train_cov_matrices or train_data must be provided")
    }
    # Calculate covariance matrices using only training data
    train_cov_matrices <- lapply(1:n_phenotypes, function(k) {
      phenotype_data <- train_data[train_data$Label == k, frequencies]
      cov(as.matrix(phenotype_data))
    })
  }
  
  # Calculate prior probabilities from training data
  priors <- table(train_data$Label) / nrow(train_data)
  
  # Process test data by subject
  grouped_data <- split(test_data, test_data$SID)
  
  for (sid in names(grouped_data)) {
    group <- grouped_data[[sid]]
    group <- group[order(group$EAR), ]
    
    # Initialize classification results for this subject
    classifications <- numeric(nrow(group))
    probabilities <- matrix(0, nrow(group), n_phenotypes)
    
    # Process each ear independently
    for (i in 1:nrow(group)) {
      y_i <- as.numeric(group[i, frequencies])
      log_probs <- numeric(n_phenotypes)
      
      # Calculate discriminant function for each phenotype
      for (k in 1:n_phenotypes) {
        mu_k <- theta[k, ]
        sigma_k <- train_cov_matrices[[k]] + diag(1e-6, length(frequencies))
        log_probs[k] <- dmvnorm(y_i, mu_k, sigma_k, log = TRUE) + log(priors[k])
      }
      
      # Convert to probabilities
      max_log_prob <- max(log_probs)
      probs <- exp(log_probs - max_log_prob)
      probs <- probs / sum(probs)
      
      classifications[i] <- which.max(probs)
      probabilities[i, ] <- probs
    }
    
    all_classifications[[sid]] <- classifications
    all_probabilities[[sid]] <- probabilities
  }
  
  return(list(
    classifications = all_classifications,
    probabilities = all_probabilities
  ))
}

improved_QDA <- function(test_data, theta, empirical_sds, frequencies, cov_components, train_data = NULL) {
  n_phenotypes <- nrow(theta)
  all_classifications <- list()
  all_probabilities <- list()
  
  # Use naive QDA for single-ear cases, so calculate individual covariance matrices
  naive_cov_matrices <- lapply(1:n_phenotypes, function(k) {
    phenotype_data <- train_data[train_data$Label == k, frequencies]
    cov(as.matrix(phenotype_data))
  })
  
  # Calculate phenotype prior probabilities from training data
  # For single ear priors
  single_ear_priors <- table(train_data$Label) / nrow(train_data)
  
  # For joint priors (k1, k2) - modified to be symmetric
  joint_priors <- matrix(0, n_phenotypes, n_phenotypes)
  paired_ears <- train_data %>%
    group_by(SID) %>%
    filter(n() == 2) %>%
    summarise(labels = list(sort(Label)), .groups = "drop")
  
  # Calculate joint prior probabilities from paired ears in training data
  if (nrow(paired_ears) > 0) {
    for (i in 1:nrow(paired_ears)) {
      labels <- unlist(paired_ears$labels[i])
      k1 <- labels[1]
      k2 <- labels[2]
      joint_priors[k1, k2] <- joint_priors[k1, k2] + 1
      if (k1 != k2) {
        joint_priors[k2, k1] <- joint_priors[k2, k1] + 1
      }
    }
    joint_priors <- joint_priors / sum(joint_priors)
  } else {
    # Fallback: use product of marginals if no pairs available
    for (k1 in 1:n_phenotypes) {
      for (k2 in 1:n_phenotypes) {
        joint_priors[k1, k2] <- single_ear_priors[k1] * single_ear_priors[k2]
      }
    }
  }
  
  # Process test data by subject
  grouped_data <- split(test_data, test_data$SID)
  
  for (sid in names(grouped_data)) {
    group <- grouped_data[[sid]]
    
    if (nrow(group) == 1) {
      # SINGLE EAR CASE
      Y_i <- as.numeric(group[, frequencies])
      log_probs <- numeric(n_phenotypes)
      
      for (k in 1:n_phenotypes) {
        mu_k <- theta[k, ]
        Sigma_k <- naive_cov_matrices[[k]]
        Sigma_reg <- Sigma_k + diag(1e-8 * max(diag(Sigma_k)), nrow(Sigma_k))
        
        # Calculate p(y|k) * p(k)
        log_probs[k] <- dmvnorm(Y_i, mu_k, Sigma_reg, log = TRUE) + log(single_ear_priors[k])
      }
      
      # Convert to probabilities (stable way)
      max_log_prob <- max(log_probs)
      if(is.finite(max_log_prob)) {
        probs <- exp(log_probs - max_log_prob)
        probs <- probs / sum(probs)
        
        classification <- which.max(probs)
        pi_ijk <- matrix(probs, 1, n_phenotypes)
        
        all_classifications[[sid]] <- classification
        all_probabilities[[sid]] <- pi_ijk
      }
      
    } else if (nrow(group) == 2 && all(c('L', 'R') %in% group$EAR)) {
      # TWO EAR CASE - Use paired ear covariance with frequency structure
      group <- group[order(group$EAR), ]  # Ensure L then R order
      Y_i <- as.matrix(group[, frequencies])
      y_combined <- c(Y_i[1,], Y_i[2,])  # [L frequencies, R frequencies]
      
      # Calculate joint posterior P(K1=k1, K2=k2 | Y)
      joint_log_probs <- matrix(-Inf, n_phenotypes, n_phenotypes)
      
      for (k1 in 1:n_phenotypes) {
        for (k2 in 1:n_phenotypes) {
          # Get full mean vector and covariance matrix for this phenotype pair
          mu_combined <- c(theta[k1,], theta[k2,])

          # Same phenotype - only need one calculation
          Sigma_combined <- build_cov_matrix(k1, k2, cov_components)

          # Calculate joint likelihood: p(y | k1, k2)
          log_likelihood <- dmvnorm(y_combined, mu_combined, Sigma_combined, log = TRUE)
          
          # Joint prior probability: p(k1, k2)
          log_prior <- log(joint_priors[k1, k2])
          
          # Joint posterior (unnormalized): p(y | k1, k2) * p(k1, k2)
          joint_log_probs[k1, k2] <- log_likelihood + log_prior
        }
      }
    
      # Convert to probabilities (using log-sum-exp trick for numerical stability)
      max_log_prob <- max(joint_log_probs)
      if(is.finite(max_log_prob)) {
        joint_probs <- exp(joint_log_probs - max_log_prob)
        joint_probs <- joint_probs / sum(joint_probs)
        
        # Find classification with highest joint probability
        max_idx <- which(joint_probs == max(joint_probs), arr.ind = TRUE)[1,]
        classifications <- c(max_idx[1], max_idx[2])
        
        # Calculate ear-specific marginal probabilities
        # For left ear: sum over all right ear possibilities
        # For right ear: sum over all left ear possibilities
        pi_ijk <- matrix(0, 2, n_phenotypes)
        for (k in 1:n_phenotypes) {
          pi_ijk[1, k] <- sum(joint_probs[k, ])  # Left ear (row k)
          pi_ijk[2, k] <- sum(joint_probs[, k])  # Right ear (column k)
        }
        
        all_classifications[[sid]] <- classifications
        all_probabilities[[sid]] <- pi_ijk
      }
    }
  }
  
  return(list(
    classifications = all_classifications,
    probabilities = all_probabilities
  ))
}
