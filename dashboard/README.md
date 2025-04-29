# EarCorrelator

## Advanced Quadratic Discriminant Analysis for Correlated Audiometric Data

This repository implements enhanced Quadratic Discriminant Analysis (QDA) algorithms that account for complex correlation structures in audiometric data. The methods are particularly designed for paired-organ data where standard independence assumptions are violated.

### Features

- **GEE-based QDA**: Implements a parsimonious correlation structure with three key parameters:
  - Within-ear correlations between frequencies
  - General between-ear correlations
  - Between-ear correlations at the same frequency

- **Mixed Model QDA**: Directly models the hierarchical structure through nested random effects:
  - Person-level random effects
  - Ear-specific random effects within individuals
  - Frequency-specific random effects within phenotypes
      
- **Missing Data Handling**: Accommodates participants with measurements from one or both ears

## Ear Type Classification Dashboard

This interactive dashboard provides real-time classification of ear types using QDA methods.

### Features
- Upload CSV data or enter values manually
- Run both Naive QDA and Improved QDA classification
- Visualize audiograms and classification probabilities
- Compare results between different classification methods

### Setup Instructions
1. Ensure R and required packages are installed
2. Run `app.R` to start the dashboard
3. Pre-trained models will be generated on first run

### Deployment
The dashboard is deployed at: https://bellayqian.shinyapps.io/EarCorrelator/

### Data Format

The package expects data in long format with:
- One row per ear per subject
- Columns for each frequency measurement
- Columns identifying the subject, ear, and phenotype

Example:
```
subject_id ear_id phenotype f250 f500 f1000 f2000 f3000 f4000 f6000 f8000
1          1      2         10   15   20    25    35    45    50    55
1          2      2         12   18   22    28    38    48    52    58
2          1      3         15   25   35    45    60    70    75    80
...
```

### Applications

Beyond audiometric data, this methodology can be extended to:

1. Multi-center clinical trials (modeling correlations between treatment arms within centers)
2. Healthcare systems research with nested data structures
3. Medical research involving paired organs (eyes, kidneys, lungs)

### Citation

If you use this software in your research, please cite:

```
Qian, Y., Chen, Y., & Wang, M. (2025). Advanced quadratic discriminant 
analysis algorithms for correlated audiometric data. Manuscript in preparation.
```

### License

This project is licensed under the MIT License - see the LICENSE file for details.

### Acknowledgments

This research was conducted at the Harvard T.H. Chan School of Public Health Department of Biostatistics.
