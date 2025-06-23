# Survival Analysis Simulation: Exponential Distribution with Censoring

## Overview

This repository contains R code for simulating survival data from a hypothetical clinical trial with two treatment arms. The simulation demonstrates:
- Exponentially distributed survival times
- Non-informative right censoring
- Kaplan-Meier survival curves
- Cox proportional hazards modeling

## Key Results from Simulation

Based on the simulation with seed 123:

### Patient Distribution
- **Control Group**: 265 patients
- **Treatment Group**: 235 patients
- **Total**: 500 patients

### Event Summary
- **Total Events**: 310 (62%)
- **Censored**: 190 (38%)
- **Control Events**: 183/265 (69.1%)
- **Treatment Events**: 127/235 (54.0%)

### Primary Outcome
- **Hazard Ratio (Treatment vs Control)**: 0.573 (95% CI: 0.456 - 0.721)
- **p-value**: < 0.001
- **Interpretation**: The treatment reduces the hazard of death by approximately 43%

### Median Survival Times
| Group | Observed Median | Theoretical Median |
|-------|-----------------|-------------------|
| Control | 6.91 | 6.93 |
| Treatment | 12.21 | 9.90 |

## Repository Structure

```
survival-analysis-simulation/
│
├── README.md               # This file
├── simulation_code.R       # Main simulation script
├── output/
│   ├── simulation_output.txt    # Console output
│   └── plots/                   # Generated plots (if saved)
└── docs/
    └── simulation_guide.md      # Detailed exercise instructions
```

## Requirements

### R Packages
- `survival` - Core survival analysis functions
- `survminer` - Enhanced visualization of survival curves
- `ggplot2` - Required by survminer
- `ggpubr` - Required by survminer

### Installation
```r
install.packages(c("survival", "survminer"))
```

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| n | 500 | Total number of patients |
| λ_control | 0.1 | Hazard rate for control group |
| λ_treatment | 0.07 | Hazard rate for treatment group |
| λ_censor | 0.05 | Censoring rate |
| Proportion treated | 0.5 | Randomization ratio |

## Key Features

1. **Random Assignment**: Uses `rbinom()` for treatment allocation
2. **Data Generation**: Exponential distributions via `rexp()`
3. **Censoring Mechanism**: Independent exponential censoring
4. **Statistical Analysis**:
   - Kaplan-Meier curves with confidence intervals
   - Log-rank test
   - Cox proportional hazards model
   - Hazard ratio estimation

## Bonus Analysis: Age Covariate

The simulation includes an optional age covariate:
- Age ~ Normal(60, 10)
- Centered for interpretation
- Shows minimal impact on treatment effect (HR remains 0.573)

## Visualizations Generated

1. **Kaplan-Meier Curves** (Base R)
2. **Enhanced KM Plot** (ggsurvplot) with:
   - Confidence intervals
   - Risk table
   - p-value
3. **Age Distribution Boxplot**
4. **Schoenfeld Residuals Plot** (PH assumption check)

## Reproducibility

The simulation uses `set.seed(123)` for complete reproducibility. Running the code will generate identical results each time.

## Educational Value

This simulation demonstrates:
- Proper survival data generation
- Realistic censoring patterns
- Agreement between theoretical and observed values
- Standard survival analysis workflow
- Model diagnostics and assumptions

## Citation

If you use this code for educational purposes, please reference:
```
Survival Analysis Simulation: Exponential Distribution with Censoring
GitHub: [your-username]/survival-analysis-simulation
```

## License

This code is provided for educational purposes. Feel free to use and modify for teaching or learning survival analysis concepts.
