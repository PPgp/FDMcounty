# Flow Difference Method

This repository contains code and data for one county to estimate the Bayesian and Deterministic versions of the Flow Difference Method (FDM) for projecting net migration by age.

## Bayesian FDM
**Steps:**

1. Estimate the model via [estimateBayesFDMwacounty.R](estimateBayesFDMwacounty.R).
2. Extract MCMC samples via [extract_rc_samples.R](extract_rc_samples.R)
3. Apply the estimates to a given total net migration via [applyFDM.R](applyFDM.R)

## Deterministic FDM
The estimation and application is shown in [applyFDM.R](applyFDM.R).

# References
H. Sevcikova, J. Raymer J., A. E. Raftery (2024). Forecasting Net Migration By Age: The Flow-Difference Approach. [arXiv:2411.09878](https://arxiv.org/abs/2411.09878).