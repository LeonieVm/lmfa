	
## `lmfa` package
 
`lmfa` allows to explore within-person changes and between person differences in measurement models in (intensive) longitudinal data by means of three-step continuous-time latent Markov factor analysis (LMFA).
 
### Installation
You can download the development version from GitHub as follow:

```javascript
install.packages("devtools"); library(devtools)

devtools::install_github("leonievm/lmfa")
```

 
### Usage
 
After successful installation, you can perform LMFA by means of the three-step estimation. The package consists of the three main functions. For details about the function arguments, you may call the function documentations with ?'functionname' 

1. The step 1 function estimates the state-specific factor analysis models by means of an expectation maximization algorithm (with or without model selection):
```javascript
step1(data,
      indicators,
      n_state = NULL,
      n_fact = NULL, 
      modelselection = FALSE, 
      n_state_range = NULL, 
      n_fact_range = NULL,
      n_starts = 25,
      n_initial_ite = 10,
      n_m_step = 10,
      em_tolerance = 1e-8, 
      m_step_tolerance = 1e-3, 
      max_iterations = 1000,
      n_mclust = 5)
```

2. The step 2 function obtains the posterior state-membership probabilities and the modal state assignments and calculates the classification error: 
```javascript
step2(data, model)
```

3. The step 3 function estimates the transitions between the states (conditional on covariates) by means of a continuous-time latent Markov model:
```javascript
step3(data,
      identifier,
      n_state,
      postprobs,
      timeintervals = NULL,
      initialCovariates = NULL,
      transitionCovariates = NULL,
      n_starts = 25,
      n_initial_ite = 10,
      method = "BFGS",
      max_iterations = 10000,
      tolerance = 1e-10,
      scaling = "proxi")
```

### Contribution

If you have any suggestions or if you found any bugs, please feel free to contact me via email.
