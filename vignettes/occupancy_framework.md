---

title: "Bayesian Occupancy Modelling Framework for Breeding Bird Indicators"
author: "Your Name"
output: rmarkdown::html_vignette
bibliography: references.bib
vignette: >
%\VignetteIndexEntry{Bayesian Occupancy Modelling Framework}
%\VignetteEngine{knitr::rmarkdown}
\usepackage[utf8]{inputenc}
---------------------------

# Introduction

This vignette describes a Bayesian occupancy modelling framework for breeding bird monitoring data. The framework uses the `spOccupancy` package to fit single-species and multi-species occupancy models with optional spatial random effects.

The aim is to estimate breeding bird indicators while accounting for imperfect detection, spatial dependence, and variation among species or ecological guilds.

# Workflow

```{mermaid}
flowchart TD
  A[Raw breeding bird survey data] --> B[Validate required columns]
  B --> C[Create detection/non-detection records]
  C --> D[Build species x site x year array]
  D --> E[Standardize survey effort]
  E --> F[Prepare occupancy covariates]
  E --> G[Prepare detection covariates]
  F --> H[Scale numeric covariates]
  G --> H
  H --> I[Assign species to guilds or species groups]
  I --> J{Model configuration}
  J --> K[PGOcc: single-species non-spatial]
  J --> L[spPGOcc: single-species spatial]
  J --> M[msPGOcc: multi-species non-spatial]
  J --> N[spMsPGOcc: multi-species spatial]
  K --> O[Posterior summaries]
  L --> O
  M --> O
  N --> O
  O --> P[Occupancy psi]
  O --> Q[Latent occupancy z]
  O --> R[Species richness]
  O --> S[Spatial effects w]
  O --> T[Diagnostics: PPC, WAIC, R-hat]
```

# Data structure

The input data are expected to contain:

* site identifier;
* species name;
* survey year;
* count or detection information;
* longitude and latitude;
* occupancy covariates;
* detection covariates.

Counts are converted to binary detection/non-detection values:

[
y_{ijk} =
\begin{cases}
1, & \text{if species } i \text{ is detected at site } j \text{ in year } k \
0, & \text{otherwise}
\end{cases}
]

where:

* (i = 1,\dots,N) indexes species,
* (j = 1,\dots,J) indexes sites,
* (k = 1,\dots,K) indexes survey years or sampling occasions.

# Survey effort standardization

A key feature of this framework is that survey effort is defined before filtering to species or guilds.

The model constructs a complete array:

[
\mathbf{Y} = {y_{ijk}}
]

with dimensions:

[
N \times J \times K
]

For each surveyed site-year combination, species that were not recorded are assigned zero rather than missing. This means that non-detections are defined consistently across single-species and multi-species analyses.

This step is important because it ensures that differences among species, guilds, or ecosystems are not artifacts of inconsistent survey effort.

# Hierarchical occupancy model

Occupancy models separate the ecological process from the observation process.

The ecological process is:

[
z_{ij} \sim \text{Bernoulli}(\psi_{ij})
]

where:

* (z_{ij}) is the latent true occupancy state;
* (\psi_{ij}) is the probability that species (i) occupies site (j).

The observation process is:

[
y_{ijk} \sim \text{Bernoulli}(z_{ij}p_{ijk})
]

where:

* (p_{ijk}) is the detection probability;
* (y_{ijk}=0) may mean either true absence or non-detection.

# Occupancy submodel

For a non-spatial single-species model:

[
\text{logit}(\psi_j) =
\beta_0 + \sum_{q=1}^{Q} x_{jq}\beta_q
]

For a multi-species model:

[
\text{logit}(\psi_{ij}) =
\beta_{0i} + \sum_{q=1}^{Q} x_{jq}\beta_{qi}
]

where:

* (x_{jq}) is occupancy covariate (q) at site (j);
* (\beta_{qi}) is the species-specific occupancy effect.

# Detection submodel

Detection probability is modelled as:

[
\text{logit}(p_{ijk}) =
\alpha_{0i} + \sum_{r=1}^{R} v_{jkr}\alpha_{ri}
]

where:

* (v_{jkr}) is detection covariate (r);
* (\alpha_{ri}) is the species-specific detection effect.

Detection covariates may include observer effort, survey year, weather, or other observation-level variables.

# Community-level hierarchical priors

For multi-species models fitted with `msPGOcc()` or `spMsPGOcc()`, species-specific coefficients are treated as exchangeable draws from community-level distributions.

For occupancy coefficients:

[
\beta_{qi} \sim \text{Normal}(\mu_{\beta_q}, \tau_{\beta_q}^{2})
]

For detection coefficients:

[
\alpha_{ri} \sim \text{Normal}(\mu_{\alpha_r}, \tau_{\alpha_r}^{2})
]

where:

* (\mu_{\beta_q}) is the community-level mean occupancy effect for covariate (q);
* (\tau_{\beta_q}^{2}) is the among-species variance in occupancy responses;
* (\mu_{\alpha_r}) is the community-level mean detection effect for covariate (r);
* (\tau_{\alpha_r}^{2}) is the among-species variance in detection responses.

The community-level priors can be written as:

[
\mu_{\beta_q} \sim \text{Normal}(0, \sigma_{\mu_\beta}^{2})
]

[
\tau_{\beta_q}^{2} \sim \text{Inverse-Gamma}(a_\beta, b_\beta)
]

[
\mu_{\alpha_r} \sim \text{Normal}(0, \sigma_{\mu_\alpha}^{2})
]

[
\tau_{\alpha_r}^{2} \sim \text{Inverse-Gamma}(a_\alpha, b_\alpha)
]

This structure allows partial pooling across species. Rare species borrow information from the community, while common species are informed more strongly by their own detection histories.

# Spatial occupancy model

Spatial models include a spatial random effect in the occupancy process:

[
\text{logit}(\psi_{ij}) =
\beta_{0i} + \sum_{q=1}^{Q} x_{jq}\beta_{qi} + w_{ij}
]

where (w_{ij}) captures spatial residual variation for species (i) at site (j).

The spatial random effects are modelled using a Gaussian process:

[
\mathbf{w}_i \sim \text{MVN}(\mathbf{0}, \mathbf{C}_i)
]

with covariance:

[
C_i(j,j') =
\sigma_i^2 \exp(-\phi_i d_{jj'})
]

where:

* (d_{jj'}) is the distance between sites (j) and (j');
* (\sigma_i^2) is the spatial variance;
* (\phi_i) is the spatial decay parameter.

The effective range is:

[
\text{Effective range}_i =
\frac{3}{\phi_i}
]

Spatial models in this framework can use nearest-neighbor Gaussian processes, which are computationally efficient for large spatial datasets.

# Model routing

The function routes models according to two user choices:

| `single_species` | `fit_spatial_model` | spOccupancy function |
| ---------------- | ------------------- | -------------------- |
| `TRUE`           | `FALSE`             | `PGOcc()`            |
| `TRUE`           | `TRUE`              | `spPGOcc()`          |
| `FALSE`          | `FALSE`             | `msPGOcc()`          |
| `FALSE`          | `TRUE`              | `spMsPGOcc()`        |

# Species grouping

The framework supports three grouping strategies:

1. **Single-species models**: each species is fitted independently.
2. **User-defined guilds**: species are grouped using a supplied guild table.
3. **Automatic grouping**: species are split into groups when the number of species exceeds `maxSpeciesInGrps`.

# Derived ecological quantities

## Occupancy probability

[
\psi_{ij}
]

is the posterior probability that species (i) occupies site (j).

## Latent occupancy state

[
z_{ij}
]

is the posterior distribution of true occurrence.

## Species richness

For multi-species models, richness at site (j) is estimated as:

[
S_j = \sum_{i=1}^{N} z_{ij}
]

Posterior summaries of (S_j) provide expected richness and credible intervals.

## Spatial random effects

[
w_{ij}
]

represents spatial structure not explained by covariates. Mapping (w_{ij}) can reveal residual spatial clustering.

# Bayesian inference

Models are fitted using MCMC. Posterior summaries include:

* posterior mean;
* posterior standard deviation;
* 95% credible interval;
* posterior probability that a coefficient is positive.

For a coefficient (\theta), the posterior probability of a positive effect is:

[
Pr(\theta > 0 \mid \mathbf{Y})
]

This quantity is useful for ecological interpretation of covariate effects.

# Model diagnostics

## Posterior predictive checks

Posterior predictive checks compare observed data to replicated data generated from the fitted model.

Using a discrepancy statistic (T(\cdot)), the Bayesian p-value is:

[
p_B =
Pr{T(\mathbf{Y}^{rep}) > T(\mathbf{Y}) \mid \mathbf{Y}}
]

Values close to 0.5 indicate adequate fit. Values close to 0 or 1 may indicate lack of fit.

In `spOccupancy`, posterior predictive checks can be computed with:

```r
ppcOcc(model, fit.stat = "freeman-tukey", group = 1)
```

## WAIC

Model comparison can be performed using WAIC:

```r
waicOcc(model)
```

Lower WAIC values indicate better expected out-of-sample predictive performance, conditional on the set of models being compared.

## MCMC convergence

Convergence is assessed using the Gelman-Rubin statistic:

[
\hat{R}
]

Values close to 1 indicate acceptable convergence.

# Output structure

The framework returns a nested R object containing:

* ecosystem name;
* model type;
* occupancy covariates;
* detection covariates;
* results by guild or species;
* combined summary tables.

The combined output may include:

* `psi`: occupancy probability summaries;
* `z`: latent occupancy summaries;
* `richness`: species richness estimates;
* `spatial_random_effects`: posterior summaries of (w);
* `spatial_params`: (\sigma^2), (\phi), and effective range;
* `occ_effects`: occupancy covariate effects;
* `det_effects`: detection covariate effects;
* `ppc`: posterior predictive checks.

# Example usage

```r
results <- fit_ecosystem_spoccupancy(
  ecosystemBreedingBird = bird_data,
  ecosystemType = "Forest",
  occ_covariates = c("elevation", "forest_cover"),
  det_covariates = c("observer_effort"),
  species_groupings = guild_table,
  single_species = FALSE,
  fit_spatial_model = TRUE,
  n.batch = 1000,
  batch.length = 25,
  n.burn = 5000,
  n.thin = 5,
  n.chains = 3
)
```

Extract occupancy estimates:

```r
results$combined$psi
```

Extract richness estimates:

```r
results$combined$richness
```

Extract covariate effects:

```r
results$combined$occ_effects
results$combined$det_effects
```

# Interpretation

Positive occupancy coefficients indicate that occupancy probability increases with the covariate. Positive detection coefficients indicate that detection probability increases with the detection covariate.

Spatial effects should be interpreted as residual spatial structure after accounting for measured covariates. Large positive spatial effects indicate sites with higher occupancy than expected from the covariates alone.

Species richness estimates summarize community-level occupancy patterns and can be used as breeding bird biodiversity indicators.

# References

Datta, A., Banerjee, S., Finley, A. O., and Gelfand, A. E. (2016). Hierarchical nearest-neighbor Gaussian process models for large geostatistical datasets. *Journal of the American Statistical Association*, 111, 800–812.

Dorazio, R. M., and Royle, J. A. (2005). Estimating size and composition of biological communities by modeling the occurrence of species. *Journal of the American Statistical Association*, 100, 389–398.

Doser, J. W., Finley, A. O., Kéry, M., and Zipkin, E. F. (2022). spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. *Methods in Ecology and Evolution*, 13, 1670–1678.

Finley, A. O., Datta, A., Banerjee, S., and Gelfand, A. E. (2019). spNNGP: Spatial modeling using nearest neighbor Gaussian processes. *Journal of Statistical Software*, 91, 1–30.

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, J. A., and Langtimm, C. A. (2002). Estimating site occupancy rates when detection probabilities are less than one. *Ecology*, 83, 2248–2255.

Polson, N. G., Scott, J. G., and Windle, J. (2013). Bayesian inference for logistic models using Pólya-Gamma latent variables. *Journal of the American Statistical Association*, 108, 1339–1349.
