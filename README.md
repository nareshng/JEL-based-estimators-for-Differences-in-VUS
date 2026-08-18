# JEL-based-estimators-for-Differences-in-VUS

# Abstract

The volume under the ROC surface (VUS) and its high-dimensional extension, the hypervolume under ROC manifolds (HUM), are key measures for assessing multiclass and high-dimensional classification performance. In this paper, we develop non-parametric inference on differences in the VUS and HUM using the methods: U-statistics and bootstrap based on kernel density and jackknife empirical likelihood (JEL). For developing JEL-based inference, we develop a general framework for JEL for multivariate multi-sample U-statistics and study its properties.  An extensive Monte Carlo simulation is done to evaluate the performance of the finite sample properties of the proposed methods. Finally, we illustrate the practical relevance of our inferential procedures using the Alzheimer’s disease data set.


# JEL-Based Estimators for Differences in VUS

This repository contains the R code used to compare confidence intervals for the difference between two volumes under the ROC surface,

$$\theta = \mathrm{VUS}_1-\mathrm{VUS}_2.$$

The methods considered are:

* Jackknife empirical likelihood (JEL)
* Normal approximation using jackknife variance
* Kernel percentile bootstrap

## Repository files

* `Pareto_Sim_Code.R` — Pareto simulation
* `MOBVE_simulation_code.R` — Marshall–Olkin simulation
* `Data analysis_Code.R` — real-data analysis

## Reproducibility settings

The Pareto simulation uses:

* Master seed: `2026`
* Monte Carlo replications: `1000` per scenario
* Bootstrap replications: `1000`
* Confidence level: `95%`
* JEL evaluation grid: `401` points
* R version: `4.4.2`

Six sample-size configurations are considered:

```text
(10, 20, 30)
(30, 20, 10)
(30, 30, 30)
(60, 70, 40)
(50, 80, 100)
(30, 120, 50)
```

Together with five parameter settings (`P1`–`P5`), these produce 30 Pareto scenarios.

Scenario-specific seeds are generated deterministically:

```text
Data seed =
2026 + 10,000,000 × setting + 100,000 × design + 10 × replication + 1

Bootstrap seed =
2026 + 10,000,000 × setting + 100,000 × design + 10 × replication + 2
```

## Required R packages

```r
install.packages(c("copula", "emplik", "openxlsx"))
```

## Running the analysis

```bash
Rscript Pareto_Sim_Code.R
Rscript MOBVE_simulation_code.R
Rscript "Data analysis_Code.R"
```


The workbook contains parameter settings, true values, sample sizes, scenario summaries, status codes, and replication-level results.

## Note on the P5 kernel results

The zero or near-zero kernel coverage under `P5` is a genuine simulation result, not missing output. All bootstrap intervals were successfully computed. The poor coverage is caused by severe smoothing bias under the extremely heavy-tailed Pareto distributions.

The finalized Marshall–Olkin results will be documented after verification of the corresponding output workbook.
