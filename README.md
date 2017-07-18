# flowLearn
flowLearn: Fast and precise identification and quality checking of cell populations in flow cytometry

## Abstract

### Motivation

Identification of cell populations in flow cytometry is a critical part of
analysis and lays the groundwork for both clinical diagnostics and research discovery. The current
paradigm of manual analysis is time consuming and subjective. If the goal is to match manual
analysis, supervised tools provide the best performance, however they require fine parameterization
to obtain the best results, Hence, there is a strong need for methods that are fast to setup,
accurate and interpretable at the same time.

### Results

flowLearn is a semi-supervised
approach for the quality-checked identification of cell populations. Using as few as one manually
gated sample, through density alignments it is able to predict gates on other samples with high
accuracy and speed. On two state-of-the-art data sets, our tool achieves `median F_1`-measures
exceeding `F_1 > 0.99` for `31%`, and `F_1 > 0.90` for `80%` of all analyzed populations.
Furthermore, users can directly interpret and adjust automated gates on new sample files to
iteratively improve the initial training.

## Using this software

- FlowLearn is very well documented and provides a quick start vignette with some example data.

- Clone the repository and build the `R` package using `devtools::install('package', build_vignettes = T)`

- Load the package using `library(flowLearn)` and view the quick-start vignette using `vignette('flowLearnVignette')`.
