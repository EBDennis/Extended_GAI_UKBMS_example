# Extended_GAI_UKBMS_example
This repository contains data and code to demonstrate the extension of the generalised abundance index (GAI) approach to incorporate the annual model, with application to [UK Butterfly Monitoring Scheme (UKBMS)](https://ukbms.org/) data. The approach is demonstrated for two butterfly species (Chalk Hill Blue and Gatekeeper).

Please see the associated paper for details: Dennis, E.B., Diana, A., Matechou, E. and Morgan, B.J.T (2023) Efficient statistical inference methods for assessing changes in species' populations using citizen science data. Under submission.

## Instructions

fit_extendedGAI.R will fit the models and produce outputs.

extendedGAI_outputs.R produces plots from the outputs.

## Data

The following data files are provided:

UKBMS_counts_2sp_1976-2022.rds - UKBMS counts for two species  for 1976-2011.

UKBMS_twostageGAI_output_2sp.rds - output from applying the two-stage GAI to UKBMS data for two species.

The [UKBMS](https://ukbms.org/) is organised and funded by Butterfly Conservation, the British Trust for Ornithology (BTO), and the Joint Nature Conservation Committee (JNCC). The UKBMS is indebted to all volunteers who contribute data to the scheme.

## References

Dennis, E.B., Morgan, B.J.T., Freeman, S.N., Brereton, T.M. & Roy, D.B. (2016). A generalized abundance index for seasonal invertebrates. Biometrics, 72, 1305–1314. https://doi.org/10.1111/biom.12506
