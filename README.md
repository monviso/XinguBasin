# Forest structure types in the southern Amazon Basin are influenced by anthropogenic edges and physiography

This repository contains the R and SLURM scripts used to reproduce the analyses presented in **“Forest structure types in the southern Amazon Basin are influenced by anthropogenic edges and physiography.”** The workflow identifies GEDI-derived forest vertical architecture classes, evaluates their robustness, independently characterizes them using optical and airborne-lidar data, and models their environmental and landscape associations.

If you use these scripts or adapt the workflow for another study, please cite the associated article and the archived dataset:

> **Data and analysis outputs:** [https://doi.org/10.6084/m9.figshare.34023453](https://doi.org/10.6084/m9.figshare.34023453)

The complete citation for the article should be added here once it is available.

## 1. Data availability

The input files required to reproduce the published analyses, together with intermediate outputs and final figures, are archived on Figshare:

[https://doi.org/10.6084/m9.figshare.34023453](https://doi.org/10.6084/m9.figshare.34023453)

Download the required files before running the scripts. In particular:

- the GEDI sampling and classification inputs are required for Steps 01–09;
- the overlapping airborne laser scanning (ALS) data are required for Steps 07–08;
- `Xingu_GAM_prepared.rds` contains the assembled modelling dataset required to run the GAM fitting and downstream analyses without repeating predictor extraction;
- intermediate result tables and figures are provided for users who wish to inspect the published outputs without rerunning the complete workflow.

The ALS data should be used in accordance with the terms and citation requirements of their original data providers, as documented in the archive.

## 2. Computing environments

Steps 01–09 were run locally on a Dell Precision 3680 with the following specifications:

- Intel Core i7-14700K processor at 3.4 GHz;
- 20 physical cores;
- 28 logical processors.

Steps 13–14 and Step 22 were run on the Northern Arizona University Monsoon high-performance computing system because of their computational and memory requirements:

[https://in.nau.edu/arc/details/](https://in.nau.edu/arc/details/)

The supplied `.sbatch` files use the SLURM workload manager. Resource requests, module loading, conda environments and filesystem paths may need to be modified for other HPC systems.

Steps 15–17 and 23–24 can be run on a suitable local workstation after the required model outputs have been generated or downloaded.

## 3. Workflow overview

| Steps | Main purpose | Recommended environment |
|---|---|---|
| 01–09 | GEDI filtering and sampling, NMF and k-means classification, robustness analyses, radiometric characterization, and independent ALS characterization | Local workstation |
| 10–12 | Predictor preparation and assembly; not distributed in this repository | See Section 5 |
| 13–14 | Generalized additive model fitting using `Xingu_GAM_prepared.rds` | HPC with SLURM |
| 15–17 | Post-fitting model summaries, diagnostics and inferential outputs | Local workstation |
| 22 | Computationally intensive extraction and plotting of model partial effects | HPC with SLURM |
| 23–24 | Final local summaries and manuscript figures | Local workstation |

Missing script numbers correspond to ancillary or development scripts whose functions were incorporated into the current workflow, or to the predictor-generation stage described below. The numbered scripts distributed here constitute the final analysis sequence used for the revised study.

## 4. Steps 01–09: classification and independent characterization

Scripts 01–09 cover:

1. preparation and prescreening of the GEDI discovery sample;
2. non-negative matrix factorization (NMF) rank selection;
3. within-component k-means subdivision and final class assignment;
4. classification stability and sensitivity analyses;
5. extraction and analysis of optical radiometric information;
6. identification of spatial overlap between GEDI footprints and ALS coverage;
7. extraction of ALS footprint profiles and structural metrics; and
8. independent characterization of GEDI architecture classes using ALS observations.

The classification labels are ordered by increasing class-average RH98: Class 1 represents the shortest architecture and Class *N* the tallest architecture for a solution containing *N* classes.

To reproduce Steps 07–08, download the ALS data overlapping GEDI observations in the Xingu Basin from the Figshare archive and update the scripts with the corresponding local directory. ALS measurements are used as an independent structural characterization rather than as a direct replication of GEDI waveform classification.

## 5. Predictor preparation

The original predictor-extraction scripts are not included in this repository. The predictors were derived from a new Amazon-wide raster product describing annual land-cover composition, forest-patch configuration and physiographic conditions. That product and its production framework are currently under review.

To make the statistical analyses reproducible without redistributing the unpublished predictor-generation workflow, the Figshare archive provides:

```text
Xingu_GAM_prepared.rds
```

This object contains the response variables, modelling partitions, sampling weights, predictor histories and ancillary information required by the GAM-fitting scripts. Consequently, users can reproduce Steps 13 onward without rerunning predictor calculation or point extraction.

## 6. Steps 13–14: GAM fitting

Steps 13–14 fit the one-versus-rest generalized additive models used to evaluate environmental, landscape and spatial associations with each forest architecture class.

On a SLURM-based HPC system, edit the supplied `.sbatch` file to specify:

- the location of `Xingu_GAM_prepared.rds`;
- the output directory;
- the R environment and required package library;
- the requested memory, runtime and array configuration; and
- any cluster-specific account or partition settings.

Submit the fitting job using the appropriate `sbatch` command. Model fits and task-specific logs are written to the configured output directories. Review the log files and confirm that every expected model task completed before running downstream summaries.

## 7. Downstream summaries and figures

Steps 15–17 summarize model performance, explained deviance, predictor contributions and statistical significance. Step 22 generates the complete set of partial-effect outputs and may require substantial memory. Steps 23–24 produce downstream summaries and publication figures that can generally be regenerated locally from the archived model outputs.

Large numerical tables are provided through the Figshare archive rather than embedded in the manuscript Supporting Material. This includes detailed robustness results, pairwise comparisons and complete model-effect tables.

## 8. Configuration and reproducibility notes

Before running any script:

1. inspect its configuration block;
2. replace Windows and Monsoon paths with paths appropriate to your system;
3. create the required output directories;
4. verify that all expected input files exist;
5. check overwrite and resume options before rerunning completed stages; and
6. retain the random seeds supplied in the scripts when exact reproducibility is required.

Do not assume that the default absolute paths will exist on another computer. Scripts that process large datasets should first be tested on a small subset before starting the complete analysis.

## 10. Citation and reuse

If you reuse the data, scripts, classification workflow or derived outputs, please cite:

1. the associated research article;
2. the Figshare archive: [https://doi.org/10.6084/m9.figshare.34023453](https://doi.org/10.6084/m9.figshare.34023453); and
3. the original GEDI, HLS, MapBiomas, TMF and ALS datasets applicable to the reused analysis.

Questions, bug reports and reproducibility issues should be submitted through the repository issue tracker, including the script number, software environment, relevant log output and the command used to run the analysis.
