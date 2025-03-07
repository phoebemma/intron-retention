# Impact of resistance exercise training on splicing efficiency


## Introduction

This study is aimed at investigating the impact of age and resistance
exercise training on splicing efficiency (SE) . It uses SpliceQ
@demelocosta2021 to quantify splicing efficiency in six datasets;

the [COPD study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8261934/)
,

[Volume study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7708234/) ,

[ContraTrain study](https://clinicaltrials.gov/study/NCT03795025) ,

[Alpha and Omega study](https://clinicaltrials.gov/study/NCT04279951) ,

[ReLiEf
study](https://www.inn.no/english/research/research-projects/relief/)

and

[SRP102542](https://pubmed.ncbi.nlm.nih.gov/28273480/)

## Repository Organisation

- `/R`This contains the following sub-folders and scripts with codes
  usedin the data analyses.

  1.  `/Data_extraction` A folder containing six different R scripts
      that shows the extraction of the data from the different RT
      studies . The scripts are named to match their corresponding RT
      intervention studies.

  2.  `/archive` Old scripts used but not wished to be discarded

  3.  `Data_exploration.R` Basic exploration of the intron retention
      data used in this study

  4.  `Figures.R` Codes used to produce the figures/images used in the
      manuscript

  5.  `Preexercise_model.R` Codes for building the baseline model

  6.  `Preexercise_model_exploration.R` Exploration of the baseline
      model outputs

  7.  `RT_effect_model_TrainOme.R` Codes for building the RT model

  8.  `RT_model_exploration.R` Exploration of the RT model

  9.  `Trainome_functions.R` Functions created for use in this study or
      other TrainOme studies

  10. `beta-reg-sim-notes.qmd` Codes used in model simulation

- `/public_data` Contains the files/documents/articles of publicly
  available data used or previously intended to be used in the study

- `/resources` Contains the reference library, csl, and template used
  for the manuscript

- `Relief_sampleIDs.xlsx` Annotation of relief samples/participants.

- `Alpha_Omega_sample_list_transcriptomics.xlsx`
