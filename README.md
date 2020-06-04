# telo-growth

### WASH Benefits Bangladesh Telomere Length and Growth Associations

This is the repository for the WASH Benefits Bangladesh [NCT01590095](#https://clinicaltrials.gov/ct2/show/NCT01590095) telomere length and growth associations substudy. The primary analysis assessed the associations between telomere length as the exposure and various measures of growth such as length-for-age Z-score as outcomes. Post hoc analyses were conducted to examine the associations in the opposite direction. This repository contains scripts for both the primary and post hoc analyses as well as for generating tables and figures.

### Associated protocols and datasets
The pre-specified analysis plan and the data required for the analysis will be available through the Open Science Framework: [https://osf.io/9snat/](https://osf.io/9snat/).

### WASH Benefits Package and tmleAb
This analysis requires [washb](https://github.com/ben-arnold/washb), a package developed for WASH Benefits analyses. 

For all scripts, you will need to change directory statements to reflect locations of files within your local directory. Changes will also need to be made to directory statements when saving output files such as results, tables, and figures.

### Directory Structure

**`0-config`**: script for configuring data directories and loading libraries

**`analysis`** : analysis scripts

* **`0-base-quartileTMLE_functions`**: quartile TMLE functions
* **`quartileTMLE_telo_growth`**: primary and posthoc analyses

**`figure scripts`** : scripts to make figures

**`figures`** : resulting figures

**`results`** : analysis results

**`tables`** : scripts to make tables and resulting tables

**`tables`** : resulting tables




