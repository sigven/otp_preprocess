## Scripts to parse JSON files from Open Targets Platform

This repository contains scripts (R) to establish global datasets from the Open Targets Platform, including:

-   gene-disease associations
-   drug-target associations
-   drug-disease indications
-   other functional target information (tractability data, function description, cancer hallmarks information etc.)

### Overview

-   Code to download/rsync raw JSON files from the Open Targets Platform and preprocess them is available in [code](code)
-   Raw JSON files are available under [data](data)
-   Processed data (RDS files, for use in R) is available in [output](output)
