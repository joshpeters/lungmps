---
layout: default
title: Processing
nav_order: 2
description: "Lung MPS"
permalink: /processing
---

## Collating and processing data

All data was downloaded using the links in `study_metadata.csv`. These raw files were downloaded to `data/raw_data` and were loaded and processed using the scripts inside `src/dataset_preparation` in R.

We attempted to harmonize several metadata variables, including `cell_barcode, study, sample_id, batch, donor, condition, disease_classification, celltype`.

This generated R objects and 10X-formatted directories for every study in this work. The raw R objects are available on Zenodo.

We omitted samples from several studies, which can be seen within the scripts or in Supplementary Table 1.



---
[View the code on Github](https://github.com/joshpeters/lungmps){: .btn .fs-2}
