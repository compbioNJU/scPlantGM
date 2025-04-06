# scPlantGM
`scPlantGM` is a hierarchical and gene module-based annotation framework optimized for plant single-cell data. With high accuracy, biological interpretability and robustness to data imbalance, `scPlantGM` has capability to conduct mutiple functions, for example, **automatic cell type annotation**, marker reference generation and cross-species conservatism analysis etc.

# Overview
<img width="779" alt="17597cda3b3a538529cfcdd23c9e301" src="https://github.com/user-attachments/assets/d1efd773-230b-4f19-8b43-c1d693b3c084" />

# Reference
Before usage, please remember to cite our work. Comprehensive Annotation of Plant Cell Atlases Across Species Based on Gene Modules.(will be a link)

# Installation
You can install `scPlantGM` through devtools as follows:

```r
library(devtools)
devtools::install_github("XuanKunXK/scPlantGM")
```

# Usage
## Prerequisite
Your data should be processed by classic Seurat pipnline until finishing clustering. Data preprocessing pipnline refers to [Seurat](https://satijalab.org/seurat/)

## Mode
In summary, `scPlantGM` owns three modes based on different reference sources. User can flexibly choose which mode is more suitable in `custom` argument.
- `None`: user won't provide any custom reference but use built-in reference (built-in reference includes arabidopisis, maize, rice)
- `Semi`: user will integrate their own reference and built-in reference for comprehensive usage.
- `All`: user will only use their provided reference for prediction.
### None mode
#### requirements
### Semi mode
#### requirements
### All mode
#### requirements

## Prediction

## Result

## Advanced usage

# Q&A

# Contact us
If you have any question, suggestion or bug found in the method. Feel free to contact us! Email: XuanKunXK@gmail.com
