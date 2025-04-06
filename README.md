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

## Data preprocession
User must preprocess Seurat data before prediction as follows:
### Prepare `seuorder`
One necessary and important object--`seuorder`: a vector shows the order of Suerat objects in the list. Strongly recommend using 'Sample' to present the order. Or user can just simply use `c(1,2,3...)` for a quick start.
### Preprocess data
```r
query <- process_obj(query, seuorder, 'query')

reference <- process_obj(reference, seuorder, 'reference')
```
### Annotate reference
If you want to use your own reference, please don't forget to add annotation infomation into reference Seurat objects manually.

```r
reference[[1]]@meta.data$scPlantGM.refanno = ref_annotation#ref_annotation is a vector carrying reference annotation
```

## Mode
In summary, `scPlantGM` owns three modes based on different reference sources. User can flexibly choose which mode is more suitable in `custom` argument.
- `None`: user won't provide any custom reference but use built-in reference (built-in reference includes Arabidopis, Maize, Rice)
- `Semi`: user will integrate their own reference and built-in reference for comprehensive usage.
- `All`: user will only use their provided reference for prediction.

## Requirements
### `None` mode
- `query`: a (list of) query Seurat object that you want to conduct prediction.
- `species`: species of scRNA-seq data. You can only choose within `Arabidopis`, `Maize`, and `Rice` in this mode.
- `organ`: organ of species. Built-in organs: Arabidopis(`Root`,`Leaf`,`Pollen`,`Inflorescence`), Maize(`Root`,`Leaf`), Rice(`Root`,`Leaf`).
- `custom` = `None`.
### `Semi` mode
- `query`: a (list of) query Seurat object that you want to conduct prediction.
- `reference`: a (list of) Seurat object as reference.
- `species`: species of scRNA-seq data. You can only choose within `Arabidopis`, `Maize`, and `Rice` in this mode.
- `organ`: organ of species. Built-in organs: Arabidopis(`Root`,`Leaf`,`Pollen`,`Inflorescence`), Maize(`Root`,`Leaf`), Rice(`Root`,`Leaf`).
- `custom` = `Semi`.
### `All` mode
- `query`: a (list of) query Seurat object that you want to conduct prediction.
- `reference`: a (list of) Seurat object with annotation as reference.
- `custom` = `All`.

## Prediction
Run the method:

```r
result <- scPlantGM(query, reference, species, organ, custom)
```

## Result
The prediction result will be stored as a data.frame as follows:

<img width="289" alt="38066e856f5742ea8bff53e7115483d" src="https://github.com/user-attachments/assets/3fb76cb9-9107-4deb-b371-eee04b3f9d46" />


## Advanced usage
- layer_info: a format-fixed data.frame includes layers information (will be ignored if custom='None'). A example shown in data/layer_info_example.csv
- layer: layer of prediction: must be an int (0,1,2...). In `None` and `Semi` mode, user can choose one of built-in layer: Arabidopis(`1`,`2`,`3`), Maize(`1`,`2`), Rice(`1`,`2`)
- p_thres: a threshold for determining what level of module purity is acceptable. The higher the value, the purer the module are.
- m_thres: a threshold for determining how many modules can be accepted. The higher the value, the more modules can be accepted.
- x_thres: a threshold for determining whether the cell type should be ‘unknown’.
- cores: numbers of CPU for running.

# Q&A

# Contact us
If you have any question, suggestion or bug found in the method. Feel free to contact us! Email: XuanKunXK@gmail.com
