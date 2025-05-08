# scPlantGM
`scPlantGM` is a hierarchical and gene module-based annotation framework optimized for plant single-cell data. With high accuracy, biological interpretability and robustness to data imbalance, `scPlantGM` demonstrates strong capabilities to conduct mutiple functions, including **automatic cell type annotation**, **marker reference generation** and **cross-species conservatism analysis** etc.

# Overview
<img width="779" alt="17597cda3b3a538529cfcdd23c9e301" src="https://github.com/user-attachments/assets/d1efd773-230b-4f19-8b43-c1d693b3c084" />

# Installation
You can install `scPlantGM` through devtools as follows:

```r
library(devtools)
devtools::install_github("compbioNJU/scPlantGM")
```



# Usage
In summary, `scPlantGM` owns three modes based on different reference sources. User can flexibly choose which mode is more suitable in `custom` argument.
- `None`: users won't provide any custom reference but use built-in reference (built-in reference includes Arabidopis, Maize, Rice)
- `Semi`: users will integrate their own reference and built-in reference for comprehensive usage.
- `All`: users will only use their provided reference for prediction.

We will introduce the three modes, seprately.

<details>
<summary>1. Annotate using built-in reference</summary>

In 'None' mode, users can directly input the query list for annotation. Please note that the implementation of scPlantGM in "None" mode is currently limited to some tissues of Arabidopis(`Root`, `Leaf`, `Pollen`, `Inflorescence`), Maize(`Root`, `Leaf`), Rice(`Root`, `Leaf`). The specific tissues included can also be viewed through the following command.
```r
data(Species_Tissue, package = "scPlantGM")
```

## Load data
We load `Inflorescence_query` as query list. Seurat object in query list needs a column named `seurat_clusters`, representing the clustering result, which can be obtained by the `FindCluster` function during [Seurat](https://satijalab.org/seurat/) pipnline.

```r
data("Inflorescence_query", package = "scPlantGM")
```


## Preprocess
Anothor necessary and important object is `seuorder`,  which is a vector shows the order of Suerat objects in the list. Strongly recommend using 'Sample' to present the order. Or user can just simply use `c(1,2,3...)` for a quick start.
```r
query <- process_obj(query, seuorder = c("SRX13437605"), type = 'query')
```


## Prediction
Run scPlantGM for cell type prediction.
```r
result <- scPlantGM(query, species="Arabidopis", organ="Inflorescence", custom="None", layer = 2)
```
- `query`: a (list of) query Seurat object that you want to conduct prediction.
- `species`: species of scRNA-seq data. You can only choose within `Arabidopis`, `Maize`, and `Rice` in this mode.
- `organ`: organ of species. Built-in organs: Arabidopis(`Root`, `Leaf`, `Pollen`, `Inflorescence`), Maize(`Root`, `Leaf`), Rice(`Root`, `Leaf`).
- `custom` = `None`.
- `layer`: layer of prediction: must be an int (0,1,2...). In `None` and `Semi` mode, user can choose one of built-in layer: Arabidopis(`1`, `2`, `3`), Maize(`1`, `2`), Rice(`1`, `2`).

The prediction result will be stored as a data.frame as follows:
```r
> head(result)
#                              Cell       Cluster            prediction1      prediction2             prediction
#1 SRX13437605@@_AAACCCAAGCCGCACT-1 SRX13437605:0          Ground tissue           Cortex                 Cortex
#2 SRX13437605@@_AAACCCAAGCTTGTTG-1 SRX13437605:1          Ground tissue           Cortex                 Cortex
#3 SRX13437605@@_AAACCCAAGGAACTAT-1          <NA>                   <NA>             <NA>                   <NA>
#4 SRX13437605@@_AAACCCAAGTATGGCG-1 SRX13437605:8        Vascular tissue Vascular cambium       Vascular cambium
#5 SRX13437605@@_AAACCCAAGTCTAGAA-1 SRX13437605:4  Shoot apical meristem  Flower meristem        Flower meristem
#6 SRX13437605@@_AAACCCACAAAGGAGA-1 SRX13437605:3 Shoot system epidermis                / Shoot system epidermis
```

</details>



<details>
<summary>2. Annotate using user-provided and built-in reference</summary>
In `Semi` mode, users need to provide pre-annotated reference datasets. User-provided reference will be integrated with built-in reference for comprehensive annotation. The input requires two lists of Seurat objects: query list (query datasets) and reference list (reference datasets). Seurat object in query list needs a column named `seurat_clusters`, representing the clustering result, which can be obtained by the `FindClusters` function during [Seurat](https://satijalab.org/seurat/) pipnline. Seurat object in reference list needs two columns named `seurat_clusters` and `scPlantGM.refanno`, seprately. The former is same as the one in query list, and the latter is cell annotation used as a reference.

## Load data
We load `Inflorescence_query` as query list and `Inflorescence_ref` as reference list.
```r
data("Inflorescence_query", package = "scPlantGM")
data("Inflorescence_ref", package = "scPlantGM")
```

Load table of cell type level of target tissue.
```r
data("layer_info_Arabidopis", package = "scPlantGM")
layer_info <- layer_info_Arabidopis %>% dplyr::filter(Organ=='Inflorescence') %>% dplyr::select(Celltype1,Celltype2)
```

## Preprocess
Another necessary and important object is `seuorder`,  which is a vector shows the order of Suerat objects in the list. Strongly recommend using 'Sample' to present the order. Or user can just simply use `c(1,2,3...)` for a quick start.
```r
query <- process_obj(query, seuorder = c("SRX13437605"), type = 'query')
reference <- process_obj(reference, seuorder = c("SRX10918085"), type = 'reference')
```


## Prediction
Run scPlantGM for cell type prediction.
```r
result <- scPlantGM(query, reference, species="Arabidopis", organ="Inflorescence", custom="Semi", layer_info, layer = 2)
```
- `query`: a (list of) query Seurat object that you want to conduct prediction.
- `reference`: a (list of) Seurat object as reference.
- `species`: species of scRNA-seq data. You can only choose within `Arabidopis`, `Maize`, and `Rice` in this mode.
- `organ`: organ of species. Built-in organs: Arabidopis(`Root`, `Leaf`, `Pollen`, `Inflorescence`), Maize(`Root`, `Leaf`), Rice(`Root`, `Leaf`).
- `custom` = `Semi`.
- `layer_info`: a format-fixed data.frame includes layers information (will be ignored if custom='None'). A example shown in data/layer_info_example.rda.
- `layer`: layer of prediction: must be an int (0,1,2...). In `None` and `Semi` mode, user can choose one of built-in layer: Arabidopis(`1`, `2`, `3`), Maize(`1`, `2`), Rice(`1`, `2`).

The prediction result will be stored as a data.frame as follows:
```r
> head(result)
#                              Cell       Cluster            prediction1      prediction2             prediction
#1 SRX13437605@@_AAACCCAAGCCGCACT-1 SRX13437605:0          Ground tissue           Cortex                 Cortex
#2 SRX13437605@@_AAACCCAAGCTTGTTG-1 SRX13437605:1          Ground tissue           Cortex                 Cortex
#3 SRX13437605@@_AAACCCAAGGAACTAT-1          <NA>                   <NA>             <NA>                   <NA>
#4 SRX13437605@@_AAACCCAAGTATGGCG-1 SRX13437605:8        Vascular tissue Vascular cambium       Vascular cambium
#5 SRX13437605@@_AAACCCAAGTCTAGAA-1 SRX13437605:4  Shoot apical meristem  Flower meristem        Flower meristem
#6 SRX13437605@@_AAACCCACAAAGGAGA-1 SRX13437605:3 Shoot system epidermis                / Shoot system epidermis
```

</details>


<details>
<summary>3. Annotate using user-provided reference</summary>

In `All` mode, users need to provide pre-annotated reference datasets. The input requires two lists of Seurat objects: `query` (query datasets) and `reference` (reference datasets). Seurat object in query list needs a column named `seurat_clusters`, representing the clustering result, which can be obtained by the `FindCluster` function during [Seurat](https://satijalab.org/seurat/) pipnline. Seurat object in reference list needs two columns named 'seurat_clusters' and 'scPlantGM.refanno', seprately. The former is same as the one in query list, and the latter is cell annotation used as a reference.

## Load data
We load `Inflorescence_query` as query list and `Inflorescence_ref` as reference list.
```r
data("Inflorescence_query", package = "scPlantGM")
data("Inflorescence_ref", package = "scPlantGM")
```

Load table of cell type level of target tissue.
```r
data("layer_info_Arabidopis", package = "scPlantGM")
layer_info <- layer_info_Arabidopis %>% dplyr::filter(Organ=='Inflorescence') %>% dplyr::select(Celltype1,Celltype2)
```


## Preprocess
Another necessary and important object is `seuorder`,  which is a vector shows the order of Suerat objects in the list. Strongly recommend using 'Sample' to present the order. Or user can just simply use `c(1,2,3...)` for a quick start.
```r
query <- process_obj(query, seuorder = c("SRX13437605"), type = 'query')
reference <- process_obj(reference, seuorder = c("SRX10918085"), type = 'reference')
```

## Prediction
Run scPlantGM for cell type prediction.
```r
result <- scPlantGM(query, reference, custom="All", layer_info, layer = 2)
```
- `query`: a (list of) query Seurat object that you want to conduct prediction.
- `reference`: a (list of) Seurat object with annotation as reference.
- `custom` = `All`.
- `layer_info`: a format-fixed data.frame includes layers information (will be ignored if custom='None'). A example shown in data/layer_info_example.rda.
- `layer`: layer of prediction: must be an int (0,1,2...). In `None` and `Semi` mode, user can choose one of built-in layer: Arabidopis(`1`, `2`, `3`), Maize(`1`, `2`), Rice(`1`, `2`).

The prediction result will be stored as a data.frame as follows:
```r
head(result)
#                              Cell       Cluster            prediction1     prediction2             prediction
#1 SRX13437605@@_AAACCCAAGCCGCACT-1 SRX13437605:0          Ground tissue          Cortex                 Cortex
#2 SRX13437605@@_AAACCCAAGCTTGTTG-1 SRX13437605:1          Ground tissue          Cortex                 Cortex
#3 SRX13437605@@_AAACCCAAGGAACTAT-1          <NA>                   <NA>            <NA>                   <NA>
#4 SRX13437605@@_AAACCCAAGTATGGCG-1 SRX13437605:8        Vascular tissue          Phloem                 Phloem
#5 SRX13437605@@_AAACCCAAGTCTAGAA-1 SRX13437605:4  Shoot apical meristem Flower meristem        Flower meristem
#6 SRX13437605@@_AAACCCACAAAGGAGA-1 SRX13437605:3 Shoot system epidermis               / Shoot system epidermis
```

</details>

## Advanced usage
More comprehensive arguments of scPlantGM
- `p_thres`: a threshold for determining what level of module purity is acceptable. The higher the value, the purer the module are.
- `m_thres`: a threshold for determining how many modules can be accepted. The higher the value, the more modules can be accepted.
- `x_thres`: a threshold for determining whether the cell type should be ‘unknown’.
- `cores`: numbers of CPU for running.


# Citation
Lu, K., et al. A gene module-based framework for plant cell atlas annotation and cross-species marker gene discovery.


# Contact us
If you have any question, suggestion or bug found in the method, feel free to contact us! Email: dijunchen@nju.edu.cn; Lukaiyann@163.com
