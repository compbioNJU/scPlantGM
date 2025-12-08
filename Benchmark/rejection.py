# TOSICA
import os
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import pandas as pd
import pyreadr
import anndata
import time
import torch
import numpy as np

organ = sys.argv[1]

# Load data and preprocess
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')
result = pyreadr.read_r("cluster.overlap.rds")
jaccard_mat = list(result.values())[0]

info1 = info[info['Organ'] == organ]
datasets = info1['Dataset'].unique()

df_data = pd.read_excel("Cell_Type.xlsx", sheet_name="ara")
data1 = df_data[df_data['Organ'] == organ]

for seldataset in datasets:

    info2 = info1[info1['Dataset'] == seldataset]
    
    if organ == "Leaf":
        refdata = data1[data1['Celltype1'].isin(["Mesophyll", "Leaf epidermis"])]
        querydata = data1[data1['Celltype1'] == "Vascular tissue"]
    elif organ == "Root":
        refdata = data1[data1['Celltype1'].isin(["Root epidermis", "Root cortex", "Root cap"])]
        querydata = data1[data1['Celltype1'] == "Root stele"]

    reftype = set(pd.concat([refdata['Celltype1'], refdata['Celltype2'], refdata['Celltype3']]).unique()) - {"/"}
    querytype = set(pd.concat([querydata['Celltype1'], querydata['Celltype2'], querydata['Celltype3']]).unique()) - {"/"}

    refcluster = info2[info2['annotation'].isin(reftype)]['Cluster'].unique().astype(str)
    querycluster = info2[info2['annotation'].isin(querytype)]['Cluster'].unique().astype(str)

    querycluster = np.intersect1d(querycluster, jaccard_mat.columns)

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].tolist()
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].tolist()

    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    ref_data = ov.pp.preprocess(ref_data, mode='shiftlog|pearson', n_HVGs=3000)
    que_data = ov.pp.preprocess(que_data, mode='shiftlog|pearson', n_HVGs=3000)

    # Retain highly variable genes
    ref_data = ref_data[:, ref_data.var.highly_variable_features]
    ref_data.var_names_make_unique()
    que_data.var_names_make_unique()
    ret_gene = list(set(que_data.var_names) & set(ref_data.var_names))
    len(ret_gene)
    que_data = que_data[:, ret_gene]
    ref_data = ref_data[:, ret_gene]

    # train TOSICA
    pro_path = 'TOSICA_' + organ + "_" + str(seldataset) + "/" + organ + "_demo"
    tosica_obj = ov.single.pyTOSICA(adata=ref_data, depth=1, label_name='annotation', project_path=pro_path, batch_size=8)
    tosica_obj.train(epochs=10)
    tosica_obj.save()
    #tosica_obj.load()

    # predict
    pred = tosica_obj.predicted(pre_adata=que_data,batch_size=50)
    acc = pred.obs[(pred.obs["Prediction"] == pred.obs["annotation"])].shape[0] / pred.obs.shape[0]

    result_path = 'TOSICA_' + organ + "_" + str(seldataset) + "_pyresult.csv"
    prediction = pred.obs["Prediction"]
    annotation = pred.obs["annotation"]
    probability = pred.obs["Probability"]

    result = pd.concat([prediction, annotation, probability], axis=1)
    result.to_csv(result_path)


#############################################################################################################################
# CellTypist
import os
import pandas as pd
import omicverse as ov
import scanpy as sc
import celltypist
import time
import sys
import numpy as np
import pyreadr

organ = sys.argv[1]

# Load data and process
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

info1 = info[info['Organ'] == organ]
datasets = info1['Dataset'].unique()

df_data = pd.read_excel("Cell_Type.xlsx", sheet_name="ara")
data1 = df_data[df_data['Organ'] == organ]

for seldataset in datasets:

    info2 = info1[info1['Dataset'] == seldataset]
    
    if organ == "Leaf":
        refdata = data1[data1['Celltype1'].isin(["Mesophyll", "Leaf epidermis"])]
        querydata = data1[data1['Celltype1'] == "Vascular tissue"]
    elif organ == "Root":
        refdata = data1[data1['Celltype1'].isin(["Root epidermis", "Root cortex", "Root cap"])]
        querydata = data1[data1['Celltype1'] == "Root stele"]

    reftype = set(pd.concat([refdata['Celltype1'], refdata['Celltype2'], refdata['Celltype3']]).unique()) - {"/"}
    querytype = set(pd.concat([querydata['Celltype1'], querydata['Celltype2'], querydata['Celltype3']]).unique()) - {"/"}

    refcluster = info2[info2['annotation'].isin(reftype)]['Cluster'].unique().astype(str)
    querycluster = info2[info2['annotation'].isin(querytype)]['Cluster'].unique().astype(str)

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].tolist()
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].tolist()

    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)

    # Training a CellTypist model.
    new_model = celltypist.train(ref_data, labels = ref_data.obs['annotation'])
    # Write out the model.
    file_path = 'CellTypist_' + organ + '_' + seldataset + '_model.pth'

    new_model.write(file_path)
    # The model summary information.
    new_model

    ncell = que_data.shape[1]
    print(ncell)
    if ncell < 100000:
        predictions = celltypist.annotate(que_data, model = new_model)
    elif 100000 <= ncell < 500000:
        predictions = celltypist.annotate(que_data, model = new_model, use_SGD = True)
    else:
        predictions = celltypist.annotate(que_data, model = new_model, use_SGD = True, mini_batch = True)

    predictions.to_table(folder = 'CellTypist_dataset', prefix =  'CellTypist_' + organ + '_' + seldataset + '_')
   


##############################################################################################################################
##############################################################################################################################
# TOSICA
import os
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import pandas as pd
import pyreadr
import anndata
import time
import torch
import numpy as np


# Load data and process
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')
result = pyreadr.read_r("cluster.overlap.rds")
jaccard_mat = list(result.values())[0]


organ1 = sys.argv[1]
organ2 = sys.argv[2]

info1 = info[info['Organ'] == organ1]
info2 = info[info['Organ'] == organ2]

datasets1 = info1['Dataset'].unique()
datasets2 = info2['Dataset'].unique()

res1 = pyreadr.read_r(f"{organ1}_id.rds")
res2 = pyreadr.read_r(f"{organ2}_id.rds")
id1 = res1[None].values.flatten()
id2 = res2[None].values.flatten()

for k in range(10):

    print(k)

    sel1 = datasets1[id1[k]]
    sel2 = datasets2[id2[k]]
    
    info11 = info1[info1['Dataset'] == sel1]
    info22 = info2[info2['Dataset'] == sel2]
    
    print(info11['annotation'].value_counts())
    print(info22['annotation'].value_counts())
    
    refcluster = pd.unique(info11['Cluster']).astype(str)
    querycluster = pd.unique(info2['Cluster']).astype(str)
    querycluster = np.intersect1d(querycluster, jaccard_mat.columns)

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].tolist()
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].tolist()

    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    ref_data = ov.pp.preprocess(ref_data, mode='shiftlog|pearson', n_HVGs=3000)
    que_data = ov.pp.preprocess(que_data, mode='shiftlog|pearson', n_HVGs=3000)

    # Retain highly variable genes
    ref_data = ref_data[:, ref_data.var.highly_variable_features]
    ref_data.var_names_make_unique()
    que_data.var_names_make_unique()
    ret_gene = list(set(que_data.var_names) & set(ref_data.var_names))
    len(ret_gene)
    que_data = que_data[:, ret_gene]
    ref_data = ref_data[:, ret_gene]

    # train TOSICA
    pro_path = 'TOSICA_' + organ1 + "2" + organ2 + "/" + organ1 + "_" + sel1 + "_" + organ2 + "_" + sel2 + "_demo"
    tosica_obj = ov.single.pyTOSICA(adata=ref_data, depth=1, label_name='annotation', project_path=pro_path, batch_size=8)
    tosica_obj.train(epochs=10)
    tosica_obj.save()
    #tosica_obj.load()

    # predict
    pred = tosica_obj.predicted(pre_adata=que_data,batch_size=50)
    acc = pred.obs[(pred.obs["Prediction"] == pred.obs["annotation"])].shape[0] / pred.obs.shape[0]

    result_path = 'TOSICA_' + organ1 + "_" + sel1 + "_" + organ2 + "_" + sel2 + "_pyresult.csv"
    prediction = pred.obs["Prediction"]
    annotation = pred.obs["annotation"]
    probability = pred.obs["Probability"]

    result = pd.concat([prediction, annotation, probability], axis=1)
    result.to_csv(result_path)



#############################################################################################################################
# CellTypist
import os
import pandas as pd
import omicverse as ov
import scanpy as sc
import celltypist
import time
import sys
import numpy as np
import pyreadr


organ1 = sys.argv[1]
organ2 = sys.argv[2]

# Load data and process
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')
info['Cluster'] = info['Cluster_1']

info1 = info[info['Organ'] == organ1]
info2 = info[info['Organ'] == organ2]

datasets1 = info1['Dataset'].unique()
datasets2 = info2['Dataset'].unique()

res1 = pyreadr.read_r(f"{organ1}_id.rds")
res2 = pyreadr.read_r(f"{organ2}_id.rds")
id1 = res1[None].values.flatten()
id2 = res2[None].values.flatten()

for k in range(10):

    print(k)

    sel1 = datasets1[id1[k]-1]
    sel2 = datasets2[id2[k]-1]
    
    info11 = info1[info1['Dataset'] == sel1]
    info22 = info2[info2['Dataset'] == sel2]
    
    print(info11['annotation'].value_counts())
    print(info22['annotation'].value_counts())
    
    refcluster = pd.unique(info11['Cluster']).astype(str)
    querycluster = pd.unique(info2['Cluster']).astype(str)

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].tolist()
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].tolist()

    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)

    # Training a CellTypist model.
    new_model = celltypist.train(ref_data, labels = ref_data.obs['annotation'])
    # Write out the model.
    file_path = 'CellTypist_' + organ1 + "_" + sel1 + "_" + organ2 + "_" + sel2 + '_model.pth'

    new_model.write(file_path)
    # The model summary information.
    new_model

    ncell = que_data.shape[1]
    print(ncell)
    if ncell < 100000:
        predictions = celltypist.annotate(que_data, model = new_model)
    elif 100000 <= ncell < 500000:
        predictions = celltypist.annotate(que_data, model = new_model, use_SGD = True)
    else:
        predictions = celltypist.annotate(que_data, model = new_model, use_SGD = True, mini_batch = True)

    predictions.to_table(folder = 'CellTypist_tissue', prefix =  'CellTypist_' + organ1 + "_" + sel1 + "_" + organ2 + "_" + sel2 + '_')

 