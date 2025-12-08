# Split by dataset-level cross-validation

import os
import runpy
globals().update(runpy.run_path("process.py"))
import omicverse as ov
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import pandas as pd
import pyreadr
import anndata
import time
import torch

organ = sys.argv[1]
k = int(sys.argv[2])

info = pd.read_csv('ara_info.csv')

ref_list = get_sample(info, organ, k, sheet_name="ara_dataset")[0]
que_list = get_sample(info, organ, k, sheet_name="ara_dataset")[1]

# Load data and process
alldata = ov.utils.read('Arabidopsis_clean.h5ad')

ref_data = alldata[alldata.obs['orig.ident'].isin(ref_list), :]
que_data = alldata[alldata.obs['orig.ident'].isin(que_list), :]

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
start_time = time.time()

pro_path = 'TOSICA_result/' + organ + "/" + organ + 'K' + str(k) + "_demo"
tosica_obj = ov.single.pyTOSICA(adata=ref_data, depth=1, label_name='annotation', project_path=pro_path, batch_size=8)
tosica_obj.train(epochs=10)
tosica_obj.save()
#tosica_obj.load()

# predict
pred = tosica_obj.predicted(pre_adata=que_data,batch_size=50)

end_time = time.time()
runtime = end_time - start_time
print(runtime)

acc = pred.obs[(pred.obs["Prediction"] == pred.obs["annotation"])].shape[0] / pred.obs.shape[0]
print(acc)

result_path = 'TOSICA_result/' + organ + "_" + 'K' + str(k) + "_pyresult.csv"
prediction = pred.obs["Prediction"]
annotation = pred.obs["annotation"]
probability = pred.obs["Probability"]

result = pd.concat([prediction, annotation, probability], axis=1)
result.to_csv(result_path)




###############################################################################################
# Split by sample-level cross-validation

import os
import runpy
globals().update(runpy.run_path("process.py"))
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
n = int(sys.argv[2])

alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')
alloc_all = pd.read_excel('ref_query.xlsx', sheet_name="ara_sample")

alldataset = alloc_all['Dataset'][np.where(alloc_all['Organ']==organ)[0]].unique()
dataset = alldataset[n-1]
ks = max(alloc_all['K'][np.where(alloc_all['Dataset']==dataset)[0]])

for k in range(ks):
    k = k + 1
    ref_list = get_sample(info, organ, k, dataset=dataset, sheet="ara_sample")[0]
    que_list = get_sample(info, organ, k, dataset=dataset, sheet="ara_sample")[1]

    # Load data and process
    ref_data = alldata[alldata.obs['orig.ident'].isin(ref_list), :]
    que_data = alldata[alldata.obs['orig.ident'].isin(que_list), :]

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
    start_time = time.time()

    pro_path = 'TOSICA_result/' + organ + "_supp/" + organ + '_' + dataset + '_K' + str(k) + "_demo" 
    tosica_obj = ov.single.pyTOSICA(adata=ref_data, depth=1, label_name='annotation', project_path=pro_path, batch_size=8)
    tosica_obj.train(epochs=10)
    tosica_obj.save()
    #tosica_obj.load()

    # predict
    pred = tosica_obj.predicted(pre_adata=que_data,batch_size=50)

    end_time = time.time()
    runtime = end_time - start_time
    print(k)
    if runtime<60:
        print(round(runtime, 4))
    else:
        print(round(runtime/60, 4))

    acc = pred.obs[(pred.obs["Prediction"] == pred.obs["annotation"])].shape[0] / pred.obs.shape[0]
    print("accuracy:" + str(acc))

    result_path = 'TOSICA_result/TOSICA_' + organ + "_" + 'K' + str(k) + "_pyresult.csv"
    prediction = pred.obs["Prediction"]
    annotation = pred.obs["annotation"]
    probability = pred.obs["Probability"]

    result = pd.concat([prediction, annotation, probability], axis=1)
    result.to_csv(result_path)
