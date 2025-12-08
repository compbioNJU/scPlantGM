# Split by dataset-level cross-validation
import os
import runpy
# runpy.run_path("process.py")
globals().update(runpy.run_path("process.py"))
import pandas as pd
import omicverse as ov
import scanpy as sc
import celltypist
import time
import sys


organ = sys.argv[1]
k = int(sys.argv[2])

info = pd.read_csv('ara_info.csv')

ref_list = get_sample(info, organ, k, sheet="ara_dataset")[0]
que_list = get_sample(info, organ, k, sheet="ara_dataset")[1]

# Load data and Preprocess
alldata = ov.utils.read('Arabidopsis_clean.h5ad')

sc.pp.normalize_total(alldata, target_sum=1e4)
sc.pp.log1p(alldata)
sc.pp.highly_variable_genes(alldata, n_top_genes=3000)

# Retain highly variable genes
alldata = alldata[:, alldata.var['highly_variable']].copy()

# Split into reference and query set
ref_data = alldata[alldata.obs['orig.ident'].isin(ref_list), :]
que_data = alldata[alldata.obs['orig.ident'].isin(que_list), :]

# data process
start_time = time.time()

# Training a CellTypist model.
new_model = celltypist.train(ref_data, labels = ref_data.obs['annotation'], check_expression = False)
# Write out the model.
new_model.write('CellTypist_result/CellTypist_' + organ + '_k' + str(k) + '_model.pkl')
# new_model = models.Model.load('CellTypist_result/CellTypist_' + organ + '_k' + str(k) + '_model.pkl')
# The model summary information.
new_model

ncell = que_data.shape[1]
if ncell < 100000:
    predictions = celltypist.annotate(que_data, model = new_model)
elif 100000 <= ncell < 300000:
    predictions = celltypist.annotate(que_data, model = new_model, use_SGD = True)
else:
    predictions = celltypist.annotate(que_data, model = new_model, use_SGD = True, mini_batch = True)

end_time = time.time()
runtime = end_time - start_time
runtime = runtime/60
# print(start_time)
# print(end_time)
print(f"Running time: {runtime:.2f} min") # 输出运算时间

predictions.to_table(folder = 'CellTypist_result', prefix =  'CellTypist_' + organ + '_k' + str(k) + '_')



###############################################################################################
# Split by sample-level cross-validation
import os
import runpy
globals().update(runpy.run_path("process.py"))
import pandas as pd
import omicverse as ov
import scanpy as sc
import celltypist
import time
import sys
import numpy as np

alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')
alloc_all = pd.read_excel('ref_query.xlsx', sheet_name="ara_sample")


organ = sys.argv[1]
n = sys.argv[2]
n = int(n)

alldataset = alloc_all['Dataset'][np.where(alloc_all['Organ']==organ)[0]].unique()
dataset = alldataset[n-1]
ks = max(alloc_all['K'][np.where(alloc_all['Dataset']==dataset)[0]])

for k in range(ks):
    
    k = k + 1

    ref_list = get_sample(info, organ, k+1, dataset=dataset, sheet="ara_sample")[0]
    que_list = get_sample(info, organ, k+1, dataset=dataset, sheet="ara_sample")[1]

    # Load data and Preprocess
    ref_data = alldata[alldata.obs['orig.ident'].isin(ref_list), :]
    que_data = alldata[alldata.obs['orig.ident'].isin(que_list), :]

    ref_data = ov.pp.preprocess(ref_data, mode='shiftlog|pearson', n_HVGs=3000)
    que_data = ov.pp.preprocess(que_data, mode='shiftlog|pearson', n_HVGs=3000)

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)

    # data process
    start_time = time.time()

    # Training a CellTypist model.
    new_model = celltypist.train(ref_data, labels = ref_data.obs['annotation'])
    # Write out the model.
    new_model.write('CellTypist_result/CellTypist_' + organ + '_' + dataset + '_k' + str(k) + '_model.pkl')
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

    end_time = time.time()
    runtime = end_time - start_time
    # print(start_time)
    # print(end_time)
    print(f"Running time: {runtime:.2f} s")

    predictions.to_table(folder = 'CellTypist_result', prefix =  'CellTypist_' + organ + '_' + dataset + '_k' + str(k) + '_')


