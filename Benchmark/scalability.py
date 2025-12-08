# Split by dataset-level cross-validation
# TOSICA
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
import pickle


organ = sys.argv[1]
percent = sys.argv[2]

# Load data and preprocess
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

for kvalue in range(5):
    # kvalue = 2
    kvalue = kvalue + 1

    rds_file1 = f"downsample_ids_{organ}_{percent}_{kvalue}_ref.rds"
    result1 = pyreadr.read_r(rds_file1)
    refcluster = list(result1.values())[0].iloc[:, 0].tolist()

    rds_file2 = f"downsample_ids_{organ}_{percent}_{kvalue}_query.rds"
    result2 = pyreadr.read_r(rds_file2)
    querycluster = list(result2.values())[0].iloc[:, 0].tolist()

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].values
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].values

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
    pro_path = 'TOSICA_' + organ + "_" + str(percent) + "/" + organ + "_" + 'K' + str(kvalue) + "_demo" 
    tosica_obj = ov.single.pyTOSICA(adata=ref_data, depth=1, label_name='annotation', project_path=pro_path, batch_size=8)
    tosica_obj.train(epochs=10)
    tosica_obj.save()
    #tosica_obj.load()

    # predict
    pred = tosica_obj.predicted(pre_adata=que_data,batch_size=50)

    result_path = 'TOSICA_' + organ + "_" + str(percent) + "_" + 'K' + str(kvalue) + "_pyresult.csv"
    prediction = pred.obs["Prediction"]
    annotation = pred.obs["annotation"]
    probability = pred.obs["Probability"]

    result = pd.concat([prediction, annotation, probability], axis=1)
    result.to_csv(result_path)




###############################################################################################################################
# scBalance
import os
import sys
sys.path.append('scBalance_result')
import runpy
globals().update(runpy.run_path("process.py"))
import omicverse as ov
import scBalance as sb
import scBalance.scbalance_IO as ss
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix, cohen_kappa_score
import torch
import torch.nn as nn
from torch.optim import Adam
import torch.utils.data as Data
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
from collections import Counter
import torch.nn.functional as F
import time
import pyreadr


organ = sys.argv[1]
percent = sys.argv[2]


alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

for kvalue in range(5):

    kvalue = kvalue + 1

    rds_file1 = f"downsample_ids_{organ}_{percent}_{kvalue}_ref.rds"
    result1 = pyreadr.read_r(rds_file1)
    refcluster = list(result1.values())[0].iloc[:, 0].tolist()

    rds_file2 = f"downsample_ids_{organ}_{percent}_{kvalue}_query.rds"
    result2 = pyreadr.read_r(rds_file2)
    querycluster = list(result2.values())[0].iloc[:, 0].tolist()

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].values
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].values

    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    #data process
    ref_label = pd.DataFrame(ref_data.obs['annotation'])
    que_label = pd.DataFrame(que_data.obs['annotation'])

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)


    gene = pd.Index.intersection(ref_data.var_names, que_data.var_names)

    ref_matrix = ref_data.to_df()[gene]
    que_matrix = que_data.to_df()[gene]

    ref_label.columns = ['Label']
    ref_label.Label.value_counts()
    status_dict = ref_label['Label'].unique().tolist()
    ref_label['transfromed']=ref_label['Label'].apply(lambda x : status_dict.index(x))
    y_train = ref_label['transfromed'].values

    X_train = ref_matrix.values
    X_test = que_matrix.values

    ref_label.Label.value_counts()

    #status_dict

    dtype = torch.float
    X_train=torch.from_numpy(X_train)
    X_test=torch.from_numpy(X_test)

    dtype = torch.float

    y_train = torch.from_numpy(y_train.to_numpy())

    class_sample_count = torch.tensor(
        [(y_train == t).sum() for t in torch.unique(y_train, sorted=True)])
    weight = 1. / class_sample_count.float()
    samples_weight = torch.tensor([weight[t] for t in y_train])
    sampler = Data.WeightedRandomSampler(samples_weight, len(samples_weight))

    train_data= Data.TensorDataset(X_train,y_train)
    #train_loader= Data.DataLoader(dataset=train_data,batch_size=128,shuffle=True,num_workers=0)
    remainder = 86657 % 128
    if remainder!=1:
        train_loader= Data.DataLoader(dataset=train_data,batch_size=128,shuffle=True,num_workers=1)

    if remainder==1:
        train_loader= Data.DataLoader(dataset=train_data,batch_size=128,shuffle=True,num_workers=1, drop_last=True)

    #### define classfier
    import classifier as cla
        
    input_size = X_train.shape[1]
    num_class = len(status_dict)
    num_epochs = 20

    model = cla.classifier(input_size, num_class)

    print(model)

    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
    total_step = len(train_loader)


    print("--------Start annotating----------")
    start = time.perf_counter()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    #device = 'cpu'
    #device = 'cuda'
    print("Computational unit be used is:", device)

    model.to(device)
    #loss_weight = loss_weight.to(device)


    model.train()
    for epoch in range(num_epochs):
        for step,(batch_x,batch_y) in enumerate(train_loader):
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            outputs = model(batch_x.type(torch.float))
            train_loss = criterion(outputs, batch_y)
            
            # Backward and optimize
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()
            
        a = "="*(epoch+1)
        b = "."*(num_epochs-epoch-1)
        c = ((epoch+1)/num_epochs)*100
        dur = time.perf_counter() - start
        print("\r{:^3.0f}%[{}->{}]{:.2f}s".format(c,a,b,dur),end = '')
        time.sleep(0.1)

    print("--------Annotation Finished----------")

    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        output = model(X_test.type(torch.float))

    temp = F.softmax(output, dim=1)
    #pred = torch.argmax(temp, dim=1).cpu()
    pred = torch.argmax(temp, dim=1).cpu()

    result = list(pred.numpy())
    cluster_result = result
    results = []
    for i in cluster_result:
        results.append(status_dict[i])

    # 保存模型
    modelpath = 'scBalance_' + organ + '_' + str(percent) + '_k' + str(kvalue) + '_model.pth'
    torch.save(model.state_dict(), modelpath)

    #result
    pred_result = pd.DataFrame(que_label)
    pred_result['scBalance_pred'] = results
    result.append(pred_result)
    pred_result.columns = ['annotation', 'scBalance_pred']

    result_path = 'scBalance_' + organ + '_' + str(percent) + '_k' + str(kvalue) + '_result.csv'
    pred_result.to_csv(result_path)



###############################################################################################################################
# CellTypist
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
import pyreadr

organ = sys.argv[1]
percent = sys.argv[2]

alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

for kvalue in range(5):

    kvalue = kvalue + 1

    rds_file1 = f"downsample_ids_{organ}_{percent}_{kvalue}_ref.rds"
    result1 = pyreadr.read_r(rds_file1)
    refcluster = list(result1.values())[0].iloc[:, 0].tolist()

    rds_file2 = f"downsample_ids_{organ}_{percent}_{kvalue}_query.rds"
    result2 = pyreadr.read_r(rds_file2)
    querycluster = list(result2.values())[0].iloc[:, 0].tolist()

    cells1 = info[info['Cluster'].isin(refcluster)]['Cell'].values
    cells2 = info[info['Cluster'].isin(querycluster)]['Cell'].values

    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)

    # Training a CellTypist model.
    new_model = celltypist.train(ref_data, labels = ref_data.obs['annotation'])
    # Write out the model.
    file_path = 'CellTypist_' + organ + '_' + str(percent) + '_k' + str(kvalue) + '_model.pth'
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

    predictions.to_table(folder = 'cluster_downsample', prefix =  'CellTypist_' + organ + '_' + str(percent)  + '_k' + str(kvalue) + '_')
    



##############################################################################################################################################
##############################################################################################################################################
# Split by sample-level cross-validation
# TOSICA
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
import pickle


organ = sys.argv[1]
percent = sys.argv[2]

alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

for kvalue in range(5):
    # kvalue = 2
    kvalue = kvalue + 1

    rds_file1 = f"dataset_downsample_ids_{organ}_{percent}_{kvalue}_ref.rds"
    result1 = pyreadr.read_r(rds_file1)
    refsample = list(result1.values())[0].iloc[:, 0].tolist()

    rds_file2 = f"dataset_downsample_ids_{organ}_{percent}_{kvalue}_query.rds"
    result2 = pyreadr.read_r(rds_file2)
    querysample = list(result2.values())[0].iloc[:, 0].tolist()

    # Load data and process
    cells1 = info[info['Sample'].isin(refsample)]['Cell'].values
    cells2 = info[info['Sample'].isin(querysample)]['Cell'].values

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
    pro_path = 'dataset_TOSICA_' + organ + "_" + str(percent) + "/" + organ + "_" + 'K' + str(kvalue) + "_demo"
    tosica_obj = ov.single.pyTOSICA(adata=ref_data, depth=1, label_name='annotation', project_path=pro_path, batch_size=8)
    tosica_obj.train(epochs=10)
    tosica_obj.save()
    #tosica_obj.load()

    # predict
    pred = tosica_obj.predicted(pre_adata=que_data,batch_size=50)

    result_path = 'dataset_TOSICA_' + organ + "_" + str(percent) + "_" + 'K' + str(kvalue) + "_pyresult.csv"
    prediction = pred.obs["Prediction"]
    annotation = pred.obs["annotation"]
    probability = pred.obs["Probability"]

    result = pd.concat([prediction, annotation, probability], axis=1)
    result.to_csv(result_path)




###############################################################################################################################
# scBalance
import os
import runpy
globals().update(runpy.run_path("process.py"))
import sys
sys.path.append('scBalance_result')
import omicverse as ov
import scBalance as sb
import scBalance.scbalance_IO as ss
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix, cohen_kappa_score
import torch
import torch.nn as nn
from torch.optim import Adam
import torch.utils.data as Data
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
from collections import Counter
import torch.nn.functional as F
import time
import pyreadr


organ = sys.argv[1]
percent = sys.argv[2]

#数据读取与预处理
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

for kvalue in range(5):

    kvalue = kvalue + 1

    rds_file1 = f"dataset_downsample_ids_{organ}_{percent}_{kvalue}_ref.rds"
    result1 = pyreadr.read_r(rds_file1)
    refsample = list(result1.values())[0].iloc[:, 0].tolist()

    rds_file2 = f"dataset_downsample_ids_{organ}_{percent}_{kvalue}_query.rds"
    result2 = pyreadr.read_r(rds_file2)
    querysample = list(result2.values())[0].iloc[:, 0].tolist()

    # 获取细胞名
    cells1 = info[info['Sample'].isin(refsample)]['Cell'].values
    cells2 = info[info['Sample'].isin(querysample)]['Cell'].values

    # 在 AnnData 中子集提取细胞
    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    #data process
    ref_label = pd.DataFrame(ref_data.obs['annotation'])
    que_label = pd.DataFrame(que_data.obs['annotation'])

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)


    gene = pd.Index.intersection(ref_data.var_names, que_data.var_names)

    ref_matrix = ref_data.to_df()[gene]
    que_matrix = que_data.to_df()[gene]

    ref_label.columns = ['Label']
    ref_label.Label.value_counts()
    status_dict = ref_label['Label'].unique().tolist()
    ref_label['transfromed']=ref_label['Label'].apply(lambda x : status_dict.index(x))
    y_train = ref_label['transfromed'].values

    X_train = ref_matrix.values
    X_test = que_matrix.values

    ref_label.Label.value_counts()

    #status_dict

    dtype = torch.float
    X_train=torch.from_numpy(X_train)
    X_test=torch.from_numpy(X_test)

    dtype = torch.float

    y_train = torch.from_numpy(y_train.to_numpy())

    class_sample_count = torch.tensor(
        [(y_train == t).sum() for t in torch.unique(y_train, sorted=True)])
    weight = 1. / class_sample_count.float()
    samples_weight = torch.tensor([weight[t] for t in y_train])
    sampler = Data.WeightedRandomSampler(samples_weight, len(samples_weight))

    train_data= Data.TensorDataset(X_train,y_train)
    #train_loader= Data.DataLoader(dataset=train_data,batch_size=128,shuffle=True,num_workers=0)
    remainder = 86657 % 128
    if remainder!=1:
        train_loader= Data.DataLoader(dataset=train_data,batch_size=128,shuffle=True,num_workers=1)

    if remainder==1:
        train_loader= Data.DataLoader(dataset=train_data,batch_size=128,shuffle=True,num_workers=1, drop_last=True)

    #### define classfier
    import classifier as cla
        
    input_size = X_train.shape[1]
    num_class = len(status_dict)
    num_epochs = 20

    model = cla.classifier(input_size, num_class)

    print(model)

    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
    total_step = len(train_loader)


    print("--------Start annotating----------")
    start = time.perf_counter()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    #device = 'cpu'
    #device = 'cuda'
    print("Computational unit be used is:", device)

    model.to(device)
    #loss_weight = loss_weight.to(device)


    model.train()
    for epoch in range(num_epochs):
        for step,(batch_x,batch_y) in enumerate(train_loader):
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            outputs = model(batch_x.type(torch.float))
            train_loss = criterion(outputs, batch_y)
            
            # Backward and optimize
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()
            
        a = "="*(epoch+1)
        b = "."*(num_epochs-epoch-1)
        c = ((epoch+1)/num_epochs)*100
        dur = time.perf_counter() - start
        print("\r{:^3.0f}%[{}->{}]{:.2f}s".format(c,a,b,dur),end = '')
        time.sleep(0.1)

    print("--------Annotation Finished----------")

    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        output = model(X_test.type(torch.float))

    temp = F.softmax(output, dim=1)
    #pred = torch.argmax(temp, dim=1).cpu()
    pred = torch.argmax(temp, dim=1).cpu()

    result = list(pred.numpy())
    cluster_result = result
    results = []
    for i in cluster_result:
        results.append(status_dict[i])

    modelpath = 'dataset_scBalance_' + organ + '_' + str(percent) + '_k' + str(kvalue) + '_model.pth'
    torch.save(model.state_dict(), modelpath)

    #result
    pred_result = pd.DataFrame(que_label)
    pred_result['scBalance_pred'] = results
    result.append(pred_result)
    pred_result.columns = ['annotation', 'scBalance_pred']

    result_path = 'dataset_scBalance_' + organ + '_' + str(percent) + '_k' + str(kvalue) + '_result.csv'
    pred_result.to_csv(result_path)



##############################################################################################################
# CellTypist
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
import pyreadr

organ = sys.argv[1]
percent = sys.argv[2]

#数据读取与预处理
alldata = ov.utils.read('Arabidopsis_clean.h5ad')
info = pd.read_csv('ara_info.csv')

for kvalue in range(5):

    kvalue = kvalue + 1

    rds_file1 = f"dataset_downsample_ids_{organ}_{percent}_{kvalue}_ref.rds"
    result1 = pyreadr.read_r(rds_file1)
    refcluster = list(result1.values())[0].iloc[:, 0].tolist()

    rds_file2 = f"dataset_downsample_ids_{organ}_{percent}_{kvalue}_query.rds"
    result2 = pyreadr.read_r(rds_file2)
    querycluster = list(result2.values())[0].iloc[:, 0].tolist()

    # 获取细胞名
    cells1 = info[info['Sample'].isin(refcluster)]['Cell'].values
    cells2 = info[info['Sample'].isin(querycluster)]['Cell'].values

    # 在 AnnData 中子集提取细胞
    ref_data = alldata[cells1]
    que_data = alldata[cells2]

    sc.pp.normalize_total(ref_data, target_sum=1e4)
    sc.pp.log1p(ref_data)

    sc.pp.normalize_total(que_data, target_sum=1e4)
    sc.pp.log1p(que_data)

    # Training a CellTypist model.
    new_model = celltypist.train(ref_data, labels = ref_data.obs['annotation'])
    # Write out the model.
    file_path = 'dataset_CellTypist_' + organ + '_' + str(percent) + '_k' + str(kvalue) + '_model.pth'
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

    predictions.to_table(folder = 'cluster_downsample', prefix =  'CellTypist_' + organ + '_' + str(percent)  + '_k' + str(kvalue) + '_')
  