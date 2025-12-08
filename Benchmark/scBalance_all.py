# Split by dataset-level cross-validation

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
import pickle

organ = sys.argv[1]
k = int(sys.argv[2])


info = pd.read_csv('ara_info.csv')

ref_list = get_sample(info, organ, k, sheet_name="ara_dataset")[0]
que_list = get_sample(info, organ, k, sheet_name="ara_dataset")[1]


alldata = ov.utils.read('Arabidopsis_clean.h5ad')

ref_data = alldata[alldata.obs['orig.ident'].isin(ref_list), :]
que_data = alldata[alldata.obs['orig.ident'].isin(que_list), :]

#data process
start_time = time.time()

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

####define classfier
import classifier as cla
    
input_size = X_train.shape[1]
num_class = len(status_dict)
num_epochs = 20

model = cla.classifier(input_size, num_class)

print(model)

criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
total_step = len(train_loader)

import time
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

end_time = time.time()
runtime = end_time - start_time
print(runtime)

modelpath = 'scBalance_result/scBalance_' + organ + '_k' + str(k) + '_model.pth'
torch.save(model.state_dict(), modelpath)

#result
result = []
pred_result = pd.DataFrame(que_label)
pred_result['scBalance_pred'] = results
result.append(pred_result)
pred_result.columns = ['annotation', 'scBalance_pred']

import pickle
filepath = '/mnt/public3/luky/plants/methods/scBalance_result/scBalance_' + organ + '_k' + str(k) + '_result.pkl'
with open(filepath, 'wb') as f:
    pickle.dump(result, f)



###################################################################################################################
# Split by sample-level cross-validation

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

    ref_data = alldata[alldata.obs['orig.ident'].isin(ref_list), :]
    que_data = alldata[alldata.obs['orig.ident'].isin(que_list), :]

    #data process
    start_time = time.time()

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

    ####define classfier
    import classifier as cla
        
    input_size = X_train.shape[1]
    num_class = len(status_dict)
    num_epochs = 20

    model = cla.classifier(input_size, num_class)

    print(model)

    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
    total_step = len(train_loader)

    import time
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

    end_time = time.time()
    runtime = end_time - start_time
    print(start_time)
    print(end_time)
    print(runtime)

    modelpath = 'scBalance_result/scBalance_sample_' + organ + '_' + dataset + '_k' + str(k) + '_model.pth'
    torch.save(model.state_dict(), modelpath)

    #result
    result = []
    pred_result = pd.DataFrame(que_label)
    pred_result['scBalance_pred'] = results
    result.append(pred_result)
    pred_result.columns = ['annotation', 'scBalance_pred']

    import pickle
    filepath = 'scBalance_result/scBalance_sample_' + organ + '_' + dataset + '_k' + str(k) + '_result.pkl'
    with open(filepath, 'wb') as f:
        pickle.dump(result, f)
