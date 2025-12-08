def get_sample(info, organ, kvalue, dataset=None, sheet=None):
    
    data = pd.read_excel("ref_query.xlsx", sheet_name=sheet)
    
    if dataset is None:
        data1 = data[(data['Organ'] == organ) & (data['K'] == kvalue)]
    else:
        data1 = data[(data['Dataset'] == dataset) & 
                     (data['Organ'] == organ) & 
                     (data['K'] == kvalue)]
    
    if data1['Dataset'].nunique() > 1:
        dataset1 = data1.loc[data1['RefQuery'] == 'ref', 'Dataset'].unique()
        dataset2 = data1.loc[data1['RefQuery'] == 'query', 'Dataset'].unique()
        
        ref_sample = (info[(info['Organ'] == organ) & 
                       (info['Dataset'].isin(dataset1))]
                  ['Sample'].unique().tolist())
        que_sample = (info[(info['Organ'] == organ) & 
                       (info['Dataset'].isin(dataset2))]
                  ['Sample'].unique().tolist())
    else:
        ref_sample = data1.loc[data1['RefQuery'] == 'ref', 'Sample'].tolist()
        que_sample = data1.loc[data1['RefQuery'] == 'query', 'Sample'].tolist()
    
    return ref_sample, que_sample