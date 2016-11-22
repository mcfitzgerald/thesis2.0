import numpy as np

def init_nanarr(shape):
    nanarr = np.empty(shape)
    nanarr[:] = np.NAN
    return nanarr

def aligner(arr_long,arr_short):
    aligned_arr = init_nanarr(arr_long.shape)
    indices = [] 
    for i in range(len(arr_short)):
        x = np.abs(arr_long - arr_short[i])
        indices.append(np.argmin(x))
        
    for i in range(len(indices)):
        aligned_arr[indices[i]] = arr_short[i]
        
    return aligned_arr
    
#this needs to be cleaned up and made explicit that long array goes into slot 1
def aligner_wrap(arr1, arr2):
    """arr1 is array of largest size in set"""
    if arr1.shape == arr2.shape:
        return arr2
    else:
        return aligner(arr1,arr2)
    
def align_datset(datset):
    '''takes a numpy array of numpy arrays and aligns all arrays to to the first
    largest array it finds in the super array'''
    lng_idx = np.argmax([i.size for i in datset]) #index for first-found longest array
    
    oth_idxs = list(range(len([i.size for i in datset]))) #list of indices
    oth_idxs.remove(lng_idx) #remove index for long array
    
    #print(lng_idx,oth_idxs) PASSED
    
    container = [aligner_wrap(datset[lng_idx],datset[i]) for i in oth_idxs] #call aligner func
    container.insert(0,datset[lng_idx])
    
    return np.array(container)