def init_nanarr(shape):
    nanarr = np.empty(shape)
    nanarr[:] = np.NAN
    return nanarr
	
def size_nanarr(arr1, arr2):
    if arr1.shape > arr2.shape:
        new_arr = init_nanarr(arr1.shape)
    else:
        new_arr = init_nanarr(arr1.shape)
        
    return new_arr
    

def aligner(arr1, arr2):
    if arr1.shape > arr2.shape:
        arr_long = arr1
        arr_short = arr2
        aligned_arr = init_nanarr(arr1.shape)
    else:
        arr_long = arr2
        arr_short = arr1
        aligned_arr = init_nanarr(arr2.shape)
        
    indices = []
    
    for i in range(arr_short.shape[0]):
        x = np.abs(arr_long - arr_short[i])
        indices.append(np.argmin(x))
        
    for i in range(len(indices)):
        aligned_arr[indices[i]] = arr_short[i]
        
    return np.array([arr_long,aligned_arr])
        