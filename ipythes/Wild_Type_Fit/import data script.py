#import tab sep values file
df = pd.read_csv('All_WT_for_coopfit_tsv.csv', sep='\t')

#put saturation data column headers into list
colist = df.columns.tolist()[1:]

#create list of dataframes for each [lig]-sat data set for each receptor conc. and ensure sort oder (to avoid potting issues)
dflist = [df[['[egf]', i]].dropna().sort_values(by='[egf]') for i in colist]

#make headers uniform 
colist_low = [colist[i].lower() for i in range(len(colist))]

#apply uniform headers to dataframe
dflist_low = [dflist[i].rename(columns={colist[i]:colist_low[i]}) for i in range(len(colist_low))]

#import total receptor concentrations
with open('rtotswt.csv', 'r') as f:
    g = csv.reader(f)
    rtots = (list(g))
    
#convert to numerical values (floats)
rtots = [float(rtots[i][0]) for i in range(len(rtots))]

EGFR_WT_RTOTS_ord = Ord(zip(colist_low,rtots))

EGFR_WT_DICT_ord = Ord(zip(colist_low,dflist_low))

keylist = list(EGFR_WT_DICT_ord.keys())

keyed_sats = np.array([EGFR_WT_DICT_ord[key][key].as_matrix() for key in keylist])

keyed_ligs = np.array([EGFR_WT_DICT_ord[key]['[egf]'].as_matrix() for key in keylist])

keyed_rtots = np.array([EGFR_WT_RTOTS_ord[key] for key in keylist])