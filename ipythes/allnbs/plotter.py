import matplotlib.pyplot as plt

def semlog(lig,sat,labs=None):
    if labs is None:
        for i in range(len(sat)):
            plt.semilogx(lig[i],sat[i],'o')
    else:
        for i in range(len(sat)):
            plt.semilogx(lig[i],sat[i],'o', label=labs[i])
            plt.legend(loc=4,prop={'size':6},numpoints=1)
    
    