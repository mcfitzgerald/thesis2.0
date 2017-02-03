i = 2
plt.semilogx(WT_09_lig[i], models.wyman_bestfit(res_wt1.x,WT_08_lig[i],WT_08_rtot[i]),'-',color=colors[i])
plt.semilogx(WT_08_lig[i], models.wyman_bestfit(res_wt1h.x,WT_08_lig[i],WT_08_rtot[i]),'--',color=colors[i])
plt.semilogx(WT_08_lig[i], WT_08_sat[i],'.',color=colors[i])

def viewf(parms1, ligs, data, rtots, index=None):
    """parms is the least_squares result.x"""
    colors = ['b','g','r','c','m','y','k']
    if index is None:
        for i in range(len(ligs)):
            plt.semilogx(ligs[i], models.wyman_bestfit(parms,ligs[i],rtots[i]),'-', color=colors[i])
            plt.semilogx(ligs[i], sats[i],'.',color=colors[i])
    else:
        plt.semilogx(ligs[index], models.wyman_bestfit(parms,ligs[index],rtots[index]),'-', color=colors[index])
        plt.semilogx(ligs[index], sats[index],'.',color=colors[index])
        
    