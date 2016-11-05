def unibind(parm,lig,rtot):

    k11 = parm[0]
    k21 = parm[1]
    k22 = parm[2]
    l20 = parm[3]

    rfree = (((-1 - k11*lig) + \
    (np.sqrt((1 + k11*lig)**2 + 8*l20*rtot*(1 + k21*lig + \
    k21*k22*(lig**2)))))/(4*l20*(1 + k21*lig + k21*k22*(lig**2))))

    bfrac = (k11*lig + l20*k21*rfree*lig + \
    2*l20*k21*k22*rfree*(lig**2))/(1 + 2*l20*rfree + k11*lig + \
    2*l20*k21*rfree*lig + 2*l20*k21*k22*rfree*(lig**2))