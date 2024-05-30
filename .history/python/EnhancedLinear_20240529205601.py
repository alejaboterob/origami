import numpy as np

def sec(x):
    # Small number to prevent division by zero
    epsilon = 1e-10  
    return 1 / (np.cos(x) + epsilon)

def EnhancedLinear(he, h0, kpi, L0, limlft, limrht):
    limlft = limlft/180*np.pi
    partl = np.pi/limlft
    limrht = limrht/180*np.pi
    partr = np.pi/(2*np.pi-limrht)
    
    if np.ndim(kpi) == 0:
        kpi = np.repeat(kpi, len(he))
    
    Rspr = np.zeros_like(he)
    Kspr = np.zeros_like(he)
    
    Lind = he < limlft
    Rind = he > limrht
    Mind = ~(Lind | Rind)
    if np.any(Lind):
    # Ensure shapes match; reshape should align with Rspr[Lind]'s 1D expected shape
        Rspr[Lind] = kpi[0,Lind] * (np.real(limlft - h0[Lind].reshape(-1))) + \
                 kpi[0,Lind] * np.tan(partl / 2 * (he[Lind].reshape(-1) - limlft)) / (partl / 2)
        Kspr[Lind] = kpi[0,Lind]*sec(partl/2*(he[Lind]-limlft))**2


    
    # Rspr[Lind] = kpi[Lind]*(np.real(limlft - h0[Lind].reshape(-1, 1))) + kpi[Lind]*np.tan(partl/2*(he[Lind].reshape(-1, 1)-limlft))/(partl/2)
    
    if np.any(Rind):
        Rspr[Rind] = kpi[0,Rind]*np.real(limrht-h0[Rind]) + kpi[Rind]*np.tan(partr/2*(he[Rind]-limrht))/(partr/2)
        Kspr[Rind] = kpi[0,Rind]*sec(partr/2*(he[Rind]-limrht))**2
    
    if np.any(Mind):
        Rspr[Mind] = kpi[0,Mind]*np.real(he[Mind][0]-h0[Mind][0])
        Kspr[Mind] = kpi[0,Mind]
    
    Rspr = L0*Rspr
    Kspr = L0*Kspr
    
    if len(Rspr.shape) > 1 and Rspr.shape[1] > 2:
        Espr = np.zeros_like(he)
        Espr[Lind] = 0.5*kpi[Lind]*np.real(h0[Lind]-limlft)**2 + kpi[Lind]*np.real(h0[Lind]-limlft)*(limlft-he[Lind]) - 4*kpi[Lind]/partl**2*np.log(np.abs(np.cos(partl/2*(limlft-he[Lind]))))
        Espr[Rind] = 0.5*kpi[Rind]*np.real(limrht-h0[Rind])**2 + kpi[Rind]*np.real(limrht-h0[Rind])*(he[Rind]-limrht) - 4*kpi[Rind]/partr**2*np.log(np.abs(np.cos(partr/2*(he[Rind]-limrht))))
        Espr[Mind] = 0.5*kpi[Mind]*np.real(he[Mind]-h0[Mind])**2
        Espr = L0*Espr
        return Rspr, Kspr, Espr
    else:
        return Rspr, Kspr

