import numpy as np
def sec(x):
    '''The function `sec(x)` calculates the secant of the input angle `x` while preventing division by
    zero.

    Parameters
    ----------
    x
        The function `sec(x)` calculates the secant of the input angle `x`. The parameter `x` represents
    the angle in radians for which you want to calculate the secant.

    Returns
    -------
        The function `sec(x)` returns the secant of the input value `x`, which is calculated as 1 divided
    by the cosine of `x` plus a small epsilon value to prevent division by zero errors.

    '''
    # Small number to prevent division by zero
    epsilon = 1e-10  
    return 1 / (np.cos(x) + epsilon)

def EnhancedLinear(he, h0, kpi, L0, limlft, limrht):
    '''The `EnhancedLinear` function calculates the spring properties (Rspr, Kspr, Espr) based on input
    parameters and conditions.
    
    Parameters
    ----------
    he
        The `he` parameter in the `EnhancedLinear` function represents the input values for which you want
    to calculate the corresponding output values. It is an array of values representing the input
    positions.
    h0
        The `h0` parameter in the `EnhancedLinear` function represents the initial height values for each
    element in the input array `he`. It is used in the calculations to determine the behavior of the
    system based on the input heights `he`.
    kpi
        The `kpi` parameter in the `EnhancedLinear` function represents the stiffness coefficient for the
    linear spring model. It is used to calculate the spring properties based on the input parameters and
    conditions provided in the function.
    L0
        It seems like the description of the parameter L0 is missing. Could you please provide more
    information about what L0 represents in the context of the EnhancedLinear function?
    limlft
        It seems like the code snippet you provided defines a function named `EnhancedLinear` that
    calculates some values based on the input parameters. The parameters used in the function are:
    limrht
        It seems like you were about to provide a description of the parameter `limrht` but the message got
    cut off. Could you please provide more information or let me know if you need assistance with
    something specific related to the `limrht` parameter?
    
    Returns
    -------
        The function `EnhancedLinear` returns three arrays: `Rspr`, `Kspr`, and `Espr`. These arrays
    contain the calculated values based on the input parameters and the logic implemented in the
    function.
    
    '''
    limlft = limlft/180*np.pi
    partl = np.pi/limlft
    limrht = limrht/180*np.pi
    partr = np.pi/(2*np.pi-limrht)
    
    if np.ndim(kpi) == 0:
        kpi = np.repeat(kpi, np.size(he))

    if np.size(h0) != np.size(he): #OJOI 
        h0 = h0[0]

    if np.size(kpi) != 0: #OJOI
        if np.size(kpi) != 1:
            kpi = kpi[0]
    
    Rspr = np.zeros(np.size(he))
    Kspr = np.zeros(np.size(he))
    Espr = np.zeros(np.size(he))

    Lind = he < limlft
    Rind = he > limrht
    Mind = ~(Lind | Rind)
    # h0[Lind]

    if np.any(Lind):
        Rspr[Lind] = kpi*(np.real(limlft - h0)) + kpi*np.tan(partl / 2 * (he - limlft)) / (partl / 2)    
        Kspr[Lind] = kpi*sec(partl/2*(he-limlft))**2
        Espr[Lind] = 0.5*kpi*np.real(h0-limlft)**2 + kpi*np.real(h0-limlft)*(limlft-he) - 4*kpi/partl**2*np.log(np.abs(np.cos(partl/2*(limlft-he))))

    
    if np.any(Rind):
        Rspr[Rind] = kpi*np.real(limrht-h0[Rind].reshape(-1)) + kpi*np.tan(partr/2*(he[Rind].reshape(-1)-limrht))/(partr/2)
        Kspr[Rind] = kpi*sec(partr/2*(he[Rind]-limrht))**2
        Espr[Rind] = 0.5*kpi*np.real(limrht-h0[Rind])**2 + kpi*np.real(limrht-h0[Rind])*(he[Rind]-limrht) - 4*kpi/partr**2*np.log(np.abs(np.cos(partr/2*(he[Rind]-limrht))))

    if np.any(Mind):
        # Rspr[Mind] = kpi[Mind]*np.real(he[Mind]-h0[Mind])
        Rspr = kpi*np.real(he-h0)
        Kspr = kpi
        Espr = 0.5*kpi*np.real(he-h0)**2

    Rspr = L0*Rspr
    Kspr = L0*Kspr
    Espr = L0*Espr
    return Rspr, Kspr, Espr



# import numpy as np

# def sec(x):
#     # Small number to prevent division by zero
#     epsilon = 1e-10  
#     return 1 / (np.cos(x) + epsilon)

# def EnhancedLinear(he, h0, kpi, L0, limlft, limrht):
#     limlft = limlft/180*np.pi
#     partl = np.pi/limlft
#     limrht = limrht/180*np.pi
#     partr = np.pi/(2*np.pi-limrht)
    
#     if np.ndim(kpi) == 0:
#         kpi = np.repeat(kpi, np.size(he))
    
#     Rspr = np.zeros_like(he)
#     Kspr = np.zeros_like(he)
    
#     Lind = he < limlft
#     Rind = he > limrht
#     Mind = ~(Lind | Rind)
#     if np.any(Lind):
#     # Ensure shapes match; reshape should align with Rspr[Lind]'s 1D expected shape
#         Rspr[Lind] = kpi[0,Lind] * (np.real(limlft - h0[Lind].reshape(-1))) + \
#                  kpi[0,Lind] * np.tan(partl / 2 * (he[Lind].reshape(-1) - limlft)) / (partl / 2)
#         Kspr[Lind] = kpi[0,Lind]*sec(partl/2*(he[Lind]-limlft))**2


    
#     # Rspr[Lind] = kpi[Lind]*(np.real(limlft - h0[Lind].reshape(-1, 1))) + kpi[Lind]*np.tan(partl/2*(he[Lind].reshape(-1, 1)-limlft))/(partl/2)
    
#     if np.any(Rind):
#         Rspr[Rind] = kpi[0,Rind]*np.real(limrht-h0[Rind].reshape(-1)) + kpi[Rind]*np.tan(partr/2*(he[Rind].reshape(-1)-limrht))/(partr/2)
#         Kspr[Rind] = kpi[0,Rind]*sec(partr/2*(he[Rind]-limrht))**2
    
#     if np.any(Mind):
#         Rspr[Mind] = kpi[0,Mind]*np.real(he[Mind][0]-h0[Mind][0])
#         Kspr[Mind] = kpi[0,Mind]
    
#     Rspr = L0*Rspr
#     Kspr = L0*Kspr
    
#     if len(Rspr.shape) > 1 and Rspr.shape[1] > 2:
#         Espr = np.zeros_like(he)
#         Espr[Lind] = 0.5*kpi[Lind]*np.real(h0[Lind]-limlft)**2 + kpi[Lind]*np.real(h0[Lind]-limlft)*(limlft-he[Lind]) - 4*kpi[Lind]/partl**2*np.log(np.abs(np.cos(partl/2*(limlft-he[Lind]))))
#         Espr[Rind] = 0.5*kpi[Rind]*np.real(limrht-h0[Rind])**2 + kpi[Rind]*np.real(limrht-h0[Rind])*(he[Rind]-limrht) - 4*kpi[Rind]/partr**2*np.log(np.abs(np.cos(partr/2*(he[Rind]-limrht))))
#         Espr[Mind] = 0.5*kpi[Mind]*np.real(he[Mind]-h0[Mind])**2
#         Espr = L0*Espr
#         return Rspr, Kspr, Espr
#     else:
#         return Rspr, Kspr

