import numpy as np

# VISCOSITY MODELS 
def bird_carreau(srate, nu0,nuInf, k, n):
    a = 2.0 # By default in OpenFOAM
    return nuInf + (nu0 - nuInf)*(1.0+(k*srate)**a)**((n-1.0)/a)

def bird_carreau_log(srate, nu0, nuInf, k, n):
    return np.log10(bird_carreau(srate, nu0, nuInf, k, n) )

def carreau(srate,nu0, k, n):
    a = 2 # By default in OpenFOAM
    return nu0*(1.0+(k*srate)**a)**((n-1.0)/a)

def carreau_log(srate,nu0, k, n):
    return np.log10(carreau(srate, nu0, k, n) )

def power_law(srate, nuMin, nuMax, k, n):
    if srate < nuMin or srate > nuMax:
        raise Exception("Shear rate provided is out of [%f, %f] interval." % (nuMin, nuMax))
    return k * srate**(n-1) 

def cross_power_law(srate, nu0, nuInf, m, n):
    return nuInf + (nu0 - nuInf)/(1.0+(m*srate)**n)

def cross_power_law_log(srate, nu0, nuInf, m, n):
    return np.log10(nuInf + (nu0 - nuInf)/(1.0+(m*srate)**n))
