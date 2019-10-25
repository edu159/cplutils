import numpy as np
from scipy.optimize import curve_fit


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
    return np.log10(cross_power_law(srate, nu0, nuInf, m, n))

def cross_power_law_temp((srate, T), nu0, nuInf, m, n, a1, a2, b1, k):
    # NOTE: T0 seems not to work well if we let it as a free parameters (i.e benzene P=1000)
    T0 = np.min(T) # Representative value
    m = 1.0
    dT = 1.0/T - 1.0/T0
    Tfactor = np.expand_dims(np.exp(a1*dT + a2*dT**2), axis=1)
    return (Tfactor*(nuInf + (nu0 - nuInf)/(1.0+(m*np.outer((b1*T)**k,np.log10(srate+1.0)))**n))).ravel()

def cross_power_law_temp_log((srate, T), nu0, nuInf, m, n, a1, a2, b1, k):
    visc = cross_power_law_temp((srate,T), nu0, nuInf, m, n, a1, a2, b1, k)
    return np.log10(visc)


NAME2MODEL = {"CrossPowerLaw": {"func": cross_power_law, "log_func": cross_power_law_log},
             "CrossPowerLawTemp": {"func": cross_power_law_temp, "log_func": cross_power_law_temp_log}}

def fit_viscosity_model(nu_data, srate, model_name, srate_range=None, gk_value=None, points=10000):
    assert len(nu_data) == len(srate)
    srate_aux = srate
    nu_data_aux = nu_data
    if gk_value is not None:
        srate_aux = np.zeros(len(srate)+1)
        nu_data_aux = np.zeros(len(srate)+1)
        srate_aux[1:] = srate
        srate_aux[0] = 1.0
        nu_data_aux[1:] = nu_data
        nu_data_aux[0] = gk_value
    if model_name not in NAME2MODEL.keys():
        raise Exception("Error: Viscosity model '{}' not knonw!".format(model_name))
    model_func = NAME2MODEL[model_name]["func"]
    model_func_log = NAME2MODEL[model_name]["log_func"]
    if srate_range is None:
        srate_range = [np.min(srate_aux), np.max(srate_aux)]
    axis = np.geomspace(srate_range[0], srate_range[1], points)
    popt, pcov = curve_fit(model_func_log, np.log10(srate_aux), np.log10(nu_data_aux), maxfev=int(1e8), method="dogbox")
    return model_func(np.log10(axis), *popt), axis, popt 


