import math
import numpy as np
from libcloudphxx import common as cm

def mole_frac_to_mix_ratio(mole_frac_g, p, M, T, rhod):
    """
    convert mole fraction [1] to mixing ratio [kg/kg dry air]
    """
    return mole_frac_g * p * M / cm.R / T / rhod

def mix_ratio_to_mole_frac(mix_r, p, M, T, rhod):
    """
    convert mixing ratio [kg/kg dry air] to mole fraction [1] 
    """
    return mix_r * cm.R * T * rhod / p / M

def rh_to_rv(RH, T, p):
    """
    convert relative humidity [%] to water vapor mixing ratio [kg/kg]
    """
    return cm.eps * RH * cm.p_vs(T) / (p - RH * cm.p_vs(T))

def rhod_calc(T, p, rv):
    """
    calculate dry air density
    """
    th_0 = T * (cm.p_1000 / p)**(cm.R_d / cm.c_pd)
    return cm.rhod(p, th_0, rv)

def rho_calc(T, p, rv):
    """
    calculate total air density (not dry air density)
    """
    rho = p / T / (rv/(1.+rv) * cm.R_v + 1./(1.+rv) * cm.R_d)
    return rho

def henry_teor(chem, p, T, vol, mixr_g, rhod, conc_H):
    """ 
    calculate theoretical number of moles of chemical species dissolved into cloud droplets - Henry law 
    (per kg of dry air)
    """
    # read the molar mass of chem species in gas phase
    M_g = getattr(cm, "M_"+chem)

    # calculate the correction due to dissociation
    if chem in ["O3", "H2O2"]:
        pH_corr = 1  # do nothing for those that don't dissociate
    elif chem == "SO2":
        K1 = getattr(cm, "K_SO2") * np.exp(getattr(cm, "dKR_SO2")  * (1./T - 1./298))
        K2 = getattr(cm, "K_HSO3")* np.exp(getattr(cm, "dKR_HSO3") * (1./T - 1./298)) 
        pH_corr = 1. + K1 / conc_H + K1 * K2 / conc_H / conc_H
    elif chem == "CO2":
        K1   = getattr(cm, "K_CO2") * np.exp(getattr(cm, "dKR_CO2") * (1./T - 1./298))
        K2   = getattr(cm, "K_HCO3")* np.exp(getattr(cm, "dKR_HCO3") * (1./T - 1./298))
        pH_corr = 1. + K1 / conc_H + K1 * K2 / conc_H / conc_H
    elif chem == "HNO3":
        K = getattr(cm, "K_HNO3") * np.exp(getattr(cm, "dKR_"+chem) * (1./T - 1./298))
        pH_corr = 1. + K / conc_H
    elif chem == "NH3":
        K = getattr(cm, "K_NH3") * np.exp(getattr(cm, "dKR_"+chem) * (1./T - 1./298))
        pH_corr = 1 + K / getattr(cm, "K_H2O") * conc_H
    else:
        raise Exception("chem should be O3, H2O2, SO3, CO2, HNO3, or NH3")

    # correction to Henry constant due to temperature and pH
    H       = getattr(cm, "H_"  +chem)
    dHR     = getattr(cm, "dHR_"+chem)
    henry_T = H * np.exp(dHR * (1./T - 1./298)) * pH_corr

    # dissolved  = partial prsessure * Henry_const * molar mass * drop volume
    partial_prs = mix_ratio_to_mole_frac(mixr_g, p, M_g, T, rhod) * p
    return partial_prs * henry_T * vol

def dissoc_teor(chem, T):
    """ 
    calculate theoretical dissociation constants (as a function of temperature)

    """
    # correction to dissoc constant due to temperature
    K        = getattr(cm, "K_"  +chem)
    dKR      = getattr(cm, "dKR_"+chem)

    return K * np.exp(dKR * (1./T - 1./298))

def log10_size_of_lnr(n_tot, mean_r, lnr, gstdev):
    """
    log-normal size distribution (defined as a function of log_10(radius))

    """
    return n_tot / math.sqrt(2 * math.pi) / math.log(gstdev, 10)\
           * math.exp(-1. * math.pow((lnr - math.log(mean_r, 10)), 2) / 2. / math.pow(math.log(gstdev,10),2))

def diag_n_OH(V, conc_H):
    """
    calculate the number of OH moles
    """
    return cm.K_H2O * V / conc_H

def diag_n_NH3_H2O(n_N3, T, conc_H):
    """
    calculate the number of NH3*H2O moles
    """
    return n_N3 / (1. + dissoc_teor("NH3", T) * conc_H / cm.K_H2O)

def diag_n_NH4(n_N3, T, conc_H):
    """
    calculate the number of NH4+ moles
    """
    return n_N3 * conc_H * dissoc_teor("NH3", T) / cm.K_H2O / (1. + dissoc_teor("NH3", T) * conc_H / cm.K_H2O)

def diag_n_HNO3(n_N5, T, conc_H):
    """
    calculate the number of HNO3 moles
    """
    return n_N5 / (dissoc_teor("HNO3", T) / conc_H + 1)
   
def diag_n_NO3(n_N5, T, conc_H):
    """
    calculate the number of NO3- moles
    """
    return dissoc_teor("HNO3", T) / conc_H * n_N5 / (dissoc_teor("HNO3", T) / conc_H + 1)

def diag_n_CO2_H2O(n_C4, T, conc_H):
    """
    calculate the number of CO2*H2O moles
    """
    return n_C4 / (1 + dissoc_teor("CO2", T) / conc_H + dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H)

def diag_n_HCO3(n_C4, T, conc_H):
    """
    calculate the number of HCO3- moles
    """
    return n_C4 * dissoc_teor("CO2", T) / conc_H \
           / (1 + dissoc_teor("CO2", T) / conc_H + dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H)

def diag_n_CO3(n_C4, T, conc_H):
    """
    calculate the number of CO3-- moles
    """
    return n_C4 * dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H \
           / (1 + dissoc_teor("CO2", T) / conc_H + dissoc_teor("CO2", T) * dissoc_teor("HCO3", T) / conc_H / conc_H)

def diag_n_SO2_H2O(n_S4, T, conc_H):
    """
    calculate the number of SO2*H2O moles
    """
    return n_S4 / (1 + dissoc_teor("SO2", T) / conc_H + dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H)

def diag_n_HSO3(n_S4, T, conc_H):
    """
    calculate the number of HSO3- moles
    """
    return n_S4 * dissoc_teor("SO2", T) / conc_H \
           / (1 + dissoc_teor("SO2", T) / conc_H + dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H)


def diag_n_SO3(n_S4, T, conc_H):
    """
    calculate the number of SO3-- moles
    """
    return n_S4 * dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H \
           / (1 + dissoc_teor("SO2", T) / conc_H + dissoc_teor("SO2", T) * dissoc_teor("HSO3", T) / conc_H / conc_H)

def diag_n_HSO4(n_S6, T, conc_H):
    """
    calculate the number of HSO4- moles
    """
    return n_S6 * conc_H / (conc_H + dissoc_teor("HSO4", T))

def diag_n_SO4(n_S6, T, conc_H):
    """
    calculate the number of SO4-- moles
    """
    return n_S6 * dissoc_teor("HSO4", T) / (conc_H + dissoc_teor("HSO4", T))
