"""
Convert GCR, Mcore to radius

@author: Eve J. Lee
May 8th 2020
"""
import numpy as np
import consts
from scipy import integrate as sint
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import pdb
import matplotlib.pyplot as pl

# Constants
gamma = 7./5.
grad_ad = (gamma-1)/gamma
mu = 2.374
rho_e = 5.5 # Earth rho

# From Rogers & Seager (2010)
alpha = 0.68
beta = 0.45
C = -7.32

# Inverse function
def Iinteg(ratio, gam, n):
    return (1./ratio - 1)**(1./(gam-1))*ratio**n

r_ratio = np.logspace(-5, np.log10(0.9999), 1000)

I2 = np.zeros(len(r_ratio))
I1 = np.zeros(len(r_ratio))
for i, r in enumerate(r_ratio):
    I1[i] = sint.quad(Iinteg, r, 1, args=(gamma,1))[0]
    I2[i] = sint.quad(Iinteg, r, 1, args=(gamma,2))[0]

def get_radius(eta, Mc, a, tKH, rhocfac):
    # etas = Matm/Mcore
    # Mc: core mass in Mearth
    # a: orbital distance in AU
    # tKH: cooling time in Myrs

    Rc = (rhocfac * Mc ** 0.25) ** (-1. / 3.) * Mc ** (1. / 3.)

    Teff = consts.Teff_sun*(consts.Rsun/a/consts.au2cm)**0.5
    Rbondi = consts.G*Mc*consts.Mearth/(consts.k*Teff/mu/consts.mH)

    LHS = (eta * Mc * consts.Mearth / 4. / np.pi) * r_ratio ** (3 + 1 / (1 - gamma)) * (Rc * consts.Rearth) ** (
                -3 + 1 / (gamma - 1)) * (grad_ad * Rbondi) ** (1. / (1 - gamma)) / I2
    RHS = (64*np.pi*consts.sb*mu*consts.mH/3/consts.k)*grad_ad
    RHS *= 10**(-C)*(mu*consts.mH/consts.k)**alpha*Teff**(3-alpha-beta)
    RHS *= (I2/I1)*(tKH*consts.s2Myr/eta/Mc/consts.Mearth)*(1/r_ratio)*Rc*consts.Rearth
    f_lhs = interp1d(r_ratio, LHS)
    f_rhs = interp1d(r_ratio, RHS**(1./(1+alpha)))

    try:
        sol = root_scalar(lambda r: f_lhs(r) - f_rhs(r), method='brentq', bracket=[0.00001, 0.9999])
        return sol.root
    except ValueError:
        print("Solution cannot be found for (eta=%4.2e, Mc=%4.2f Mearth, a=%4.2f au, tKH=%i Myrs!" % (eta, Mc, a, tKH))
        return 0

def phot_corr(r_ratio, eta, Mc, a, rhocfac):
    Rc = (rhocfac * Mc ** 0.25) ** (-1./3.) * Mc ** (1./3.)

    Teff = consts.Teff_sun * (consts.Rsun / a / consts.au2cm) ** 0.5
    Rbondi = consts.G * Mc * consts.Mearth / (consts.k * Teff / mu / consts.mH)

    rho_rcb = (eta * Mc * consts.Mearth / 4. / np.pi) * r_ratio ** (3 + 1 / (1 - gamma)) * (Rc * consts.Rearth) ** (
                -3 + 1 / (gamma - 1)) * (grad_ad * Rbondi) ** (1. / (1 - gamma)) / Iinteg(r_ratio, gamma, 2)

    g = consts.G * consts.Mearth * Mc / (Rc * consts.Rearth / r_ratio) ** 2
    kappa = (10 ** C * (2 * g / 3.) ** alpha * Teff ** beta) ** (1. / (1 + alpha))
    rho_ph = (2. / 3.) * mu * consts.mH * g / consts.k / Teff / kappa

    corr_add = np.log(rho_rcb / rho_ph) * consts.k * Teff / mu / consts.mH / g
    if type(corr_add) == type(np.array([])):
        corr_add = np.where(corr_add < 0, 0, corr_add)
    else:
        if corr_add < 0:
            print(eta, Mc)
            corr_add = 0
    Rphot = Rc / r_ratio + corr_add / consts.Rearth

    return Rphot



