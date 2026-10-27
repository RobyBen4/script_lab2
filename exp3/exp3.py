import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

#----CREAZIONE DEL PLOT DELL'ESPERIMENTO----
def attenuante_r(x):
    return (2*np.pi*x*R*C)/(np.sqrt(1+(2*np.pi*x*R*C)^2))
def attenuante_c(x):
    return 1/(np.sqrt(1+(2*np.pi*x*R*C)^2))
def phi_c(x):
    return np.arctan(-2*R*C*np.pi*x)
def phi_r(x):
    return np.arctan(1/(2*R*C*np.pi*x))


figure, (ax1,ax2) = plt.subplots(2,1, figsize=(8,10), sharex=False)

x_fit_A = np.logspace(A_r.min(), A_r.max(), 200)
y_fit_A_c = attenuante_c(x_fit_A)
y_fit_A_r = attenuante_r(x_fit_A)

ax1.plot(x_fit_A,  y_fit_A_r, "-", color="blue", linewidth=1.5)
ax1.scatter(f_ic, A_r,  "-", color="red", linewidth=1.5)

ax1.plot(x_fit_A,  y_fit_A_c, "-", color="blue", linewidth=1.5)
ax1.scatter(f_ic, A_c,  "-", color="red", linewidth=1.5)

y_fit_phi_c = phi_c(x_fit_A)
y_fit_phi_r = phi_r(x_fit_A)

ax1.plot(x_fit_A,  y_fit_phi_r, "-", color="blue", linewidth=1.5)
ax1.scatter(f_ic, A_r,  "-", color="red", linewidth=1.5)

ax2.plot(x_fit_A,  y_fit_phi_c, "-", color="blue", linewidth=1.5)
ax2.scatter(f_ic, A_c,  "-", color="red", linewidth=1.5)

ax1.set_xscale('log')
ax2.set_xscale('log')

ax1.set_ylabel(nomi[1])
ax1.set_title("grafico delle attenuazioni")
ax1.grid(True, alpha=0.3)
ax1.legend()

ax2.set_xlabel(nomi[0])
ax2.set_ylabel("Residui normalizzati")
ax2.set_title("grafico dello sfasamento")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()