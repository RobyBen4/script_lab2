import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

def risonanza(R,L,C,w):
    return R/(np.sqrt(R**2+(L*w-1/(C*w))**2))

def ampiezza_ris(R,L,C,w):
    return -np.arctan((L*w-1/(C*w))/R)


figure, (ax1,ax2) = plt.subplots(2,1, figsize=(8,10), sharex=False)

x_parab = np.logspace(np.log10(f_irfit.min()),np.log10(f_irfit.max()),200)
x_retta1 = np.logspace(np.log10(f_irfit.min()),np.log10(f_irfit.max()),200)
x_fit_A = np.logspace(np.log10(f_ir.min()),np.log10(f_ir.max()),200)
y_fit = risonanza(R,L,C,2*np.pi*f_i)
y_parab = parabola((A,B,C),x_parab)
y_retta=retta((m,q),x_retta)

ax1.errorbar(f_i/f_t, A_r
         ,yerr=dA_r,
        fmt="o", color="blue", ecolor="black",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Dati (conf. a valle)"
        )
        
ax1.plot(x_fit_A/f_t,  y_fit, "-", color="blue", linewidth=1.5, label="Attenuazione")
ax1.scatter(f_i/f_t, A_r,  color="blue", linewidth=1.5 )
ax1.plot(x_parab/f_t, y_parab, "k-", linewidth=2, label= "Parabola nell'intorno di f_t")

ax2.errorbar(f_i/f_t, phi_r
         ,yerr= dphi_r,
        fmt="o", color="blue", ecolor="black",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Punti sperimentali in CR"
        )
ax2.plot(x_fit_A/f_t,  y_fit_phi_r, "-", color="blue", linewidth=1.5,label="Sfasamento")
ax2.scatter(f_i/f_t, phi_c, color="blue", linewidth=1.5 )
ax2.plot(x_retta/f_t, y_retta, "k-", linewidth=2, label= "Retta nell'intorno di f_t")

ax1.set_xscale('log')
ax2.set_xscale('log')

ax1.set_ylabel("A")
ax1.set_title("grafico delle attenuazioni")
ax1.grid(True, alpha=0.3)
ax1.legend()

ax2.set_xlabel("f/f_t")
ax2.set_ylabel("phi")
ax2.set_title("grafico dello sfasamento")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()