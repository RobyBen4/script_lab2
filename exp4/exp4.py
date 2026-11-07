import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

(V_r, dVr, V_i, dV_i, dt, Tr, f_i, Dt, DT), nomi= ott_dati("/mnt/c/Users/giova/OneDrive - Università di Napoli Federico II/exp4.xlsx")

#DATI DEL LAB
R= #da inserire domani
Cap= #da inserire domani

f_i= f_i*1000

f_ur=(2*np.pi)/Tr
df_ur =(2*np.pi*DT_r)/(Tr**2)

A_r=(V_r/V_i)
dA_r = np.sqrt(
    (dVr/V_i)**2 +
    (V_r*dV_i/V_i**2)**2
)

phi_r=(-(f_ur)*dtr)

dphi_r = np.sqrt(
    (dt*df_ur)**2 +
    (f_ur*Dtr)**2
)

f_t_teo = 10.8*1000 # DA DETERMINARE IN LAB
n_punti = 10

indice_centrale = np.argmin(np.abs(f_i - f_t_teo))

inizio = max(0, indice_centrale - n_punti // 2)
fine = min(len(f_i), indice_centrale + (n_punti + 1) // 2)

maschera = np.zeros_like(f_i, dtype=bool)
maschera[inizio:fine] = True

A_rfit=A_r[maschera]
f_ifit=f_i[maschera]
phi_rfit=phi_r[maschera]

(A,B,C),cov,(d_A,d_B,d_C)=fit(parabola,f_ifit,A_rfit,dA_r)
m,q,dm,dq,cov_p= fit_lin(f_ifit,phi_rfit,dphi_r)




A_m=((4*A*B)-(B**2))/(4*C)
w_0_Att=(-B)/(2*C)
Q_Att=(B)/(np.sqrt(2*((B**2)-(4*A*C))))
Q_phi=B/2
w_0_phi=-A/2

d_Am_quadro = (
    d_A**2 
    + (B**2 / (4 * C**2)) * d_B**2 
    + (B**4 / (16 * C**4)) * d_C**2 
    - (B / C) * cov[1,0] 
    + (B**2 / (2 * C**2)) * cov[2,0] 
    - (B**3 / (4 * C**3)) * cov[2,1]
)

d_Am=np.sqrt(d_Am_quadro)

d_Q_quadro_Att = (1 / (B**2 - 4 * A * C)**3) * (
    2 * B**2 * A**2 * d_C**2 
    + 8 * A**2 * C**2 * d_B**2 
    + 2 * B**2 * C**2 * d_A**2 
    + 4 * C * B**2 * A * cov[2,0] 
    - 8 * C * B * A**2 * cov[2,1] 
    - 8 * C**2 * B * A * cov[1,0]
)

d_Q_Att=np.sqrt(d_Q_quadro_Att)

d_w0_Att_quadro = (
    (B**2 / (4 * C**4)) * d_C**2 
    + (1 / (4 * C**2)) * d_B**2 
    - (B / (2 * C**3)) * cov[2,1]
)

d_w0_Att=np.sqrt(d_w0_Att_quadro)

d_w0_quadro_phi = (
    (1 / m**2) * dq**2 
    + (q**2 / m**4) * dm**2 
    + 2 * (-1 / m) * (q / m**2) * cov_p
)

d_w0_phi=np.sqrt(d_w0_quadro_phi)

d_Q_quadro_phi = (dq**2) / 4

d_Q_phi=np.sqrt(d_Q_quadro_phi)

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






