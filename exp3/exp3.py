import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

(V_c, dVc, V_ic, dtc, Tc, f_ic), nomi= ott_dati(/mnt/c/Users/giova/OneDrive - Università di Napoli Federico II/exp3.xlsx)
(V_r, dVr, V_ir, dtr, Tr, f_ir), nomi= ott_dati(/mnt/c/Users/giova/OneDrive - Università di Napoli Federico II/exp3.xlsx,1)

R=1 #da inserire domani
C=1 #da inserire domani
f_uc=(2*np.pi)/Tc
f_ur=(2*np.pi)/Tr
A_c=(V_c/V_ic)
dA_c=(dVc/V_ic)
dA_r=(dVr/V_ir)
A_r=(V_r/V_ir)

phi_c=(-(f_uc)*dtc)
phi_r=(-(f_ur)*dtr)

f_t_teo = 5.0
n_punti = 10

indice_centrale = np.argmin(np.abs(x - f_t_teo))

inizio = max(0, indice_centrale - n_punti // 2)
fine = min(len(x), indice_centrale + (n_punti + 1) // 2)

maschera = np.zeros_like(x, dtype=bool)
maschera[inizio:fine] = True

A_rfit=A_r[maschera]
f_irfit=f_ir[maschera]

A_cfit=A_c[maschera]
f_icfit=f_irfit[maschera]

# 1 è R e 2 è C
A_1, B_1, dA_1, dB_1, cov1 = fit_lin(f_irfit, A_rfit, dA_r) 
A_2, B_2, dA_2, dB_2, cov2= fit_lin(f_icfit, A_cfit, dA_c)

f_t=((B_2-B_1)/(A_1-A_2))



#funzione per avere l'errore su f_t
def calcola_dft_quadrato(A_1, A_2, B_1, B_2, var_A_1, var_A_2, var_B_1, var_B_2, cov_A1_B1, cov_A2_B2):
    
    num = B_2 - B_1
    den = A_1 - A_2
    
    # Primo termine della somma
    termine_1 = (num / (den**2))**2 * (var_A_1 + var_A_2)
    
    # Secondo termine della somma
    termine_2 = (1 / den)**2 * (var_B_1 + var_B_2)
    
    # Terzo termine della somma
    termine_3 = 2 * (num / (den**3)) * (cov_A1_B1 + cov_A2_B2)
    
    # Risultato finale
    dft_squared = termine_1 + termine_2 + termine_3
    
    return dft_squared

df_t=np.sqrt(calcola_dft_quadrato(A_1,A_2,B_1,B_2,dA_1,dA_2,dB_1,dB_2,cov1,cov2))

