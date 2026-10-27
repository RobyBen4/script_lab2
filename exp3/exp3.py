import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

(V_c, dVc, V_ic, dtc, Tc, f_ic), nomi= ott_dati("/mnt/c/Users/robyb/OneDrive - Università di Napoli Federico II/exp3.xlsx")
(V_r, dVr, V_ir, dtr, Tr, f_ir), nomi= ott_dati("/mnt/c/Users/robyb/OneDrive - Università di Napoli Federico II/exp3.xlsx",1)

#DATI DEL LAB
R=1 #da inserire domani
C=1 #da inserire domani
A_c_vero= 1/np.sqrt(2)
A_r_vero = 1/np.sqrt(2)
phi_c_vero = -np.pi/4
phi_r_vero = np.pi/4


f_uc=(2*np.pi)/Tc
f_ur=(2*np.pi)/Tr
A_c=(V_c/V_ic)
dA_c=(dVc/V_ic)
dA_r=(dVr/V_ir)
A_r=(V_r/V_ir)

phi_c=(-(f_uc)*dtc)
phi_r=(-(f_ur)*dtr)

f_t_teo = 5.0  # DA DETERMINARE IN LAB
n_punti = 10

indice_centrale = np.argmin(np.abs(f_ic - f_t_teo))

inizio = max(0, indice_centrale - n_punti // 2)
fine = min(len(f_ic), indice_centrale + (n_punti + 1) // 2)

maschera = np.zeros_like(f_ic, dtype=bool)
maschera[inizio:fine] = True



A_rfit=A_r[maschera]
f_irfit=f_ir[maschera]



A_cfit=A_c[maschera]
f_icfit=f_irfit[maschera]



# 1 è R e 2 è C
A_1, B_1, dA_1, dB_1, cov1 = fit_lin(f_irfit, A_rfit, dA_r) 
A_2, B_2, dA_2, dB_2, cov2= fit_lin(f_icfit, A_cfit, dA_c)


#---TEST CHI2 PER LE RETTE---
retta1 = retta((A_1, B_1), f_ir)
retta2 = retta((A_2,B_2), f_ic)

chi_1 = chi_quadro(retta1, A_rfit, dA_r)
ngl = np.size(A_cfit)-2
chi_rid, esito_test=valuta_fit(chi_1,ngl,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ridotto ha valore {chi_rid}")
else:
    print(f"Il fit non è buono, il chi quadro ridotto ha valore {chi_rid}")

chi_2 = chi_quadro(retta2, A_cfit, dA_c)
ngl = np.size(R)-2
chi_rid, esito_test=valuta_fit(chi_2,ngl,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ridotto ha valore {chi_rid}")
else:
    print(f"Il fit non è buono, il chi quadro ridotto ha valore {chi_rid}")

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

print(f"La frequenza dove avviene l'intersezione è ({f_t:.3f}+-{df_t:.3f}) ")

#STIMA DI A E PHI CON FT GIUSTA

n_punti = 10
indice_centrale = np.argmin(np.abs(f_ic - f_t))

inizio = max(0, indice_centrale - n_punti // 2)
fine = min(len(f_ic), indice_centrale + (n_punti + 1) // 2)

mask = np.zeros_like(f_ic, dtype=bool)
mask[inizio:fine] = True

phi_r_intorno = phi_r[mask]
A_r_intorno=A_r[mask]
phi_r_s = np.mean(phi_r_intorno)
A_r_s = np.mean(A_r_intorno)

phi_c_intorno = phi_c[mask]
A_c_intorno=A_c[mask]
A_c_s = np.mean(A_c_intorno)
phi_c_s = np.mean(phi_c_intorno) 

print(f"Nel circuito RC i valori sperimentali di A e phi sono: {A_c_s:.3f} e {phi_c_s:.3f}" )
print(f"Nel circuito CR i valori sperimentali di A e phi sono: {A_r_s:.3f} e {phi_r_s:.3f}" )

# CONSISTENZA
l1 = consistenza(A_c_s, A_c_vero, dA_c, 0)
if l1<2:
    print("i due valori della attenuazione di RC sono compatibili")
elif l1>2 and l1<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")

l2= consistenza(A_r_s, A_r_vero, dA_r, 0)
if l2<2:
    print("i due valori della attenuazione di CR sono compatibili")
elif l2>2 and l2<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")


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

x_fit_A = np.logspace(A_r.min(), A_r.max(), 200)/f_t
y_fit_A_c = attenuante_c(x_fit_A)
y_fit_A_r = attenuante_r(x_fit_A)
x_retta1 = np.logspace(A_rfit.min(), A_rfit.max(), 100)/f_t
x_retta2 = np.logspace(A_cfit.min(), A_cfit.max(), 100)/f_t
y_retta1 = retta((A_1, B_1), x_retta1)
y_retta2 = retta((A_2,B_2), x_retta2)


ax1.errorbar(f_ir, A_r
         ,yerr=dA_r,
        fmt="o", color="blue", ecolor="black",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Dati (conf. a valle)"
        )
ax1.plot(x_fit_A,  y_fit_A_r, "-", color="blue", linewidth=1.5, label="Attuenuazione in CR")
ax1.scatter(f_ic, A_r,  "-", color="blue", linewidth=1.5, label= "Punti sperimentali in CR")
ax1.plot(x_retta1, y_retta1, "k-", linewidth=2, label= "Retta nell'intorno di f_t in CR")
ax1.plot(x_retta2, y_retta2, "k-", linewidth=2, label= "Retta nell'intorno di f_t in RC")
ax1.scatter(f_t, A_c_s, color='black')
plt.annotate(f'({f_t:.2f}, {A_c_s:.2f})',
             (f_t, A_c_s),
             xytext=(10,10),
             textcoords='intersezione')

ax1.errorbar(
        f_ic, A_c,
        yerr=dv_v,
        fmt="o", color="red", ecolor="red",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Dati (conf. a valle)"
        )
ax1.plot(x_fit_A,  y_fit_A_c, "-", color="red", linewidth=1.5, label="Attenuazione in RC")
ax1.scatter(f_ic, A_c,  "-", color="red", linewidth=1.5, label="Punti sperimentali in RC")

y_fit_phi_c = phi_c(x_fit_A)
y_fit_phi_r = phi_r(x_fit_A)

ax1.plot(x_fit_A,  y_fit_phi_r, "-", color="blue", linewidth=1.5,label="Sfasamento in CR")
ax1.scatter(f_ic, A_r,  "-", color="blue", linewidth=1.5)

ax2.plot(x_fit_A,  y_fit_phi_c, "-", color="red", linewidth=1.5, label="Sfasamento in RC")
ax2.scatter(f_ic, A_c,  "-", color="red", linewidth=1.5)
ax2.scatter(f_t, phi_r_s, color='black', label="sfasamento massimo in CR")
ax2.scatter(f_t, phi_c_s, color='black', label = "sfasamento masssico in RC")

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