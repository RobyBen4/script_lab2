import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import ott_dati, consistenza, retta, residui_lin, fitODR, x_quadro, valuta_fit

excel="/mnt/c/Users/giova/OneDrive - Università di Napoli Federico II/Lab2exp2.xlsx"
(I_r,V_r,dI_r, dV_r),nomi=ott_dati(excel, foglio=1)

(A,B), cov, (dA,dB)=fitODR(retta, I_r, V_r, dI_r, dV_r)

#ZERO DELLA FUNZIONE
I_cc_r=-B/A

dI_cc_r=np.sqrt((1/A**2)*(dB**2)+(((B**2)/(A**4))*(dA**2)))+2*(-1/A)*(B/A**2)*cov[0,1]

#best fit
retta_best=retta((A,B),I_r)
#Test Chi Quadro
chi=x_quadro(I_r, V_r, dI_r,dV_r, A, B)
ngl = np.size(I_r)-2
chi_rid, esito_test=valuta_fit(chi,ngl,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ridotto ha valore {chi_rid}")
else:
    print(f"Il fit non è buono, il chi quadro ridotto ha valore {chi_rid}")


V_0 =10 #Volt
R_x =150 #ohm
dr_x = 1 #da vedere nella parte del ponte

print(f"il generatore ideale nella configurazione con resistenza interna ε_T vale: ({B:.3f} +- {dB:.3f})V")
print(f"la resistenza totale nella configurazione con resistenza interna R_T vale: ({A:.3f} +- {dA:.3f})Ω")
print(f"la corrente di cortocircuito vale: ({I_cc_r:.3f}+-{dI_cc_r})A")


l1 = consistenza(B, V_0, dB, 0)
if l1<2:
    print("i due valori dei generatori sono compatibili")
elif l1>2 and l1<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")

l2 = consistenza(A, R_x, dA, dr_x)
if l2<2:
    print("i due valori della resistenza sono compatibili")
elif l2>2 and l2<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")

#RESIDUI
residui=residui_lin(I_r,V_r,dI_r,dV_r,A,B)

#GRAFICO PLOT
fig, (ax1,ax2) = plt.subplots(2, 1, figsize = (8,10), sharex=False)

ax1.errorbar(
    I_r, V_r,
    xerr=dI_r, yerr=dV_r,
    fmt="o", color="blue", ecolor="blue",
    capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
    label= "Fit V_r Vs I_r  (Resistenza interna)"
)

x_fit = np.linspace(0, I_r.max(), 200)
y_fit = retta((I_r,V_r), x_fit)
ax1.plot(x_fit,y_fit, "b-", linewidth =1.5)

#GRAFICO RESIDUI
ax2.axhline(0, color="black", linewidth=1)

ax2.scatter(I_r, residui, color="blue", label="Residui del fit")

ax2.set_ylim(-10,10)
ax2.set_xlabel(nomi[0])
ax2.set_ylabel("Residui normalizzati")
ax2.set_title("Grafico dei residui")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()
