import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import ott_dati, consistenza, retta, residui_y, fit_lin, chi_quadro, valuta_fit, dev_std_unif

excel="/mnt/c/Users/robyb/OneDrive - Università di Napoli Federico II/exp_2.xlsx"
(R,I,dI),nomi=ott_dati(excel)

A, B, dA, dB, cov=fit_lin(R, I, dev_std_unif(dI))

R_star=-B/A

dR_star=np.sqrt((1/A*2)*(dB**2)+(((B**2)/(A**4))*(dA**2)))+2*(-1/A)*(B/A**2)*cov

R_1=218.4
R_2=117.7
dR_1= 1
dR_2= 1.8

Rx=(R_2/R_1)*R_star
d_rel_R1 = dR_1/R_1 
d_rel_R2 = dR_2/R_2
d_rel_R_star = dR_star/R_star
dRx=np.sqrt((((d_rel_R1)**2)+(d_rel_R2)**2)+((d_rel_R_star)**2))*np.abs(Rx)

print(R_star, dR_star)
#best fit
retta_best=retta((A,B),R)
#Test Chi Quadro
chi=chi_quadro(retta_best, I, dev_std_unif(dI))
ngl = np.size(R)-2
chi_rid, esito_test=valuta_fit(chi,ngl,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ridotto ha valore {chi_rid}")
else:
    print(f"Il fit non è buono, il chi quadro ridotto ha valore {chi_rid}")


print(f"l'errore relativo di R1 è {d_rel_R1:.3f}")
print(f"l'errore relativo di R2 è {d_rel_R2:.3f}")
print(f"l'errore relativo di R_star è {d_rel_R_star:.3f}")
print(f"la resistenza vale: ({Rx:.3f} +- {dRx:.3f})")


R_vera = 148 #DA VEDERE IN LAB
dR_vera = 0.1
l = consistenza(R_star,R_vera, dR_star, dR_vera)
if l<2:
    print("i due valori sono compatibili")
elif l>2 and l<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")

#RESIDUI
residui=residui_y(I, retta_best, dev_std_unif(dI))

#GRAFICO PLOT
fig, (ax1,ax2) = plt.subplots(2,1, figsize = (8,10), sharex=False)

ax1.errorbar(
    R, I,
    yerr=dI,
    fmt="o", color="red", ecolor="black",
    capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
    label= "Fit I Vs R"
)
ax1.axvline(R_star, color= "black", linestyle= "--")
ax1.axhline(0,color="black", linestyle= "--")
x_fit = np.linspace(R.min(), R.max(), 200)
y_fit = retta((A,B), x_fit)
ax1.plot(x_fit,y_fit, "b-", linewidth =1.5, zorder=1)
ax1.scatter(R_star , retta((A,B),R_star), color="black", s=60 , label="Punto R*", zorder=3)
ax1.set_xlabel(nomi[0])
ax1.set_ylabel(nomi[1])
ax1.set_title("Fit I Vs R")
ax1.grid(True, alpha =0.3)
ax1.legend()

#GRAFICO RESIDUI
ax2.axhline(0, color="black", linewidth=1)

ax2.scatter(R, residui, color="red", label="Residui del fit")

ax2.set_ylim(-2,2)
ax2.set_xlabel(nomi[0])
ax2.set_ylabel("Residui normalizzati")
ax2.set_title("Grafico dei residui")
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()
