import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2

from exp1.script_lab2 import *

excel="/mnt/c/Users/giova/OneDrive - Università di Napoli Federico II/Lab2exp2.xlsx"
(R,I,dI),nomi=ott_dati(excel)

(A,B), cov, (dA,dB)=fit(retta, R, I, dI)

R_star=-B/A

dR_star=np.sqrt((1/A**2)*(dB**2)+(((B**2)/(A**4))*(dA**2)))+2*(-1/A)*(B/A**2)*cov[0,1]

R_1=220
R_2=120
dR_1=1 #DA INSERIRE IN LAB
dR_2=1 #DA INSERIRE IN LAB

Rx=(R_2/R_1)*R_star

dRx=np.sqrt((((dR_1/R_1)**2)+((dR_2/R_2)**2))+((dR_star/R_star)**2))*R_x

#best fit
retta_best=retta((A,B),R)
#Test Chi Quadro
chi=chi_quadro(retta_best, I, dI)
nu=len(R)-1
chi_rid, esito_test=valuta_fit(chi,nu,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ridotto ha valore {chi_rid}")
else:
    print(f"Il fit non è buono, il chi quadro ridotto ha valore {chi_rid}")

#RESIDUI
residui=residui_lin(R,I,0,dI,A,B)

#GRAFICO RESIDUI
ax2.axhline(0,color="black",linewidth=1)
ax2.scatter(I, residui, color="r")
ax2.set_ylim(-10,10)
ax2.set_xlabel(nomi[0])
ax2.set_ylabel("Residui normalizzati")
ax2.set_title("Grafico dei residui")
ax2.grid(True, alpha=0.3)
ax2.legend()
plt.tight_layout()
plt.show()

#GRAFICI
grafico=plt.plot(R,I,color="b", label="R3 vs I")
grafico.errorbar(yerr=dI)
plt.grid()
plt.show()





