import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

(V_r, dVr, V_i, dV_i, dt, Tr, f_i, Dt, DT), nomi= ott_dati("/mnt/c/Users/Robyb/OneDrive - Università di Napoli Federico II/exp4.xlsx")

#DATI DEL LAB
R= 99.30#da inserire domani
Cap= 216.6*1e-9#da inserire domani
L=9.942*1e-3
R_l= 220

f_i= f_i*1000
Tr = Tr * 1e-6
DT= DT*1e-6
dt=dt*1e-6
Dt=Dt*1e-6

f_ur=(2*np.pi)/Tr
df_ur =(2*np.pi*DT)/(Tr**2)

A_r=(V_r/V_i)
dA_r = np.sqrt(
    (dVr/V_i)**2 +
    (V_r*dV_i/V_i**2)**2
)

phi_r=(-(f_ur)*dt)

dphi_r = np.sqrt(
    (dt*df_ur)**2 +
    (f_ur*Dt)**2
)

f_t_teo = 3.3*1000 # DA DETERMINARE IN LAB
n_punti = 10

indice_centrale = np.argmin(np.abs(f_i - f_t_teo))

inizio = max(0, indice_centrale - n_punti // 2)
fine = min(len(f_i), indice_centrale + (n_punti + 1) // 2)

maschera = np.zeros_like(f_i, dtype=bool)
maschera[inizio+1:fine] = True

A_rfit=A_r[maschera]
f_ifit=f_i[maschera]
phi_rfit=phi_r[maschera]

(A,B,C),cov,(d_A,d_B,d_C)=fit(parabola,f_ifit,A_rfit,dA_r[maschera])
m,q,dm,dq,cov_p= fit_lin(f_ifit,phi_rfit,dphi_r[maschera])

parab = parabola((A,B,C),f_ifit)
retta_phi= retta((m,q),f_ifit)

chi_1 = chi_quadro(parab, A_rfit, dA_r[maschera])
ngl = np.size(A_rfit)-3
chi_rid, esito_test=valuta_fit(chi_1,ngl,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ha valore {chi_1} con gradi di libertà:{ngl}")
else:
    print(f"Il fit non è buono, il chi quadro ha valore {chi_1}, con gradi di liberà :{ngl}")

chi_2 = chi_quadro(retta_phi, phi_rfit, dphi_r[maschera])
ngl = np.size(phi_rfit)-2
chi_rid, esito_test=valuta_fit(chi_2,ngl,alpha=0.05)
if esito_test:
    print(f"Il fit è buono, il chi quadro ha valore {chi_2} con gradi di libertà:{ngl}")
else:
    print(f"Il fit non è buono, il chi quadro ha valore {chi_2}, con gradi di liberà :{ngl}")

print(f"I parametri del fit parabolico sono: ({A}+-{d_A});({B}+-({d_B});({C:.3f}+-{d_C:.3f}))")
print(f"I parametri del fit rettilineo sono:({m}+-{dm});({q:.3f}+-{dq:.3f})")

A_m = C - B**2/(4*A)
w_0_Att = -B/(2*A)
Q_Att = B/(np.sqrt(2*(B**2 - 4*A*C)))

Q_phi=q/2
w_0_phi=-q/m

d_Am_quadro = (
    (B**4/(16*A**4))*d_A**2
    + (B**2/(4*A**2))*d_B**2
    + d_C**2
    - (B**3/(4*A**3))*cov[0,1]
    + (B**2/(2*A**2))*cov[0,2]
    - (B/A)*cov[1,2]
)

d_Am = np.sqrt(d_Am_quadro)


d_Am=np.sqrt(d_Am_quadro)

d_Q_quadro_Att = (1 / (B**2 - 4 * A * C)**3) * (
    2 * B**2 * A**2 * d_C**2 
    + 8 * A**2 * C**2 * d_B**2 
    + 2 * B**2 * C**2 * d_A**2 
    + 4 * C * B**2 * A * cov[2,0] 
    - 8 * C * B * A**2 * cov[2,1] 
    - 8 * C**2 * B * A * cov[1,0]
)

d_Q_Att=(B)/(np.sqrt(2*(B**2-4*A*C)))

d_w0_Att_quadro = (
    (B**2/(4*A**4))*d_A**2
    + (1/(4*A**2))*d_B**2
    - (B/(2*A**3))*cov[0,1]
)

d_w0_Att = np.sqrt(d_w0_Att_quadro)

d_w0_Att=np.sqrt(d_w0_Att_quadro)

d_w0_quadro_phi = (
    (1 / m**2) * dq**2 
    + (q**2 / m**4) * dm**2 
    + 2 * (-1 / m) * (q / m**2) * cov_p
)

d_w0_phi=np.sqrt(d_w0_quadro_phi)

d_Q_quadro_phi = (dq**2) / 4

d_Q_phi=np.sqrt(d_Q_quadro_phi)

l1 = consistenza(w_0_Att, w_0_phi, d_w0_Att, d_w0_phi)
if l1<2:
    print("i due valori della frequenza di risonanza sono compatibili")
elif l1>2 and l1<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")

l1 = consistenza(Q_Att, Q_phi, d_Q_Att, d_Q_phi)
if l1<2:
    print("i due valori del guadagno sono compatibili")
elif l1>2 and l1<3:
    print("la compatibilità è dubbia")
else:
    print("valori non compatibili")

pesi=np.array([peso(d_w0_Att),peso(d_w0_phi)])
x= np.array([w_0_Att,w_0_phi])

f_t = np.sum(pesi * x) / np.sum(pesi)

# errore sulla media pesata
df_t = np.sqrt(1 / np.sum(pesi))

pesi_q =np.array([peso(d_Q_Att),peso(d_Q_phi)])
x_q= np.array([Q_Att,Q_phi])

Q_t = np.sum(pesi_q * x_q) / np.sum(pesi_q)

# errore sulla media pesata
dQ_t = np.sqrt(1 / np.sum(pesi_q))

R_l_vero = R/A_m -R
dR_l_vero = (R/A_m**2)*d_Am

print(f"l'attenuazione vale{A_m}+-{d_Am}, il guadagno:{Q_Att}+- {d_Q_Att}, la pulsazione: {w_0_Att}+-{d_w0_Att}")
print(f"il guadagno nel fit rettilineo vale: {Q_phi}+-{d_Q_phi}; la pulsazione nel fit rettilineo: {w_0_phi}+-{d_w0_phi}")
print(f"Il valore del guadagno pesato è: ({Q_t:.3f}+-{dQ_t:.3f})")
print(f"Il valore della frequenza di risonananza pesata è: ({f_t}+-{df_t})")
print(f"la resistenza parassita è:({R_l_vero}+-{dR_l_vero})")

def risonanza(R,L,C,w):
    return R/(np.sqrt(R**2+(L*w-1/(C*w))**2))

def phi_ris(R,L,C,w):
    return -np.arctan((L*w-1/(C*w))/R)


figure, (ax1,ax2) = plt.subplots(2,1, figsize=(8,10), sharex=False)

x_parab = np.logspace(np.log10(f_ifit.min()),np.log10(f_ifit.max()),200)
x_retta = np.logspace(np.log10(f_ifit.min()),np.log10(f_ifit.max()),200)
x_fit_A = np.logspace(np.log10(f_i.min()),np.log10(f_i.max()),200)
y_fit = risonanza(R,L,C,x_fit_A)
y_fit_phi = phi_ris(R,L,C,x_fit_A)
y_parab = parabola((A,B,C),x_parab)
y_retta=retta((m,q),x_retta)

ax1.errorbar(f_i/f_t_teo, A_r
         ,yerr=dA_r,
        fmt="o", color="blue", ecolor="black",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Punti sperimentali della attenuazione"
        )
        

ax1.scatter(f_i/f_t_teo, A_r,  color="blue", linewidth=1.5 )
ax1.plot(x_parab/f_t_teo, y_parab, "k-", linewidth=2, label= "Parabola nell'intorno di f_t")

ax2.errorbar(f_i/f_t_teo, phi_r
         ,yerr= dphi_r,
        fmt="o", color="red", ecolor="black",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Punti sperimentali dello sfasamento"
        )
ax2.scatter(f_i/f_t_teo, phi_r, color="red", linewidth=1.5 )
ax2.plot(x_retta/f_t_teo, y_retta, "k-", linewidth=2, label= "Retta nell'intorno di f_t")

ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_yscale('log')

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






