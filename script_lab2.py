import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2


def ott_dati(excel,foglio):
    df = pd.read_excel(excel, sheet_name=foglio)
    x, y, dx, dy = df.to_numpy(dtype= "float").T
    labels = list(df.columns)
    return x, y, dx, dy, labels
    
def retta(param, x):
    m , q = param
    return m*x +q

#---CHI2 CON ERRORI SULLA X----
def x_quadro(x, y, dx, dy, m ,q ):
    residui = (y-m*x-q)**2
    errore = dy**2+(m**2)*dx**2
    return np.sum(residui/errore)


#----UTILE PER OTTENERE SUBITO IL GRAFICO DEI REISDUI DI UN FIT LINEARE CON ERRORI SIA SU X CHE SU Y----
def residui_lin(x, y, dx, dy, m ,q):
    return (y-m*x-q)/np.sqrt(dy**2+(m**2)*dx**2)    

#----GRAFICO RESIDUI DI UN FIT CON ERRORE SOLO SU Y----
def residui_y(y, f_teo, dy):
    return (y-f_teo)/dy



#---MODALITA DI FIT LINEARE NEL CASO DI ERRORI SIA SU X CHE SU Y----
def fitODR(x, y , dx, dy, stima=(1,0)):
    modello= Model(retta)
    dati = RealData(x, y, sx=dx, sy=dy)
    odr = ODR(dati, modello, beta0=stima)
    out = odr.run()
    m, q = out.beta
    dm, dq = out.sd_beta
    cov = out.cov_beta
    return m, q, dm, dq, cov

#---MODALITA DI FIT LINEARE CON ERRORI SOLO SU Y-----   
def fit_lin(f, x, y,dy, stima=[1,0]):
    popt, pcov = curve_fit(f, x, y, stima, dy)
    m, q = popt
    dm,dq = np.sqrt(np.diag(pcov))
    return m , q, dm, dq

#---CHI2 CON SOLO ERRORI SU Y----
def chi_quadro(f_teo, y, dy):
    return np.sum((y-(f_teo))**2/dy**2)

#---TEST DI IPOTESI PER ACCETTARE O RIGETTARE IL FIT---
def valuta_fit(chi2_val, nu, alpha=0.05):
    chi2_ridotto = chi2_val / nu
    p_value = chi2.sf(chi2_val, nu)
    if p_value< alpha:
        return chi2_ridotto, False
    else:
        return chi2_ridotto, True
  
#---CALCOLO DELLA DEVIAZIONE STANDARD DI UNA DISTRUBUZIONE UNIFORME---
def dev_std_unif(dx):
    return dx/np.sqrt(12)

#--- TEST DI CONSISTENZA---
def consistenza(x1,x2,dx1,dx2):
    return (np.abs(x1-x2))/np.sqrt(dx1**2+dx2**2)

def media_pesata(x1,x2,dx1,dx2):
    valore = (x1/dx1**2+x2/dx2**2)/(1/dx1**2+1/dx2**2)
    errore = (1/dx1**2+1/dx2**2)**(-1/2)
    return valore, errore


if __name__ == "__main__":
    
    
    #IMPORTIAMO I DATI DA EXCEL
    i_m,v_m,di_m,dv_m , nomi= ott_dati("/mnt/c/Users/robyb/OneDrive - Università di Napoli Federico II/exp1_lab2.xlsx", 1)  
    
    i_v, v_v, di_v, dv_v, nomi = ott_dati("/mnt/c/Users/robyb/OneDrive - Università di Napoli Federico II/exp1_lab2.xlsx", 2)
    
    #EFFETTUIAMO LA STIMA DEI PARAMETRI
    R_m, q_m, dR_m, dq_m, cov_m= fitODR(i_m, v_m, dev_std_unif(di_m), dev_std_unif(dv_m))
    print(f"il parametro R nel fit a monte risulta essere:({R_m:.3f}±{dR_m:.3f})kΩ\\L'intercetta risulta essere:({q_m:.3f}±{dq_m:.3f})V")
    
    R_v, q_v, dR_v, dq_v, cov_v= fitODR(i_v, v_v, dev_std_unif(di_v), dev_std_unif(dv_v))
    print(f"il parametro R nel fit a monte risulta essere:({R_v:.3f}±{dR_v:.3f})kΩ\\L'intercetta risulta essere:({q_v:.3f}±{dq_v:.3f})V")
    
    #STIMIAMO IL CHI2 E LA BONTA DEL FIT E COEFFICIENTE DI CORRELAZIONE
    chi_monte = x_quadro(i_m, v_m, dev_std_unif(di_m), dev_std_unif(dv_m), R_m, q_m)
    ngl = np.size(i_m)-2
    chi_m_rid, validità = valuta_fit(chi_monte, ngl)
    print(f"il valore del chi quadro ridotto risulta essere: {chi_m_rid:.3f}")
    if validità:
        print("Fit accettabile")
    else:
        print("Fit non accettabile")
        
    
    chi_valle = x_quadro(i_v, v_v, dev_std_unif(di_v), dev_std_unif(dv_v), R_v, q_v)
    ngl = np.size(i_v)-2
    chi_v_rid, validità = valuta_fit(chi_valle, ngl)
    print(f"il valore del chi quadro ridotto risulta essere: {chi_v_rid:.3f}")
    if validità:
        print("Fit accettabile")
    else:
        print("Fit non accettabile")
        
    
    #---COEFFICIENTI DI CORRELAZIONE TRA I E V---
    rho_m = np.corrcoef(i_m, v_m)[0,1]
    rho_v = np.corrcoef(i_v, v_v)[0,1]

    print(f"Correlazione monte: {rho_m:.3f}")
    print(f"Correlazione valle: {rho_v:.3f}")

    # RESIDUI
        
    residui_monte = residui(i_m,v_m,dev_std_unif(di_m),dev_std_unif(dv_m),R_m, q_m)
        
    residui_valle = residui(i_v, v_v, dev_std_unif(di_v), dev_std_unif(dv_v), R_v, q_v)
    
    
    #STIMA DEL VALORE VERO DELLA RESISTENZA USATA
    
    rv = 40                 
    ra = 0.588
    r_vero_monte= (R_m*rv)/(rv-R_m)
    dr_v_monte= (rv**2/(rv-R_m)**2)*dR_m
    r_vero_valle = (R_v - ra)
    dr_v_valle = dR_v
    
    print(f"La resistenza R_x nella configurazione a monte vale: {r_vero_monte:.3f}±{dr_v_monte:.3f} ")
    print(f"La resistenza R_x nella configurazione a valle vale: {r_vero_valle:.3f}±{dr_v_valle:.3f} ")

    #CONSISTENZA
    
    l = consistenza(r_vero_monte, r_vero_valle, dr_v_monte, dr_v_valle)
    if l<2:
        print("i due valori sono compatibili")
    elif l>2 and l<3:
        print("la compatibilità è dubbia")
    else:
        print("valori non compatibili")
        
    # media pesata e errore pesato
    R_x , e_R_x = media_pesata(r_vero_monte, r_vero_valle, dr_v_monte, dr_v_valle)
    print(f"il valore medio della resistenza è: ({R_x:.3f} ± {e_R_x:.3f})kΩ")
    
    

    
    #GRAFICI
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=False)

    # --- GRAFICO FIT ---
    ax1.errorbar(
        i_m, v_m,
        xerr=di_m, yerr=dv_m,
        fmt="o", color="blue", ecolor="blue",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Dati (conf. a monte)"
        )
    x_fit_m = np.linspace(0, i_m.max(), 200)
    y_fit_m = retta((R_m, q_m), x_fit_m)
    ax1.plot(x_fit_m, y_fit_m, "-", color="blue", linewidth=1.5)
    
    ax1.errorbar(
        i_v, v_v,
        xerr=di_v, yerr=dv_v,
        fmt="o", color="red", ecolor="red",
        capsize=2, elinewidth=0.8, markersize=5, alpha=0.85,
        label="Dati (conf. a valle)"
        )
    x_fit_v = np.linspace(0, i_v.max(), 200)
    y_fit_v = retta((R_v, q_v), x_fit_v)
    ax1.plot(x_fit_v, y_fit_v, "-", color="red", linewidth=1.5)
    
    y_fit_vero = retta((R_x, 0), x_fit_v)
    
   # --- STIMA DELL'INTERCETTA VERA CON PENDENZA NOTA R_x ---

    # Calcolo media dei punti teorici
    x_mean = np.mean(x_fit_v)
    y_mean = np.mean(y_fit_vero)

    # Intercetta
    q = y_mean - R_x * x_mean

    # Errore dell'intercetta (formula analitica con pendenza nota)
    N = len(x_fit_v)
    sigma_y = np.std(y_fit_vero, ddof=1)
    sigma_x = np.std(x_fit_v, ddof=1)
    
    dq = np.sqrt((sigma_y**2)/N + (R_x**2)*(sigma_x**2)/N)
    
    print(f"Intercetta teorica q = ({q:.3f} ± {dq:.3f})V")


    
    ax1.plot(x_fit_v, y_fit_vero, color= "green", linewidth =1, label = "grafico R_x")
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)

    ax1.set_ylabel(nomi[1])
    ax1.set_title("Stime di R")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # --- GRAFICO RESIDUI ---
    ax2.axhline(0, color="black", linewidth=1)
    
    ax2.scatter(i_m, residui_monte, color="blue", label="Residui (monte)")
    ax2.scatter(i_v, residui_valle, color="red", label="Residui (valle)")
    
    ax2.set_ylim(-10,10)
    ax2.set_xlabel(nomi[0])
    ax2.set_ylabel("Residui normalizzati")
    ax2.set_title("Grafico dei residui")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.show()

    
    


