import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
from exp1.script_lab2 import *

#----CREAZIONE DEL PLOT DELL'ESPERIMENTO----
figure, (ax1,ax2) = plt.subplots(2,1, figsize(8,10), sharex=False)
print("ciao")