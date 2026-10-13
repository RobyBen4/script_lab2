import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2

from script_lab2 import *

excel="/mnt/c/Users/giova/OneDrive - Università di Napoli Federico II/Lab2exp2.xlsx"
(R,I,dI),nomi=ott_dati(excel)

(A,B), cov, (dA,dB)=fit(retta, R, I, dI)

