import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2
import sys
import os
# percorso assoluto della cartella exp1
path_exp1 = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "exp1"))
sys.path.append(path_exp1)

import script_lab2