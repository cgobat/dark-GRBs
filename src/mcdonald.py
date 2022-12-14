import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fluxtools import effective_wavelength

for filt in "ugriz":
    filepath = f"../data/McDonald/{filt}.csv"
    data = pd.read_csv(filepath)
    print("%s: %.3f" % (filt,effective_wavelength(data,show_plot=True,dl=10)))