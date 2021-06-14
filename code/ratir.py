import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fluxtools import effective_wavelength

for filt in "riZYJH":
    filepath = f"../data/RATIR/{filt}.csv"
    data = pd.read_csv(filepath,header=None,names=["Wavelength","Transmission"])
    data["Wavelength"] *= 10
    data["Transmission"] *= 100
    print("%s: %.3f" % (filt,effective_wavelength(data,show_plot=True)))