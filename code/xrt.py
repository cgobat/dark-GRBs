import pandas as pd
from selenium import webdriver
from selenium.webdriver.support.wait import WebDriverWait
from .uncertainty import AsymmetricUncertainty

def XRT_lightcurve(burst_id,lookuptable):
    """
    Function for retrieving X-ray observations (flux values) for a given gamma-ray burst.
    
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    burst_id : string
        GRB ID/name in the form YYMMDDx
    lookuptable : pandas DataFrame
        the reference table to get the TriggerNumber

    Returns
    -------
    fluxdata : pandas DataFrame
        table containing the time series fluxes in the XRT 0.3-10 keV band. Columns are Time, Tpos, Tneg, Flux, Fluxpos, and Fluxneg.

    Raises
    ------
    
    """
    trigger = lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"]
    lightcurveURL = f"https://www.swift.ac.uk/xrt_curves/{int(trigger):0>8}/"
    
    fireFoxOptions = webdriver.FirefoxOptions()
    fireFoxOptions.headless = True
    with webdriver.Firefox(options=fireFoxOptions) as browser:
        browser.get(lightcurveURL)
        browser.find_element_by_id('flux_makeDownload').click() # find the link to the data file and virtually click it
        WebDriverWait(browser, 30).until(lambda page: ".qdp" in page.current_url) # wait for the click to go through/data file to load
        lightcurveURL = browser.current_url # update the URL with the new page location (the actual data file)

    fluxdata = pd.read_table(lightcurveURL, header=1).apply(pd.to_numeric, errors="coerce").dropna().reset_index().apply(pd.to_numeric)
    fluxdata.columns = ["Time","Tpos","Tneg","Flux","Fluxpos","Fluxneg"]
    fluxdata["GRB"] = [burst_id]*len(fluxdata)
    #print("Retrieved",burst_id)
    return fluxdata


def get_photonIndex(burst_id,lookuptable):
    """
    Function for retrieving the X-ray spectral index for a given gamma-ray burst.
    
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    burst_id : string
        GRB ID/name in the form YYMMDDx
    lookuptable : pandas DataFrame
        the reference table to get the TriggerNumber

    Returns
    -------
    gamma : AsymmetricUncertainty
        value of the spectral index in the form (value (pos_err, neg_err))

    Raises
    ------
    
    """
    trigger = lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"]
    spectrumURL = f"https://www.swift.ac.uk/xrt_spectra/{int(trigger):0>8}/"
    
    spectra_tables = pd.read_html(spectrumURL)
    PC_table = spectra_tables[len(spectra_tables)-2]
    photon_index = PC_table.loc[PC_table[0]=="Photon index",1].values
    (Gamma, Gammapos, Gammaneg) = (float(num) for num in "".join([char for char in str(photon_index[0]) if char not in "[]()+-,"]).split())
    gamma = AsymmetricUncertainty(Gamma, Gammapos, Gammaneg)
    return gamma
