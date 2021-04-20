import pandas as pd
from selenium import webdriver
from selenium.webdriver.support.wait import WebDriverWait

def XRT_lightcurve(burst_id,lookuptable):
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
    print("Retrieved",burst_id)
    return fluxdata


def get_BetaX(burst_id,lookuptable):
    trigger = lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"]
    spectrumURL = f"https://www.swift.ac.uk/xrt_spectra/{int(trigger):0>8}/"
    
    spectra_tables = pd.read_html(spectrumURL)
    PC_table = spectra_tables[len(spectra_tables)-2]
    photon_index = PC_table.loc[PC_table[0]=="Photon index",1].values
    (Gamma, Gammapos, Gammaneg) = (float(num) for num in "".join([char for char in str(photon_index[0]) if char not in "[]()+-,"]).split())
    
    return Gamma, Gammapos, Gammaneg
