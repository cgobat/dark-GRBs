import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup as bs
from astropy.coordinates import SkyCoord
from astropy.table import Table
# from selenium import webdriver
# from selenium.webdriver.support.wait import WebDriverWait
from .uncertainty import AsymmetricUncertainty

grb_list = pd.read_table("https://www.swift.ac.uk/xrt_curves/grb.list",
                         sep=" |\t",header=None,engine="python",
                         names=["_","GRB","TriggerNumber"]).drop("_",axis=1)

def XRT_lightcurve(burst_id,lookuptable=grb_list):
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
    IndexError
        if the table could not be retrieved
    
    """
    trigger = lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"]
#     lightcurveURL = f"https://www.swift.ac.uk/xrt_curves/{int(trigger):0>8}/"
    
#     fireFoxOptions = webdriver.FirefoxOptions()
#     fireFoxOptions.headless = True
#     with webdriver.Firefox(options=fireFoxOptions) as browser:
#         browser.get(lightcurveURL)
#         browser.find_element_by_id('flux_makeDownload').click() # find the link to the data file and virtually click it
#         WebDriverWait(browser, 30).until(lambda page: ".qdp" in page.current_url) # wait for the click to go through/data file to load
#         lightcurveURL = browser.current_url # update the URL with the new page location (the actual data file)
    lightcurveURL = f"https://www.swift.ac.uk/xrt_curves/{int(trigger):0>8}/flux_incbad.qdp"
    fluxdata = pd.DataFrame(columns=['Time', 'Time_perr', 'Time_nerr', 'Flux', 'Flux_perr', 'Flux_nerr'])
    i = 0
    while True:
        try:
            current = Table.read(lightcurveURL,format="ascii.qdp",table_id=i,names=["Time","Flux"]).to_pandas()
            fluxdata = fluxdata.append(current,ignore_index=True)
            i += 1
        except:
            if i>0:
                break
            else:
                raise IndexError
    fluxdata.columns = ["Time","Tpos","Tneg","Flux","Fluxpos","Fluxneg"]
    fluxdata = fluxdata.apply(pd.to_numeric)
    fluxdata["GRB"] = [burst_id]*len(fluxdata)
    #print("Retrieved",burst_id)
    return fluxdata


def get_columnDensity(burst_id,lookuptable=grb_list):
    """
    Function for retrieving the neutral hydrogen column density in the direction of a given gamma-ray burst.
    
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    burst_id : string
        GRB ID/name in the form YYMMDDx
    lookuptable : pandas DataFrame
        the reference table to get the TriggerNumber

    Returns
    -------
    N_H : AsymmetricUncertainty
        value of the neutral hydrogen column density in the form (value (pos_err, neg_err)) [units of cm^-2]

    Raises
    ------
    
    """
    try:
        trigger = int(lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"])
    except ValueError:
        trigger = int(grb_list.loc[grb_list["GRB"] == burst_id, "TriggerNumber"])
    spectrumURL = f"https://www.swift.ac.uk/xrt_spectra/{trigger:0>8}/"
    
    page = requests.get(spectrumURL)
    soup = bs(page.content,"html.parser")
    assert page.status_code != 404, "404 Error."
    assert len(soup.findAll("table",{"summary":"Model fitted to interval0 data"}))>0, "No tables found."

    tables = {}
    for section in soup.findAll("div",{"class":"fitres"}):
        mode = section.find("h3").text[:2]
        current = section.find("table",{"summary":"Model fitted to interval0 data"})
        if current is not None:
            tables[mode] = current
    if "PC" in tables.keys():
        table = tables["PC"] # use PC mode info if available
        used_mode = "PC"
    else:
        table = tables["WT"] # if not, use WT mode and flag as such
        used_mode = "WT"

    for entry in table.children:
        label = entry.find("th")
        if isinstance(label,int):
            pass
        elif label.text == "NH (intrinsic)":
            NH_str = entry.find("td").text
            
    mult = NH_str.index("×")
    power = NH_str[mult+4:mult+6]
    vals = NH_str[:mult-1]+NH_str[mult+6:]
    nominal, plus, minus = vals.split()[0], vals[vals.index("+")+1:vals.index(",")], vals[vals.index("-")+1:vals.index(")")]
    N_H = AsymmetricUncertainty(float(nominal)*10**int(power), float(plus)*10**int(power), float(minus)*10**int(power))
    
    return N_H, used_mode


def get_photonIndex(burst_id,lookuptable=grb_list):
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
    try:
        trigger = int(lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"])
    except ValueError:
        trigger = int(grb_list.loc[grb_list["GRB"] == burst_id, "TriggerNumber"])
    spectrumURL = f"https://www.swift.ac.uk/xrt_spectra/{trigger:0>8}/"
    
    page = requests.get(spectrumURL)
    soup = bs(page.content,"html.parser")
    assert page.status_code != 404, "404 Error."
    assert len(soup.findAll("table",{"summary":"Model fitted to interval0 data"}))>0, "No tables found."

    tables = {}
    for section in soup.findAll("div",{"class":"fitres"}):
        mode = section.find("h3").text[:2]
        current = section.find("table",{"summary":"Model fitted to interval0 data"})
        if current is not None:
            tables[mode] = current
    if "PC" in tables.keys():
        table = tables["PC"] # use PC mode info if available
        used_mode = "PC"
    else:
        table = tables["WT"] # if not, use WT mode and flag as such
        used_mode = "WT"

    for entry in table.children:
        label = entry.find("th")
        if isinstance(label,int):
            pass
        elif label.text == "Photon index":
            photon_index = entry.find("td").text
            
    (Gamma, Gammapos, Gammaneg) = (float(num) for num in "".join([char for char in photon_index if char not in "[]()+-,"]).split())
    gamma = AsymmetricUncertainty(Gamma, Gammapos, Gammaneg)
    return gamma, used_mode

def get_temporalIndex(burst_id,query_time,lookuptable=grb_list):
    """
    Function for retrieving the temporal index (power-law slope) of a gamma-ray burst at a given time.
    
    Author: Caden Gobat, George Washington University

    Parameters
    ----------
    burst_id : string
        GRB ID/name in the form YYMMDDx
    query_time : numeric
        time (in seconds) at which to retrieve the temporal index
    lookuptable : pandas DataFrame
        the reference table to get the TriggerNumber

    Returns
    -------
    alpha : AsymmetricUncertainty
        value of the temporal index at the specified time in the form (value (pos_err, neg_err)) [units of cm^-2]

    Raises
    ------
    
    """
    trigger = lookuptable.loc[lookuptable["GRB"] == burst_id, "TriggerNumber"]
    livecatURL = f"https://www.swift.ac.uk/xrt_live_cat/{int(trigger):0>8}/"
    
    livecat_tables = pd.read_html(livecatURL)
    slopes_table = livecat_tables[2]
    assert len(slopes_table.columns)==2
    
    breaktimes = [slopes_table.iloc[i,1][:-2] for i in slopes_table.index if "T" in slopes_table.iloc[i,0]]
    for i,time in enumerate(breaktimes):
        if "×" in time:
            vals,power = time.split('×')
            power = power.split()[0][2:]
            nominal, plus, minus = vals.split()[0], vals[vals.index("+")+1:vals.index(",")], vals[vals.index("-")+1:vals.index(")")]
            time = AsymmetricUncertainty(float(nominal), float(plus), float(minus)) * 10**int(power)
        else:
            time = AsymmetricUncertainty(time)
        breaktimes[i] = time
    breaktimes += [np.inf]
    
    alphas = [AsymmetricUncertainty(slopes_table.iloc[i,1]) for i in slopes_table.index if "α" in slopes_table.iloc[i,0]]
    
    assert len(breaktimes)+len(alphas) == len(slopes_table.index)+1, "Mismatch error in parsing table rows"
    
    for i,time in enumerate(breaktimes):
        if query_time < time:
            alpha = alphas[i]
            break
        
    return alpha