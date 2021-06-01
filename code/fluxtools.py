import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate

def effective_wavelength(filter_response,show_plot=False): # pass a dataframe with columns Wavelength (in Ang), Transmission (in %)
    vega_spec = pd.read_table("http://svo2.cab.inta-csic.es/svo/theory/fps3/morefiles/vega.dat",
                              delimiter=" ",header=None,names=["Wavelength","Flux"])
    vega_function = interpolate.interp1d(vega_spec["Wavelength"],vega_spec["Flux"])
    response_function = interpolate.interp1d(filter_response["Wavelength"],filter_response["Transmission"])
    
    dl = 0.1 # Angstrom
    
    domain = np.arange(filter_response["Wavelength"].min(),filter_response["Wavelength"].max(),dl)
    numerator = np.sum([domain*vega_function(domain)*response_function(domain)*dl])
    denominator = np.sum([vega_function(domain)*response_function(domain)*dl])
    lambda_eff = numerator/denominator
    if show_plot:
        plt.plot(vega_spec.Wavelength,vega_spec.Flux,"b")
        plt.yscale("log")
        plt.ylabel("Vega flux [erg/cm$^2$/s/Ang]")
        plt.xlabel(r"$\lambda$ [Ang]")
        plt.twinx()
        plt.plot(filter_response.Wavelength,filter_response.Transmission,color="dimgrey",linestyle=":")
        #plt.xscale("log")
        plt.ylabel("Filter transmission [%]")
        plt.vlines([lambda_eff],0,100,color="g",linestyle="--")
        plt.xlim(domain.min(), domain.max())
        plt.ylim(0,100)
        plt.show()

    return lambda_eff

def lightcurve(grb,band="optical",xlimits=False,ylimits=False,**kwargs):
    if band == "xray":
        subset = xrt_data.loc[xrt_data["GRB"]==grb]
        neg_err = [0.4*flux.value if np.isinf(flux.minus) else flux.minus for flux in subset["SpecFlux"]]
        pos_err = [0.4*flux.value if np.isinf(flux.plus) else flux.plus for flux in subset["SpecFlux"]]
        plt.errorbar(subset.Time,[flux.value for flux in subset.SpecFlux],
                     xerr=np.array(subset.Tneg,subset.Tpos).T,yerr=np.array((neg_err,pos_err)),
                     linestyle="",capthick=0,**kwargs)
        plt.xscale("log")
        plt.yscale("log")
        plt.title(f"Swift XRT Lightcurve for GRB {grb}")
        plt.xlabel("Time since trigger [s]")
        plt.ylabel("Flux (~1 keV) [Jy]")
        plt.grid(linestyle="--")
        if xlimits:
            plt.xlim(xlimits)
        if ylimits:
            plt.ylim(ylimits)
    else:
        if band == "optical":
            subset = all_optical.loc[(all_optical["GRB"]==grb) & (all_optical["λ_eff"]>=3000) & (all_optical["λ_eff"]<=8000)]
        elif band == "UV":
            subset = all_optical.loc[(all_optical["GRB"]==grb) & (all_optical["λ_eff"]<3000)]
        elif band == "IR":
            subset = all_optical.loc[(all_optical["GRB"]==grb) & (all_optical["λ_eff"]>8000)]
        else:
            print("Unrecognized band selected")
            return None
        
        fig,ax = plt.subplots()
        neg_err = [0.4*flux.value if np.isinf(flux.minus) else flux.minus for flux in subset["Flux (Jy)"]]
        pos_err = [0.4*flux.value if np.isinf(flux.plus) else flux.plus for flux in subset["Flux (Jy)"]]
        ax.errorbar(subset["Time (s)"],[flux.value for flux in subset["Flux (Jy)"]],marker=".",linestyle="",capthick=0,
                    yerr=np.array((neg_err,[point.plus for point in subset["Flux (Jy)"]])),
                    uplims=[np.isinf(point.minus) for point in subset["Flux (Jy)"]],lolims=[np.isinf(point.plus) for point in subset["Flux (Jy)"]],**kwargs)
        ax.grid(linestyle="--")
        ax.set(xscale="log",yscale="log",xlabel="Time (s)",ylabel=f"Flux ({band}) [Jy]")
        if xlimits:
            ax.set(xlim=xlimits)
        if ylimits:
            ax.set(ylim=ylimits)
    plt.show()

def simulate_spectrum(idx):
    inter_freqs = np.linspace(results.loc[idx,"nu_o"],0.3/4.135667696e-18,100)
    xray_freqs = np.linspace(0.3,10,100)/4.135667696e-18
    
    ox_spec = inter_freqs**(-results.loc[idx,"B_ox"].value)
    ox_spec = ox_spec * (results.loc[idx,"F_o"]/ox_spec[0])
    x_spec = xray_freqs**(-results.loc[idx,"B_x"])
    x_spec = x_spec * (ox_spec[-1]/x_spec[0])
    plt.plot(inter_freqs,ox_spec,label=r"$\beta_{ox}=%f$"%(-results.loc[idx,"B_ox"]))
    plt.plot(xray_freqs,x_spec,label=r"$\beta_x=%f$"%(-results.loc[idx,"B_x"]))
    plt.vlines([10**np.mean((np.log10(0.3),np.log10(10))) / 4.135667696e-18],0,10**np.log10((x_spec.max()+x_spec.min())/2),"k",linestyle="--")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    #plt.gca().set_yticklabels([])
    plt.show()