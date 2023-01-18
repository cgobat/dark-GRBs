import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate, integrate

def effective_wavelength(filter_response, show_plot=False): # pass a dataframe with columns Wavelength (in Ang), Transmission (in %)
    filter_response.sort_values(by="Wavelength",inplace=True)

    vega_spec = pd.read_table("http://svo2.cab.inta-csic.es/svo/theory/fps3/morefiles/vega.dat",
                              delimiter=" ",header=None,names=["Wavelength","Flux"])
    vega_func = interpolate.interp1d(vega_spec["Wavelength"],vega_spec["Flux"],
                                     bounds_error=False,fill_value=0)
    λ = filter_response["Wavelength"]
    T = filter_response["Transmission"]
    Vg = vega_func(λ)

    numerator = integrate.trapz(y=λ*T*Vg, x=λ) # ∫ λ T(λ) Vg(λ) dλ
    denominator = integrate.trapz(y=T*Vg, x=λ) # ∫ T(λ) Vg(λ) dλ
    
    lambda_eff = numerator/denominator
    
    if show_plot:
        fig,ax = plt.subplots()
        ax.plot(vega_spec.Wavelength,vega_spec.Flux,"b",label="Vega spectrum")
        ax.set_yscale("log")
        ax.set_xlim(filter_response["Wavelength"].min(), filter_response["Wavelength"].max())
        ax.set_ylabel("Flux [erg/cm$^2$/s/Ang]")
        ax.set_xlabel(r"$\lambda$ [Ang]")
        ax2 = ax.twinx()
        ax2.plot(λ,T,color="dimgrey",linestyle=":",label="Filter response")
        ax2.set_ylabel("Transmission [%]")
        ax2.axvline(lambda_eff,color="g",linestyle="--",label=r"$\lambda_\mathrm{eff}$")
        ax2.set_ylim(0,100)
        fig.legend(loc="upper right")
        plt.show()

    return lambda_eff

def lightcurve(grb, band="optical", xlimits=False, ylimits=False, **kwargs):
    """Intended for use in an environment where DataFrames called `xrt_data` and `all_optical`
    already exist and contain afterglow optical/X-ray flux data over time, respectively. This
    function uses that information to plot the light curve of the burst in a band of the user's
    choice."""

    global xrt_data
    global all_optical
    if band == "xray":
        subset = xrt_data.loc[xrt_data["GRB"]==grb]
        neg_err = [0.4*flux.value if np.isinf(flux.minus) else flux.minus for flux in subset["SpecFlux"]]
        pos_err = [0.4*flux.value if np.isinf(flux.plus) else flux.plus for flux in subset["SpecFlux"]]
        plt.errorbar(subset.Time,[flux.value for flux in subset.SpecFlux],
                     xerr=np.array(subset.Tneg,subset.Tpos).T, yerr=np.array((neg_err,pos_err)),
                     linestyle="", capthick=0, **kwargs)
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
        ax.errorbar(subset["Time (s)"],[flux.value for flux in subset["Flux (Jy)"]], marker=".", linestyle="", capthick=0,
                    yerr=np.array((neg_err, [point.plus for point in subset["Flux (Jy)"]])),
                    uplims=[np.isinf(point.minus) for point in subset["Flux (Jy)"]],
                    lolims=[np.isinf(point.plus) for point in subset["Flux (Jy)"]], **kwargs)
        ax.grid(linestyle="--")
        ax.set(xscale="log",yscale="log",xlabel="Time (s)",ylabel=f"Flux ({band}) [Jy]")
        if xlimits:
            ax.set(xlim=xlimits)
        if ylimits:
            ax.set(ylim=ylimits)
    plt.show()

def simulate_spectrum(idx):
    """Intended for use in an environment where a DataFrame called `results` already
    exists and contains matched optical/X-ray data points. This function uses that
    information to plot the broadband "spectrum" at a selected point."""

    global results
    inter_freqs = np.linspace(results.loc[idx,"nu_o"],0.3/4.135667696e-18,100)
    xray_freqs = np.linspace(0.3,10,100)/4.135667696e-18
    
    ox_spec = inter_freqs**(-results.loc[idx,"B_ox"].value)
    ox_spec *= (results.loc[idx,"F_o"]/ox_spec[0])
    x_spec = xray_freqs**(-results.loc[idx,"B_x"])
    x_spec *= (ox_spec[-1]/x_spec[0])
    plt.plot(inter_freqs,ox_spec,label=r"$\beta_{ox}=%f$"%(-results.loc[idx,"B_ox"]))
    plt.plot(xray_freqs,x_spec,label=r"$\beta_x=%f$"%(-results.loc[idx,"B_x"]))
    plt.vlines([10**np.mean((np.log10(0.3),np.log10(10))) / 4.135667696e-18],0,
                10**np.log10((x_spec.max()+x_spec.min())/2),"k",linestyle="--")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    #plt.gca().set_yticklabels([])
    plt.show()