# Optical darkness in short-duration &#x1D6FE;-ray bursts

[![doi:10.1093/mnras/stad1189](https://img.shields.io/badge/DOI-10.1093%2Fmnras%2Fstad1189-informational)](https://doi.org/10.1093/mnras/stad1189)&ensp;[![arXiv:2304.09122](https://img.shields.io/badge/arXiv-2304.09122-b31b1b)](https://arxiv.org/abs/2304.09122)

<b>[Caden Gobat](https://github.com/cgobat) [<img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" height=14px/>](https://orcid.org/0000-0003-1268-8845),<sup>1,2</sup> [Alexander van der Horst](https://github.com/ajvanderhorst) [<img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" height=14px/>](https://orcid.org/0000-0001-9149-6707),<sup>1</sup> and [David Fitzpatrick](https://github.com/djfitz3999)<sup>3,4</sup></b></br>
<sup>1</sup> Department of Physics, George Washington University, 725 21st St NW, Washington, DC 20052, U.S.A.</br> 
<sup>2</sup> Department of Space Operations, Southwest Research Institute, 1050 Walnut Street, Suite 300, Boulder, CO 80302, U.S.A.</br>
<sup>3</sup> Department of Aerospace Engineering Sciences, University of Colorado Boulder, 3775 Discovery Dr, Boulder, CO 80303, U.S.A.</br>
<sup>4</sup> Department of Physics, Georgetown University, 37th \& O St NW, Washington, DC 20007, U.S.A.

---

Published in [*Montly Notices of the Royal Astronomical Society*](https://academic.oup.com/mnras/article/523/1/775/7136159) on April 21, 2023.

<details>
<summary><b>Abstract</b></summary>
<p>
Gamma-ray bursts (GRBs) categorically produce broad-band afterglow emission, but in some cases, emission in the optical band is dimmer than expected based on the contemporaneously observed X-ray flux. This phenomenon, aptly dubbed "optical darkness", has been studied extensively in long GRBs (associated with the explosive deaths of massive stars), with possible explanations ranging from host environment extinction to high redshift to possibly unique emission mechanisms. However, investigations into optical darkness in short GRBs (associated with the mergers of compact object binaries) have thus far been limited. This work implements a procedure for determining the darkness of GRBs based on spectral indices calculated using temporally-matched *Swift* X-ray Telescope data and optical follow-up observations; presents a complete and up-to-date catalogue of known short GRBs that exhibit optical darkness; and outlines some of the possible explanations for optically dark short GRBs. In the process of this analysis, we developed versatile and scalable data processing code that facilitates reproducibility and reuse of our pipeline. These analysis tools and resulting complete sample of dark short GRBs enable a systematic statistical study of the phenomenon and its origins, and reveal that optical darkness is indeed quite rare in short GRBs, and highly dependent on observing response time and observational effects.
</p>
</details>

---

## Data and results

The optical data compiled for and used in this work can be found in [`products/all_optical.csv`](./products/all_optical.csv). This catalog is a compilation of data from [`data/newData.xlsx`](./data/newData.xlsx) (compiled for this work from GCNs, misc. publications, etc.),  [`data/OpticalData.csv`](./data/OpticalData.csv) (from [Fong *et al.* 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...815..102F)) and [`data/Rastinejad_Table1.csv`](./data/Rastinejad_Table1.csv) (from [Rastinejad *et al.* 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...916...89R)).

Our X-ray data is primarily sourced from the [UK Swift Science Data Centre](https://www.swift.ac.uk), and time-resolved flux data can be found compiled in [`products/Swift_XRT_lightcurves.csv`](./products/Swift_XRT_lightcurves.csv). We also retrieve X-ray spectral indices ($\beta_\text{x} = \Gamma - 1$), and include these in our main catalog of short GRBs ([`products/Swift_sGRB_catalog.csv`](./products/Swift_sGRB_catalog.csv)).  We supplement these data from the UKSSDC with Fong *et al.* (2015)'s X-ray data, which can be found in [`data/XRayData.csv`](./data/XRayData.csv) (fluxes) and [`data/BetaXData.csv`](./data/BetaXData.csv) (X-ray spectral indices).

[`products/TableA1.csv`](./products/TableA1.csv) is the digital version of Table A1 in the paper. For each GRB in our sample, it contains a listing of how many optical and X-ray data points we had available for that burst (from the data described above), how many temporal matches we were able to make using those data, and how many of those temporal matches qualify as optically dark. For any GRB with one or more optically dark data points, a plot of its optical/X-ray lightcurve can be found in [`products/dark lightcurves/`](./products/dark%20lightcurves/).

## Code Overview

Most of the front-facing code is housed in Jupyter notebooks in this top-level directory, notably [`analysis.ipynb`](./analysis.ipynb) and [`pipeline.ipynb`](./pipeline.ipynb) (`analysis` calls `pipeline` to do most of its heavy lifting). Supporting code, scripts, and utility functions reside in the [`src`](./src/) directory.

Steps to get up and running:
1. Clone this repository using Git, or just download it as a .zip and extract it.
2. Install necessary python packages using `pip install -r requirements.txt` from a command line running within this repository folder.
3. Install the [`asymmetric_uncertainty` package](https://github.com/cgobat/asymmetric_uncertainty) by following its [installation instructions](https://github.com/cgobat/asymmetric_uncertainty#installation).
4. Run the [`analysis.ipynb`](./analysis.ipynb) Jupyter notebook, which in turn calls the [`pipeline.ipynb`](./pipeline.ipynb) notebook to compile data from the various sources, process it for analysis, and perform temporal matching and calculation of $\beta_\text{ox}$. Most plotting/visualization work is done in [`analysis.ipynb`](./analysis.ipynb).

Alternatively, GitHub will render Jupyter notebooks, so they can also just be viewed/inspected here directly.

The [`src/xrt.py`](./src/xrt.py) module mostly contains functions for querying the [UKSSDC](https://www.swift.ac.uk/index.php) to retrieve *Swift* X-Ray Telescope data, incuding afterglow lightcurves, spectral parameters, temporal behavior, and related information like galactic column densities ($N_H$).

### Legacy code

This work has heritage in the research done by [David Fitzpatrick](https://github.com/djfitz3999) for his [bachelor's thesis (2020)](./pub/Fitzpatrick%20thesis%202020.pdf). The following tools were originally developed for that work and are no longer used in this codebase, but are included for posterity.

- [`src/legacy/Calculation Code/calc_beta_ox.py`](./src/legacy/Calculation%20Code/calc_beta_ox.py) and [`src/legacy/Calculation Code/calc_beta_ox.cpp`](./src/legacy/Calculation%20Code/calc_beta_ox.cpp) (compiled with Cygwin on Windows at [`src/legacy/Calculation Code/Automating the Calculation of Beta_OX.exe`](./src/legacy/Calculation%20Code/Automating%20the%20Calculation%20of%20Beta_OX.exe), or recompile it yourself) both load files from [`data/legacy/`](./data/legacy)

- [`src/legacy/Graphing Code/Graphing_Beta_OX.py`](./src/legacy/Graphing%20Code/Graphing_Beta_OX.py) loads the files generated by one of [the two aforementioned scripts](./src/legacy/Calculation%20Code) from [`products/Generated Files (C++)/`](./products/Generated%20Files%20(C%2B%2B)/)

Both tools have been updated to include GUI-based filesystem interaction and configuration.
