# Optical darkness in short-duration $\gamma$-ray bursts

[![arXiv:2304.09122](https://img.shields.io/badge/arXiv-2304.09122-b31b1b)](https://arxiv.org/abs/2304.09122)

<b><a href="https://github.com/cgobat">Caden Gobat</a> <a href="https://orcid.org/0000-0003-1268-8845"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" height=14px/></a>,<sup>1,2</sup> <a href="https://github.com/ajvanderhorst">Alexander van der Horst</a> <a href="https://orcid.org/0000-0001-9149-6707"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" height=14px/></a>,<sup>1</sup> and <a href="https://github.com/djfitz3999">David Fitzpatrick</a><sup>3,4</sup></b></br>
<sup>1</sup> Department of Physics, George Washington University, 725 21st St NW, Washington, DC 20052, U.S.A.</br> 
<sup>2</sup> Department of Space Operations, Southwest Research Institute, 1050 Walnut Street, Suite 300, Boulder, CO 80302, U.S.A.</br>
<sup>3</sup> Department of Aerospace Engineering Sciences, University of Colorado Boulder, 3775 Discovery Dr, Boulder, CO 80303, U.S.A.</br>
<sup>4</sup> Department of Physics, Georgetown University, 37th \& O St NW, Washington, DC 20007, U.S.A.

---

Accepted for publication in [*Montly Notices of the Royal Astronomical Society*](https://academic.oup.com/mnras) on April 14, 2023.

**Abstract:** Gamma-ray bursts categorically produce broadband afterglow emission, but in some cases, emission in the optical band is dimmer than expected based on the contemporaneously observed X-ray flux. This phenomenon, aptly dubbed "optical darkness", has been studied extensively in long GRBs (associated with the explosive deaths of massive stars), with possible explanations ranging from host environment extinction to high redshift to possibly unique emission mechanisms. However, investigations into optical darkness in short GRBs (associated with the mergers of compact object binaries) have thus far been limited. This work implements a procedure for determining the darkness of GRBs based on spectral indices calculated using temporally-matched *Swift*-XRT data and optical follow-up observations; presents a complete and up-to-date catalog of known short GRBs that exhibit optical darkness; and outlines some of the possible explanations for optically dark short GRBs. In the process of this analysis, we developed versatile and scalable data processing code that facilitates reproducibility and reuse of our pipeline. These analysis tools and resulting complete sample of dark short GRBs enable a systematic statistical study of the phenomenon and its origins, and reveal that optical darkness is indeed quite rare in short GRBs, and highly dependent on observing response time and observational effects.

## Overview

Most of the front-facing code is housed in Jupyter notebooks in this top-level directory, notably [`analysis.ipynb`](./analysis.ipynb) and [`pipeline.ipynb`](./pipeline.ipynb) (`analysis` calls `pipeline` to do most of its heavy lifting). Supporting code, scripts, and utility functions reside in the [./src](./src/) directory.

Steps to get up and running:
1. Clone this repository using Git, or just download it as a .zip and extract it.
2. Install necessary python packages using `pip install -r requirements.txt` from a command line running within this repository folder.
3. Install the [`asymmetric_uncertainty` package](https://github.com/cgobat/asymmetric_uncertainty) by following its [installation instructions](https://github.com/cgobat/asymmetric_uncertainty#installation).
4. Run the [`analysis.ipynb`](./analysis.ipynb) Jupyter notebook, which in turn calls the [`pipeline.ipynb`](./pipeline.ipynb) notebook to compile data from the various sources, process it for analysis, and perform temporal matching and calculation of $\beta_\text{ox}$. Most plotting/visualization work is done in [`analysis.ipynb`](./analysis.ipynb).

Alternatively, GitHub will render Jupyter notebooks, so they can also just be viewed/inspected here directly.

The [`./src/xrt.py`](./src/xrt.py) module mostly contains functions for querying the [UKSSDC](https://www.swift.ac.uk/index.php) to retrieve *Swift* X-Ray Telescope data, incuding afterglow lightcurves, spectral parameters, temporal behavior, and related information like galactic column densities ($N_H$).

### Legacy code

This work has heritage in the research done by [David Fitzpatrick](https://github.com/djfitz3999) for his [bachelor's thesis (2020)](./pub/Fitzpatrick%20thesis%202020.pdf). The following tools were originally developed for that work and are no longer used in this codebase, but are included for posterity.

- [`calc_beta_ox.py`](./src/legacy/Calculation%20Code/calc_beta_ox.py) and [`calc_beta_ox.cpp`](./src/legacy/Calculation%20Code/calc_beta_ox.cpp) (compiled with Cygwin on Windows at [`Automating the Calculation of Beta_OX.exe`](./src/legacy/Calculation%20Code/Automating%20the%20Calculation%20of%20Beta_OX.exe), or recompile it yourself) both load files from [`./data/legacy/`](./data/legacy)

- [`Graphing_Beta_OX.py`](./src/legacy/Graphing%20Code/Graphing_Beta_OX.py) loads the files generated by one of [the two aforementioned scripts](./src/legacy/Calculation%20Code) from [`./products/Generated Files (C++)/`](./products/Generated%20Files%20(C%2B%2B))

- Both tools have been updated to include GUI-based filesystem interaction and configuration.

Track progress on this project's [board](https://github.com/cgobat/dark-GRBs/projects/1). (no longer current)
