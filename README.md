# XP_Extinction_Toolkit

[English](README.md) | [中文](README-zh.md)

`XP_Extinction_Toolkit` is a Python library designed for astronomy researchers. It provides a set of extinction correction tools for spectroscopic and photometric data, based on **the median extinction curve** calculated from about 370,000 high-quality samples in Gaia DR3 and LAMOST DR7. More information can be found in our article (Zhang et al., under review). 
The functions include:

1. **Function `ext_curve`**  
   Obtain the high-precision median extinction curve for the Milky Way.

2. **Function `deredden`**  
   Performs extinction correction on input spectra, suitable for the range of 2900 - 53414 Å.

3. **Function `Cal_E4455`**  
   A tool that uses empirical formulas to convert E(B-V) and E(B-V)_SFD into E(440-550).

4. **Function `star_ext`**  
    Calculates the precise **extinction** and **effective wavelength** for input **bands** of stars with given Teff, logg, and reddening. The input reddening can be E(440-550), E(B-V), or E(B-V)_SFD (derived from the SFD full-sky two-dimensional dust reddening map [Schlegel et al. (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract)). This function use the XP extinction curve (default) or the various Rv-dependent extinction curve models provided by the [`dust_extinction`](https://dust-extinction.readthedocs.io/en/stable/index.html) package. 
    Supports the input of additional filter passbands stored as `.dat` files from the [SVO Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps/index.php).

5. **Function `star_reddening`** 
    Calculates the precise **reddenting** for input **colors** of stars with given Teff, logg, and reddening, using the XP extinction curve (default) or other extinction models. The input reddening can be E(440-550), E(B-V), or E(B-V)_SFD.

6. **Function `model_extinction`** 
    Computes simulated extinction for given bands, using the XP extinction curve (default) or other extinction models. It returns a DataFrame of modeled extinction values across different Teff, E(B-V), logg for each band.

Functions 4, 5, and 6 are based on the BOSZ spectral library, filter passbands, and extinction curves for estimating extinction (or reddening) in photometric bands (or colors). The results may differ from empirical measurements. When calculating extinction (or reddening) for passbands of GALEX, PS1, SDSS, Gaia, 2MASS, and WISE survey, we recommend using another Python package based on empirical reddening measurement, [`extinction_coefficient`](https://github.com/vnohhf/extinction_coefficient), for extinction correction. This package is mostly valid in the E(B-V) range of 0-0.5 mag and the temperature range of 4000-10000 K.


# How to Install
### Using pip
~~~
# from PyPI (recommmand)
pip install XP_Extinction_Toolkit

# from the master trunk on the repository, considered developmental code
pip install git+https://github.com/vnohhf/XP_Extinction_Toolkit.git
~~~

### From source
extinction_coefficient can be installed from the source code after downloading it from the git repo (https://github.com/vnohhf/XP_Extinction_Toolkit/):
~~~
python setup.py install
~~~

# Usage Guidelines
## 1. To obtain the XP extinction curve and plot
~~~python
import numpy as np
import matplotlib.pyplot as plt
from ExtinctionToolkit import ExtinctionToolkit
# Initialize the ExtinctionToolkit
ExtTool = ExtinctionToolkit()
# Obtain the XP extinction curve
λ, k_x_55 = ExtTool.ext_curve()

# plot
fig,ax = plt.subplots(2,1, figsize=(5,6), tight_layout=True, dpi=300)
ax[0].plot(λ[λ<10200], k_x_55[λ<10200], lw=1)
ax[1].plot(1e4/λ, k_x_55, lw=1)
ax[0].set(title='XP extinction curve (optical part)', ylabel='E(λ-55)/E(44-55)', xlabel='wavelength [Å]')
ax[1].set(title='XP extinction curve (optical and infrared)', ylabel='E(λ-55)/E(44-55)', xlabel='wavelength number [1/μm]')
for axi in ax:
    axi.minorticks_on()
    axi.grid(True,ls=':',color='dimgray',lw=0.2,which='major',zorder=1)
~~~

## 2.

## 3. Transfer E(B-V) or E(B-V)_SFD to E(440-550)
~~~python
# (1) input E(B-V)
ebv = [0.1, 0.3, 0.5]
Teff = [4000, 6000, 8000]
E4455 = ExtTool.Cal_E4455(ebv=ebv,Teff=Teff)

# (2) input E(B-V)_SFD
sfdebv = 0.2
E4455 = ExtTool.Cal_E4455(sfdebv=sfdebv)
~~~

## 4. Calculate the extinction of 1,000 stars in four bands.
~~~python
import numpy as np
from ExtinctionToolkit import ExtinctionToolkit
ExtTool = ExtinctionToolkit()
Band = ['2MASS.Ks','PS1.i','SDSS.g','WISE.W23']
# Generate random values for Teff and E(440-550)
np.random.seed(1)
Teff = np.random.uniform(low=4000, high=10000, size=1000)
E4455 = np.random.uniform(low=0, high=2, size=1000)
# Calculate the extinction for the defined bands
Extinction = ExtTool.star_ext(Band, E4455=E4455, Teff=Teff)
~~~

## 5. Calculate the reddening for given Teff and sfd E(B-V)_SFD
~~~python
# Define E(B-V)_SFD and Teff (Here can use array, tuple, etc.)
sfdebv = [0.3, 0.7, 1, 1.3]
Teff = [6542, 6000, 3900, 7659]
# Define the color indices for reddening calculation
Color = ['Johnson.B-Johnson.V','GAIA3.G-GAIA3.Grp']
# Calculate the reddening
Reddening = ExtTool.star_reddening(Color, sfdebv=sfdebv, Teff=Teff)
~~~

## 4. Initialize ExtinctionToolkit with the new filter storage directory
~~~python
# Please replace '/path/to/new/filters/' with the actual absolute path where the filters are stored.
ExtTool = ExtinctionToolkit(filter_dir='/path/to/new/filters/')
Band = ['New_Band']
Extinction = ExtTool.star_ext(Band, E4455=0.1, Teff=5000)
~~~


## 6. Obtain the model extinction of XP curve or input extinction curves
~~~python
# (1) Initialize the ExtinctionToolkit with different models
ExtTool_XP = ExtinctionToolkit(model='XP')
Band = ['GALEX.NUV','PS1.z','SDSS.u','WISE.W4']
# Calculate the model extinction
Model_Ext, _ = ExtTool_XP.model_extinction(Band,cal_lamb=False,energy=False)

# (2) Set 'cal_lamb' to True to calculate the effective wavelength as well. 
# The unit of the return value are angstrom (Å).
Model_Ext, Model_λeff = ExtTool_XP.model_extinction(['GALEX.NUV','PS1.z','SDSS.u'],cal_lamb=True,energy=False)

# (3) Some of the project filters are Photon counter, e.g. Gaia, GALEX, PAN-STARRS, SDSS, 
# and some are Energy counter, e.g. 2MASS, WISE. 
# When using Energy counter type filters, set 'energy' to True to get the correct effective wavelength.
Model_Ext, Model_λeff = ExtTool_XP.model_extinction(['2MASS.H','WISE.W3'],cal_lamb=True,energy=True)

# (4) Use different extinction models from the dust_extinction package
from dust_extinction.parameter_averages import CCM89,F99,F19,G23
# Initialize ExtinctionToolkit with F19 model
ExtTool_F19 = ExtinctionToolkit(model=F19)
Model_Ext, Model_λeff = ExtTool_F19.model_extinction(['PS1.z','SDSS.u'],cal_lamb=True,energy=False)

# Use G23 model with a different Rv value (default is 3.1).
ExtTool_G23 = ExtinctionToolkit(model=G23, Rv=2.5)
Model_Ext, Model_λeff = ExtTool_G23.model_extinction(['2MASS.H','WISE.W4'],cal_lamb=True,energy=True)
~~~
