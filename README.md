# XP_Extinction_Toolkit

[English](README.md) | [中文](README-zh.md)

`XP_Extinction_Toolkit` is a Python library designed for astronomy researchers. It provides a set of extinction correction tools for spectroscopic and photometric data, based on **the median extinction curve** calculated from about 370,000 high-quality samples in Gaia DR3 and LAMOST DR7. More information can be found in our article (Zhang et al., publishing). 
The functions are as follows. See the [Usage Guide](#Usage_Guide) for examples of use.

1. ### **Function `ext_curve`**
   Obtain the high-precision median extinction curve for the Milky Way.

2. ### **Function `deredden`**  
   Performs extinction correction on input spectra, suitable for the range of 2900 - 53414 Å.

3. ### **Function `Cal_E4455`**  
   A tool that uses empirical formulas to convert E(B-V) and E(B-V)_SFD into E(440-550).

4. ### **Function `star_ext`**  
    Calculates the precise **extinction** and **effective wavelength** for input **bands** of stars with given Teff, logg, and reddening. The input reddening can be E(440-550), E(B-V), or E(B-V)_SFD (derived from the SFD full-sky two-dimensional dust reddening map [Schlegel et al. (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract)). This function use the XP extinction curve (default) or the various Rv-dependent extinction curve models provided by the [`dust_extinction`](https://dust-extinction.readthedocs.io/en/stable/index.html) package. 
    Supports the input of additional filter passbands stored as `.dat` files from the [SVO Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps/index.php).

5. ### **Function `star_reddening`**
    Calculates the precise **reddenting** for input **colors** of stars with given Teff, logg, and reddening, using the XP extinction curve (default) or other extinction models. The input reddening can be E(440-550), E(B-V), or E(B-V)_SFD.

6. ### **Function `model_extinction`**
    Computes simulated extinction for given bands, using the XP extinction curve (default) or other extinction models. It returns a DataFrame of modeled extinction values across different Teff, E(B-V), logg for each band.

Functions 4, 5, and 6 are based on the BOSZ spectral library, filter passbands, and extinction curves for estimating extinction (or reddening) in photometric bands (or colors). The results may differ from empirical measurements. When calculating extinction (or reddening) for passbands of GALEX, PS1, SDSS, Gaia, 2MASS, and WISE survey, we recommend using another Python package based on empirical reddening measurement, [`extinction_coefficient`](https://github.com/vnohhf/extinction_coefficient), for extinction correction. This package is mostly valid in the E(B-V) range of 0 - 0.5 mag and the temperature range of 4000 - 10000 K.


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

# Usage Guide
## 1. To obtain the XP extinction curve and plot
~~~python
import numpy as np
import matplotlib.pyplot as plt
from xp_extinction_toolkit import ExtinctionToolkit

# Initialize the ExtinctionToolkit
ExtTool = ExtinctionToolkit()
# Obtain the XP extinction curve
λ, k_x_55 = ExtTool.ext_curve()

# plot
fig,ax = plt.subplots(2,1, figsize=(5,6), tight_layout=True, dpi=300)
ax[0].plot(λ[λ<10200], k_x_55[λ<10200], lw=1)
ax[1].plot(1e4/λ, k_x_55, lw=1)
ax[0].set(title='XP extinction curve (optical part)', ylabel='E(λ-55)/E(44-55)', xlabel='Wavelength [Å]')
ax[1].set(title='XP extinction curve (optical and infrared)', ylabel='E(λ-55)/E(44-55)', xlabel='Wavelength number [1/μm]')
for axi in ax:
    axi.minorticks_on()
    axi.grid(True,ls=':',color='dimgray',lw=0.2,which='major',zorder=1)
~~~

## 2. Extinction correction of input spectra
### 2.1 single spectrum
~~~python
ExtTool = ExtinctionToolkit()
# Obtain the BOSZ spectrum
BOSZ = np.load(ExtTool.start_path+'/BOSZ_spectra.npy',allow_pickle=1) 

# Get the wavelength, which have been modified to match the XP extinction curves.
λ, _ = ExtTool.ext_curve()

# Take a random one as example, get the original spectrum and its Teff
ind = np.random.choice(range(len(BOSZ)))
original_spec = BOSZ[ind]
teff = ExtTool.teff_list[ind] 

# Redden the original spectrum with E(B-V) of 0.1.
# Input ebv conversion to E4455 requires temperature information
reddened_spec = original_spec / ExtTool.deredden(ebv=0.1, Teff=teff, wave=λ)

# Deredden the reddened spectrum
# (Here intrinsic_spec is equal to original_spec, showing only the usage of dereddening)
intrinsic_spec = reddened_spec * ExtTool.deredden(ebv=0.1, Teff=teff, wave=λ) #

# plot
fig,ax = plt.subplots(1,1, figsize=(8,4), tight_layout=True, dpi=300)
ax.plot(λ, intrinsic_spec, lw=0.9, label='Intrinsic spectrum')
ax.plot(λ, reddened_spec, lw=0.9, label='Reddened spectrum')
ax.set(title=f'Teff = {teff} [K]', ylabel='Flux [$cm^{-2} s^{-1} Å^{-1}$]', xlabel='Wavelength [Å]')
ax.minorticks_on()
ax.legend()
~~~

### 2.2 Reddening/de-reddening of multiple spectra
~~~python
ebvlist = np.random.uniform(low=0, high=2, size=len(BOSZ))
reddened_spec_matrix = 1 / ExtTool.deredden(ebv=ebvlist, Teff=ExtTool.teff_list, wave=λ) * BOSZ
~~~

## 3. Transfer E(B-V) or E(B-V)_SFD to E(440-550)
### 3.1 input E(B-V)
~~~python
from xp_extinction_toolkit import ExtinctionToolkit
ebv = [0.1, 0.3, 0.5]
Teff = [4000, 6000, 8000]
E4455 = ExtTool.Cal_E4455(ebv=ebv,Teff=Teff)
~~~

### 3.2 input E(B-V)_SFD
~~~python
sfdebv = 0.2
E4455 = ExtTool.Cal_E4455(sfdebv=sfdebv)
~~~

## 4. Calculate the extinction of four bands for 1,000 stars
~~~python
import numpy as np
from xp_extinction_toolkit import ExtinctionToolkit

ExtTool = ExtinctionToolkit()

# Check bulit-in passbands
print(ExtTool.Filters.keys())

# Generate random values for Teff and E(440-550)
Band = ['2MASS.Ks', 'PS1.i', 'SDSS.g', 'WISE.W2']
Teff = np.random.uniform(low=4000, high=10000, size=1000)
E4455 = np.random.uniform(low=0, high=2, size=1000)

# Calculate the extinction for the defined bands
Extinction = ExtTool.star_ext(Band, E4455=E4455, Teff=Teff)
~~~

## 5. Calculate the reddening for stars with given Teff and E(B-V)
### 5.1 Using built-in filters
We have built-in transmission for GCPD_Johnson, GALEX, PS1, SDSS, Gaia3, 2MASS, and WISE.
~~~python
# Define E(B-V)_SFD and Teff (Here can use array, tuple, etc.)
sfdebv = [0.3, 0.7, 1, 1.3]
Teff = [6542, 6000, 3900, 7659]

# Define the color indices
Color = ['Johnson.B-Johnson.V','GAIA3.G-GAIA3.Grp']

# Calculate the reddening
Reddening = ExtTool.star_reddening(Color, sfdebv=sfdebv, Teff=Teff)
~~~

### 5.2 Using external filters
Please replace '/path/to/new/filters/' with the actual absolute path where the filters are stored.
~~~python
ExtTool = ExtinctionToolkit(filter_dir='/path/to/new/filters/')
Band = ['New_Band']
Extinction = ExtTool.star_ext(Band, E4455=0.1, Teff=5000)
~~~


## 6. Calculation of simulated extinction at different SEDs
### 6.1 Using XP extinction curve
~~~python
ExtTool_XP = ExtinctionToolkit(model='XP')
Band = ['GALEX.NUV','PS1.z','SDSS.u','WISE.W4']

# Calculate the simulated extinction
Model_Ext, _ = ExtTool_XP.model_extinction(Band,cal_lamb=False,energy=False)
~~~

### 6.2 Obtain the effective wavelength
Set 'cal_lamb' to True to calculate the effective wavelength (unit: angstrom (Å)).

~~~python
# Calculate the simulated extinction and corresponding effective wavelength
Model_Ext, Model_λeff = ExtTool_XP.model_extinction(['GALEX.NUV','PS1.z','SDSS.u'],cal_lamb=True,energy=False)
~~~

Some of the project filters are Photon counter (e.g. Gaia, GALEX, PAN-STARRS, SDSS) and some are Energy counter (e.g. 2MASS, WISE). 
When using Energy counter type filters, set 'energy' to True to get the correct effective wavelength.

~~~python
Model_Ext, Model_λeff = ExtTool_XP.model_extinction(['2MASS.H','WISE.W3'],cal_lamb=True,energy=True)
~~~

### 6.3 Using Rv-dependent models provided in [`dust_extinction`](https://dust-extinction.readthedocs.io/en/stable/index.html)
~~~python
from dust_extinction.parameter_averages import CCM89,F99,F19,G23

# Use F19 extinction curve
ExtTool_F19 = ExtinctionToolkit(model=F19)
Model_Ext, Model_λeff = ExtTool_F19.model_extinction(['PS1.z','SDSS.u'],cal_lamb=True,energy=False)

# Use G23 model with a different Rv value (default is 3.1).
ExtTool_G23 = ExtinctionToolkit(model=G23, Rv=2.5)
Model_Ext, Model_λeff = ExtTool_G23.model_extinction(['2MASS.H','WISE.W4'],cal_lamb=True,energy=True)
~~~
