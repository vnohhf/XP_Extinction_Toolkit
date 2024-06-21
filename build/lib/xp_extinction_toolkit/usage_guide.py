'''
Author       : vnohhf
Date         : 2024-06-10 20:11
LastEditTime : 2024-06-21 22:34
E-mail       : zry@mail.bnu.edu.cn
Description  : Copyright© vnohhf. ALL RIGHTS RESERVED.
'''
#%%
import numpy as np
import matplotlib.pyplot as plt
from ExtinctionToolkit import ExtinctionToolkit

#%% To obtain the XP extinction curve and plot
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


#%% Extinction correction of the input spectrum
ExtTool = ExtinctionToolkit()
# Obtain the BOSZ spectrum. The wavelengths of the spectra have been modified to match the XP extinction curves.
BOSZ = np.load(ExtTool.start_path+'/BOSZ_spectra.npy',allow_pickle=1) 
# Get the wavelength
λ, _ = ExtTool.ext_curve()
# Take a random one as example, get the original spectrum and its Teff
ind = np.random.choice(range(len(BOSZ)))
original_spec = BOSZ[ind]
teff = ExtTool.teff_list[ind]
# Redden the original spectrum with E(B-V) of 0.1.
reddened_spec =  ExtTool.deredden(ebv=0.1, Teff=teff, wave=λ) * original_spec
# Deredden the reddened spectrum
# (Here intrinsic_spec is original_spec, showing only the usage of dereddening)
intrinsic_spec = 1 / ExtTool.deredden(ebv=0.1, Teff=teff, wave=λ) * reddened_spec

# plot
fig,ax = plt.subplots(1,1, figsize=(8,4), tight_layout=True, dpi=300)
ax.plot(λ, intrinsic_spec, lw=0.9, label='Intrinsic spectrum')
ax.plot(λ, reddened_spec, lw=0.9, label='Reddened spectrum')
ax.set(title=f'Teff = {teff} [K]', ylabel='Flux [$cm^{-2} s^{-1} Å^{-1}$]', xlabel='Wavelength [Å]')
ax.minorticks_on()
ax.legend()

# Reddening/de-reddening of multiple spectra
ebvlist = np.random.uniform(low=0, high=2, size=len(BOSZ))
reddened_spec_matrix = 1 / ExtTool.deredden(ebv=ebvlist, Teff=ExtTool.teff_list, wave=λ) * BOSZ

#%% Transfer E(B-V) or E(B-V)_SFD to E(440-550)
# (1) input E(B-V)
ebv = [0.1, 0.3, 0.5]
Teff = [4000, 6000, 8000]
E4455 = ExtTool.Cal_E4455(ebv=ebv,Teff=Teff)

# (2) input E(B-V)_SFD
sfdebv = 0.2
E4455 = ExtTool.Cal_E4455(sfdebv=sfdebv)


#%% Calculate the extinction of four bands for 1,000 stars
ExtTool = ExtinctionToolkit()
Band = ['2MASS.Ks','PS1.i','SDSS.g','WISE.W2']
# Generate random values for Teff and E(440-550)
np.random.seed(1)
Teff = np.random.uniform(low=4000, high=10000, size=1000)
E4455 = np.random.uniform(low=0, high=2, size=1000)
# Calculate the extinction for the defined bands
Extinction = ExtTool.star_ext(Band, E4455=E4455, Teff=Teff)


#%% Calculate the reddening for stars with given Teff and E(B-V)
ExtTool = ExtinctionToolkit()
# Define E(B-V)_SFD and Teff (Here can use array, tuple, etc.)
sfdebv = np.array([0.3, 0.7, 1, 1.3])
Teff = [6542, 6000, 3900, 7659]
# Define the color indices for reddening calculation
Color = ['Johnson.B-Johnson.V','GAIA3.G-GAIA3.Grp']
# Calculate the reddening
Reddening = ExtTool.star_reddening(Color, sfdebv=sfdebv, Teff=Teff)

#%% Using external filters
# Initialize ExtinctionToolkit with the new filter storage directory
# Please replace '/path/to/new/filters/' with the actual absolute path where the filters are stored.
ExtTool = ExtinctionToolkit(filter_dir='/path/to/new/filters/')
Band = ['New_Band']
Extinction = ExtTool.star_ext(Band, E4455=0.1, Teff=5000)


#%% Calculation of simulated extinction at different SEDs
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
