'''
Author       : vnohhf
Date         : 2024-05-20 12:22
LastEditTime : 2024-06-21 22:30
E-mail       : zry@mail.bnu.edu.cn
Description  : Copyright© vnohhf. ALL RIGHTS RESERVED.
'''
#%%
import os
import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy import units as u
from astropy import constants as const
from scipy.interpolate import interp1d,interp2d

class ExtinctionToolkit:
    def __init__(self, filter_dir='', model='XP', Rv=3.1):
        """
        Initializes the ExtinctionToolkit with the specified model, filter directory, and Rv (R55) value.

        Parameters:
            filter_dir (str): The directory where filter data files are stored.
            model (str): The extinction model to use ('XP' by default).
            Rv (float): The R value for the extinction curve.
        """
        self.start_path = os.path.dirname(os.path.abspath(__file__))
        # The extinction model in use
        self.model_func = model
        self.Rv = Rv
        # XP spectral extinction curve
        self.Ecurve = {}
        self.λ, self.Ecurve['k(λ-55)'] = self.ext_curve()
        self.λ_intvl = np.diff(self.λ)
        self.λ_intvl = np.append(self.λ_intvl,self.λ_intvl[-1])
        # Transform the form of the XP_curve
        self.XPEC_Rv = 3.073
        self.XPEC_R55 = 2.730
        self.Ecurve['A(λ)/A(55)'] = self.Ecurve['k(λ-55)'] / self.XPEC_R55 + 1
        # Filter data
        if filter_dir != '':
            self.user_filter_dir = filter_dir
        self.Filters = self.load_filters()
        # Model parameter grid
        self.logg_list, self.teff_list, self.EBV_list = self.para_grid()
        self.E4455_list = self.EBV_list

    #% Check whether the entered value is valid
    def Check_and_Expand(self, inputLi): 
        """
        Validates the input list dimensions and expands single values into arrays of consistent length.

        Parameters:
            inputLi (list, array): The input list to validate and possibly expand.

        Returns:
            np.ndarray: Expanded and validated list of inputs.
        """
        inputLi2 = np.array(inputLi, dtype=object)
        lengthLi = np.empty_like(inputLi2, dtype=np.float64)
        midLi = np.empty_like(inputLi2)
        outputLi = np.empty_like(inputLi2)
        
        if len(inputLi2.shape) > 1:
            raise ValueError("Input must be a 1-dimensional array or list")
    
        for i, input in enumerate(inputLi2):
            midLi[i] = np.atleast_1d(input).astype(np.float64)
            lengthLi[i] = len(midLi[i])

        # If the length of the input items, is not 1, then they need to be consistent.
        not1item = lengthLi[lengthLi!=1]
        if (len(not1item) >= 2) & (len(set(not1item)) > 1):
            raise ValueError('Input sequences must be the same size')
        
        # Expanding single values into arrays of equal length, except for the first one.
        if np.max(lengthLi) > 1:
            for i, (l, item) in enumerate(zip(lengthLi, midLi)):
                if l != 1:
                    outputLi[i] = item
                else:
                    outputLi[i] = np.full(int(np.max(lengthLi)), item[0])
            return outputLi
        else:
            return midLi
    

    #% Provide XP extinction curves
    def ext_curve(self, *wave):
        """
        Provides the extinction curves using the XP model.

        Parameters:
            wave (list, array): Optional wavelengths at which to evaluate the extinction curves.

        Returns:
            list: Default wavelengths and extinction curve at default wavelengths or input wavelengths (if have).
        """
        XPEC = pd.read_csv(self.start_path+'/XP_extinction_curve.csv')
        λ = np.array(XPEC.wave[XPEC.wave<5341.4])*10 # Wavelengths (unit:Å)
        k_x_55 = np.array(XPEC.ext)[XPEC.wave<5341.4]
        if wave == ():
            return λ, k_x_55
        else:
            return λ, interp1d(λ, k_x_55)(wave)


    #% Deredden the input spectrum
    def deredden(self, sfdebv=[], ebv=[], E4455=[], Teff=[], wave=[]):
        """
        Deredden the input spectrum based on the XP extinction curve.

        Parameters:
            E4455 (list, array): Extinction at 4455 Å.
            wave (tuple): Wavelengths for which dereddening is required (unit: Å).

        Returns:
            ndarray: Dereddened spectrum values.
        """
        if len(np.atleast_1d(E4455)) == 0:
            E4455 = self.Cal_E4455(sfdebv=sfdebv, ebv=ebv, Teff=Teff)
        else:
            E4455 = np.atleast_1d(E4455).astype(np.float64)
        if (np.min(wave) < self.λ[0]) | (np.max(wave) > 53414):
            raise ValueError('Input wavelength over range: 2900 - 53514 Å')
        # Transforming the form of XP_curve, interpolate and deredden (XP_curve的形式转换，插值，去红化)
        if len(wave) == 0:
            if len(E4455)==1:
                return 10**(-0.4 * self.XPEC_R55 * E4455 * self.Ecurve['A(λ)/A(55)'])
            else: # return a matrix
                return 10**(-0.4 * self.XPEC_R55 * E4455[:, np.newaxis] * self.Ecurve['A(λ)/A(55)'])
        else:
            f_interp = interp1d(self.λ, self.Ecurve['A(λ)/A(55)'])
            if len(E4455)==1:
                return 10**(0.4 * self.XPEC_R55 * E4455 * f_interp(wave).T)
            else: # return a matrix
                return 10**(0.4 * self.XPEC_R55 * E4455[:, np.newaxis] * f_interp(wave).T)
            


    #% Reading Filter Data 
    def load_filters(self):
        """
        Load passband data from the specified filter directory.

        Returns:
            DataFrame: A DataFrame containing processed filter data, which has been differenced to self.λ.
        """
        #* 整理名称等信息 - Organizing name and other information --------------------------------------------------------------------
        filter_dir = self.start_path+'/filters/'
        filelist = glob.glob(filter_dir + '*.dat')
        if hasattr(self,'user_filter_dir'):
            filelist = np.append(filelist, glob.glob(self.user_filter_dir + '*.dat'))
        Filters = pd.DataFrame()

        #* 读取 - Reading --------------------------------------------------------------------
        for i, fpath in enumerate(filelist):
            # 名称信息 - Name information
            fname = fpath.split('/')[-1]
            _, survey = fname.split('.')[0].split('_')
            band = fname.split('.')[1]
            # 通带 - Bandpass
            table = pd.read_csv(fpath, header=None, delimiter=" ")
            throughput_wave = table[0].to_numpy()
            temp_tran = table[1].to_numpy()
            # 特殊波段 - Special bands
            if (survey == 'GALEX') and (band in ['FUV', 'NUV']):
                temp_tran = temp_tran / (np.pi * 25**2)
            # 储存 - Storing
            Filters[survey + '.' + band] = np.interp(self.λ, throughput_wave, temp_tran, left=0, right=0)

        #* 添加440、550nm的窄通带数据 - Adding data for narrow bands at 440, 550nm -----------------------------------------------
        for i, cw, B in zip(len(Filters) + np.array([0, 1]), [4400, 5500], ['440', '550']):
            nrrw_lamb = np.linspace(cw - 80, cw + 80, 20)  # unit:ang
            nrrw_tran = np.ones_like(nrrw_lamb)
            Filters[B] = np.interp(self.λ, nrrw_lamb, nrrw_tran, left=0, right=0)

        return Filters


    #% Establish the model parameter grid
    def para_grid(self):
        """
        Establishes a grid of model parameters for teff, logg, and EBV values.

        Returns:
            list: Lists of logg, teff, and EBV values.
        """
        EBV_list = np.hstack(([0.001, 0.005], np.around(np.arange(0.01, 0.2, 0.02), decimals=2),
                              np.around(np.arange(0.2, 1, 0.05), decimals=2),
                              np.around(np.arange(1, 3, 0.1), decimals=2)))
        logg_list, teff_list = np.vstack((np.array([[g, t] for t in np.arange(3750, 6001, 250) for g in np.arange(0, 3, 0.5)]),
                                          np.array([[g, t] for t in np.arange(3750, 7500, 250) for g in np.arange(3, 5.1, 0.5)]),
                                          np.array([[g, t] for t in np.arange(7500, 20000, 500) for g in np.arange(3, 5.1, 0.5)]))).T
        return logg_list, teff_list, EBV_list


    #% Calculate the simulated extinction for the input band 
    # 计算输入波段的模拟消光
    def model_extinction(self,Band,cal_lamb=False,energy=False):
        """
        Calculates simulated extinction for input bands using XP extinction curve (default) or other input extinction models.

        Parameters:
            Band (str): The photometric band for which extinction is calculated. Format: photometric system + '.' + Passband name.
            cal_lamb (bool): Whether to calculate effective wavelength.
            energy (bool): Whether to calculate using energy instead of photon counts.

        Returns:
            list of DataFrame: Extinction values and effective wavelengths.
        """
        # BOSZ spectra, wavelength range consistent with XP_extinction_curve (unit:Å)
        self.BOSZ = np.load(self.start_path+'/BOSZ_spectra.npy',allow_pickle=1) 

        #* Determine zero-magnitude stars in the AB system by the number of photons through each filter
        #* 确定AB系统0等星，通过各滤光片的光子数 -> photon_0m_ab 
        photon_0m_ab = {}
        light_speed = const.c.value # type: ignore
        freqLi = light_speed*1e10/self.λ
        f0ν = 10**(48.6/-2.5) * np.ones(len(freqLi))   #AB system zero magnitude flux (unit: erg cm-2 s-1 Hz-1)
        f0λ = light_speed*1e10/self.λ**2 * f0ν   # (unit: erg cm-2 s-1 ang-1)
        for B in Band:
            if B in self.Filters.keys():
                photon_0m_ab[B] = sum( f0λ * self.λ_intvl *1e-7 / (const.h.value * light_speed*1e10/self.λ) * self.Filters[B] ) # type: ignore
            else:
                raise ValueError(f"The input filter '{B}' does not exist in the presets. Please specify its directory with the 'filter_dir' field.")

        #* Calculate extinction and corresponding effective wavelength for each band
        #* 计算各波段消光和对应有效波长 -> A_model,λ_eff ---------------------------------
        multiindex = pd.MultiIndex.from_arrays([self.teff_list,self.logg_list], names=['teff','logg'])
        A_model,λ_eff = [pd.DataFrame(index=multiindex, columns=Band).astype(np.float32) for _ in range(2)]
        A_model = A_model.applymap(lambda x: np.empty(len(self.EBV_list))) # type: ignore
        λ_eff   = λ_eff  .applymap(lambda x: np.empty(len(self.EBV_list))) # type: ignore
        for g,t,flux in zip(self.logg_list,self.teff_list,self.BOSZ):
            # Reddening spectra -> extincLi
            if self.model_func == 'XP':
                extincLi = np.array([flux * 10**(-0.4* self.XPEC_R55*E * self.Ecurve['A(λ)/A(55)']) for E in np.append(0,self.E4455_list)])
            else:
                extincLi = np.array([flux * self.model_func(Rv=self.Rv).extinguish(1e4/self.λ /u.micron, Ebv=ebv)  for ebv in np.append(0,self.EBV_list)]) # type: ignore
            for B in Band:
                photon = np.sum( extincLi * self.λ_intvl *1e-7 / (const.h.value * light_speed * 1e10/self.λ) * np.array(self.Filters[B]) , axis=1) # type: ignore
                # Calculate extinction for each band -> A_model (unit: mag)
                magLi = np.log10(photon/photon_0m_ab[B])/-0.4
                A_model[B][t,g] = magLi[1:] - magLi[0]
                # Calculate effective wavelength for each band -> λ_eff (unit: Å)
                if cal_lamb:
                    if energy:
                        energy_counter = np.array(self.Filters[B])
                    else: # If the input filter in photon counter, convert to energy counter.
                        energy_counter = self.λ * np.array(self.Filters[B])
                    λ_eff[B][t,g] = np.sum( self.λ * extincLi[1:] * energy_counter * self.λ_intvl, axis=1) / np.sum( extincLi[1:] * energy_counter * self.λ_intvl, axis=1)
        return A_model,λ_eff


    #% Transfer E(B-V) or E(B-V)_SFD to E(440-550)
    def Cal_E4455(self,sfdebv=[],ebv=[],Teff=[]):
        """
        Transforms E(B-V) to E(440-550) using polynomial coefficients based on Teff and E(B-V); 
        or transforms E(B-V)_SFD to E(440-550) by an empirical funciton.

        Parameters:
            sfdebv (list, array): E(B-V) values from SFD map to transform.
            ebv (list, array): Standard E(B-V) values to transform.
            Teff (list, array): Effective temperatures corresponding to each E(B-V) value.

        Returns:
            array: Calculated E(440-550) values.
        """
        EBV = np.atleast_1d(ebv).astype(np.float64)
        sfdEBV = np.atleast_1d(sfdebv).astype(np.float64)
        TEFF = np.atleast_1d(Teff).astype(np.float64)
        if len(EBV) != 0:
            if len(TEFF) == 0:
                raise ValueError('Need to enter Teff when transfer E(B-V)')
            e4455_ebv_coefficient = np.array(
            [[3.82467764e+03, 4.65580735e-03, 1.11698239e+00],
             [4.08076994e+03, 3.60913231e-03, 1.12497787e+00],
             [4.33438645e+03, 2.81428447e-03, 1.12902479e+00],
             [4.58757770e+03, 4.50136493e-03, 1.12681132e+00],
             [4.84263674e+03, 6.50734691e-03, 1.12192864e+00],
             [5.09507736e+03, 7.43614474e-03, 1.11563331e+00],
             [5.34859998e+03, 8.58562934e-03, 1.10751863e+00],
             [5.59908991e+03, 9.83053797e-03, 1.09840152e+00],
             [5.85440666e+03, 1.09358805e-02, 1.08969143e+00],
             [6.11403494e+03, 1.02493077e-02, 1.08408783e+00],
             [6.36621983e+03, 1.15420271e-02, 1.07692928e+00],
             [6.61538283e+03, 1.11477107e-02, 1.07287315e+00],
             [6.87523478e+03, 1.12571069e-02, 1.06889168e+00],
             [7.50551511e+03, 1.07156778e-02, 1.06519772e+00],
             [8.69368952e+03, 9.31582914e-03, 1.06126303e+00],
             [1.00988273e+04, 9.93438541e-03, 1.05267772e+00],
             [1.15021637e+04, 1.05846406e-02, 1.04563394e+00],
             [1.28971461e+04, 1.05749522e-02, 1.04202309e+00],
             [1.42992234e+04, 1.06019595e-02, 1.03975732e+00]])
            # Find the Teff closest to the grid point
            E4455 = np.empty_like(EBV)
            for i,(e,t) in enumerate(zip(EBV,TEFF)):
                ind = np.abs(e4455_ebv_coefficient[:,0] - t).argmin()
                E4455[i] = np.polyval(e4455_ebv_coefficient[ind,1:], e) * e
        elif len(sfdEBV) != 0:
            if len(np.where(sfdEBV>1.5)[0]) > 0:
                print('Warning: E(B-V)_SFD greater than 1.5 will be set to 1.5.')
            sfdEBV[sfdEBV>1.5] = 1.5
            E4455 = 3.03018594*sfdEBV + -3.6425772*np.log(1+np.exp(sfdEBV))+2.50136596
        else:
            raise ValueError('Need to enter one of E(B-V) or E(B-V)_SFD')
        return E4455


    #% Calculate the accurate extinction for given Teff and E(440-550) 
    # 计算给定 Teff 和 E(440-550) 时准确的消光
    def star_ext(self,Band,ebv=[],sfdebv=[],E4455=[],Teff=None,Logg=4.5,cal_lamb=False,energy=False):
        """
        Calculates the precise extinction and effective wavelength for input bands of stars with given Teff, logg, and reddening, using the XP extinction curve (default) or other extinction models. The input reddening can be E(440-550), E(B-V), or E(B-V)_SFD, derived from the SFD full-sky two-dimensional dust reddening map by Schlegel et al. (1998).
        
        Parameters:
            Band (list of str): The photometric band for which extinction is calculated.
            ebv (list, array): List of E(B-V) values (unit: mag).
            sfdebv (list, array): List of E(B-V)_SFD values (unit: mag).
            E4455 (list, array): E(440-550) values (unit: mag).
            Teff (list, array): Effective temperature (unit: mag).
            Logg (list, array): Logarithm of surface gravity.
            cal_lamb (bool): Whether to calculate effective wavelength.
            energy (bool): Whether to use energy for calculations, only matters when calculating effective wavelength, i.e. 'cal_lamb' is True.

        Returns:
            DataFrame: Calculated extinction values.
        """
        # Get E4455
        if self.model_func=='XP':
            if len(np.atleast_1d(E4455)) == 0:
                E4455 = self.Cal_E4455(sfdebv=sfdebv,ebv=ebv,Teff=Teff)
        else: 
            if len(np.atleast_1d(ebv)) != 0:
            # If other extinction curves are used, use ebv instead of E4455
            # Must enter ebv
                E4455 = ebv
            else:
                raise ValueError("Need to enter E(B-V) when using models other than 'XP'")
        # Input Check
        E4455,Teff,Logg = self.Check_and_Expand([E4455,Teff,Logg])
        amount = len(Teff)
        Band = np.atleast_1d(Band)
        # Obtain model extinction and effective wavelength
        A_model,λ_eff = self.model_extinction(Band,cal_lamb,energy)
        # Start star-by-star calculation
        Ext,LAMB = pd.DataFrame(columns=Band,dtype=np.float64),pd.DataFrame(columns=Band,dtype=np.float64)
        if len(Band)>1: 
            band_iter = tqdm(Band)
            band_iter.set_description("Progress of Bands")
        else: 
            band_iter = Band
        for B in band_iter:
            ExtLi,LAMBLi = np.empty(amount),np.empty(amount)
            # 检查是否为Band中的第一个元素
            if B == Band[0]:
                star_iter = tqdm(zip(Logg, Teff, E4455))
                star_iter.set_description("Processing Stars in the first band")
            else:
                star_iter = zip(Logg, Teff, E4455)
            for i,(g,t,e) in enumerate(star_iter):
                if g!=4.5: # Find the logg closest to the grid point
                    the_logg = self.logg_list[np.abs(self.logg_list - g).argmin()]
                else:
                    the_logg = g
                # When the amount is small, use Teff, ebv two-dimensional interpolation, which is slower
                if amount<1e6:
                    mini_A = A_model[B][:,the_logg]
                    f_interp2d = interp2d(self.E4455_list, mini_A.index, np.vstack(mini_A), kind='linear')
                    ExtLi[i] = f_interp2d(e,t)[0]
                    # Calculate effective wavelength
                    if cal_lamb:
                        mini_λ = λ_eff[B][:,the_logg]
                        f_interp2d = interp2d(self.E4455_list, mini_λ.index, np.vstack(mini_λ), kind='linear')
                        LAMBLi[i] = f_interp2d(e,t)[0]
                # When the amount is large, use ebv one-dimensional interpolation, which is faster
                elif amount>=1e6: 
                    # Find the Teff closest to the grid point
                    the_teff = self.teff_list[np.abs(self.teff_list - t).argmin()]
                    mini_A = np.array(A_model[B][the_teff,the_logg])
                    f_interp1d = interp1d(self.E4455_list, mini_A, kind='linear', fill_value="extrapolate")
                    ExtLi[i] = f_interp1d(e)
                    # Calculate effective wavelength
                    if cal_lamb:
                        mini_λ = np.array(λ_eff[B][the_teff,the_logg])
                        f_interp1d = interp1d(self.E4455_list, mini_λ, kind='linear',fill_value="extrapolate")
                        LAMBLi[i] = f_interp1d(e)
            Ext[B],LAMB[B] = ExtLi,LAMBLi
        if cal_lamb:
            return Ext,LAMB
        else: 
            return Ext


    #% Calculate the color excess (reddening) for given color, Teff, and E(440-550) 
    # 计算给定 Teff 和 E(440-550) 时准确的红化
    def star_reddening(self,Color,sfdebv=[],ebv=[],E4455=[],Teff=None,Logg=4.5):  
        """
        Calculates the precise reddening for given Teff, logg, and reddening in input colors using the XP extinction curve (default) or other extinction models. The input reddening can be E(440-550), E(B-V), or E(B-V)_SFD (derived from the SFD full-sky two-dimensional dust reddening map Schlegel et al. (1998)).

        Parameters:
            Color (str): Color index for which reddening is calculated.
            sfdebv (list, array): List of SFD E(B-V) values.
            ebv (list, array): List of E(B-V) values.
            E4455 (list, array): E(440-550) values.
            Teff (list, array): Effective temperature.
            Logg (list, array): Logarithm of surface gravity.

        Returns:
            DataFrame: Calculated color excess values.
        """
        # Splitting colors into bands
        Color = np.atleast_1d(Color)
        bandpairs = [c.split('-') for c in Color]
        if np.array(bandpairs).shape[1] != 2: 
            raise ValueError("Each element of the input color list should be split between two bands with a '-'")
        Band = np.unique(np.concatenate(bandpairs))
        
        # Get E4455, add bands, and obtain extinction
        if self.model_func=='XP':
            if len(np.atleast_1d(E4455)) == 0:
                E4455 = self.Cal_E4455(sfdebv=sfdebv,ebv=ebv,Teff=Teff)
            if '440' not in Band: Band = np.append(Band,['440'])
            if '550' not in Band: Band = np.append(Band,['550'])
            Ext = self.star_ext(Band, E4455=E4455, Teff=Teff, Logg=Logg)
            CE = pd.DataFrame(columns=Color)
            for (B1, B2), C in zip(bandpairs, Color):
                CE[C] = (Ext[B1] - Ext[B2]) / (Ext['440'] - Ext['550']) * E4455
        
        # If other extinction curves are used, use ebv instead of E4455
        else: 
            if len(np.atleast_1d(ebv)) == 0:
                # Must enter ebv
                raise ValueError("Need to enter E(B-V) when using models other than 'XP'")
            else:
                if 'Johnson.B' not in Band: Band = np.append(Band,['Johnson.B'])
                if 'Johnson.V' not in Band: Band = np.append(Band,['Johnson.V'])
                Ext = self.star_ext(Band, ebv=ebv, Teff=Teff, Logg=Logg)
                CE = pd.DataFrame(columns=Color)
                for (B1, B2), C in zip(bandpairs, Color):
                    CE[C] = (Ext[B1] - Ext[B2]) / (Ext['Johnson.B'] - Ext['Johnson.V']) * ebv
        
        return CE
