# XP_Extinction_Toolkit

[English](README.md) | [中文](README-zh.md)

`XP_Extinction_Toolkit` 是一个为天文研究者设计的 Python 库。它基于 Gaia DR3 和 LAMOST DR7 中大约37万个高质量样本计算得到的中值消光曲线，提供了一组适用于光谱和测光数据的消光改正工具。这些工具可以在已知消光的情况下提供对天体光谱和任意波段测光数据的消光改正。更多信息见我们的文章 [(Zhang et al., 2024)](https://iopscience.iop.org/article/10.3847/1538-4357/ad613e)。
包含以下功能(后附使用说明)：

1. ### **函数 `ext_curve`**
   获取银河系的高精度中值消光曲线。

2. ### **函数 `deredden`**
   对输入光谱进行消光改正，适用范围是 2900 - 53414 Å。

3. ### **函数 `Cal_E4455`**
   将 E(B-V) 或 E(B-V)_SFD 通过经验公式转换成 E(440-550)。

4. ### **函数 `star_ext`**
   可以计算给定 Teff、logg 和红化时，输入**波段**的精确**消光**和**有效波长**。输入红化可以是E(440-550)、E(B-V) 或者 E(B-V)_SFD（指取自SFD全天二维尘埃红化图 [Schlegel et al. (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract)）。本函数可使用 XP 消光曲线（默认）或`dust_extinction` 包提供的各类 Rv 依赖的消光曲线模型。

5. ### **函数 `star_reddening`**
   可以计算给定 Teff、logg 和红化时，输入**颜色**的精确**红化**，使用 XP 消光曲线（默认）或其他消光模型。输入红化可以是E(440-550)、E(B-V) 或者 E(B-V)_SFD。

6. ### **函数 `model_extinction`**
   计算给定输入波段的模拟消光。使用 XP 消光曲线（默认）或其他消光模型。返回一个不同 Teff、E(B-V)、logg 时的各波段模拟消光值的 DataFrame。

其中函数4、5、6对于测光波段（颜色）的消光（红化）估计是基于BOSZ光谱库，滤光片通带和实测消光曲线的，属于理论与实测结合的计算，可能会和纯实测结果有所不同。当需要计算 GALEX，PS1，SDSS，Gaia，2MASS 和 WISE 波段的消光（或红化）时，我们推荐使用另一个基于实测的消光系数包[`extinction_coefficient`](https://github.com/vnohhf/extinction_coefficient) 来进行消光改正。该包在 0 - 0.5 mag 的E(B-V)范围和 4000 - 10000 K 的温度范围内大多有效。


# 如何安装
### 使用 pip
~~~
# PyPI (推荐)
pip install extinction_coefficient

# 从本 github 库
pip install git+https://github.com/vnohhf/extinction_coeffcient.git
~~~

### 从源代码安装
从 git 仓库中下载 extinction_coefficient 后，可以从源代码中安装
(https://github.com/vnohhf/extinction_coeffcient/):
~~~
python setup.py install
~~~

# 使用指南
## 1. 获取 XP 消光曲线并绘图
~~~python
import numpy as np
import matplotlib.pyplot as plt
from xp_extinction_toolkit import ExtinctionToolkit

# 初始化 ExtinctionToolkit
ExtTool = ExtinctionToolkit()
# 获取 XP 消光曲线
λ, k_x_55 = ExtTool.ext_curve()

# 绘图
fig,ax = plt.subplots(2,1, figsize=(5,6), tight_layout=True, dpi=300)
ax[0].plot(λ[λ<10200], k_x_55[λ<10200], lw=1)
ax[1].plot(1e4/λ, k_x_55, lw=1)
ax[0].set(title='XP extinction curve (optical part)', ylabel='E(λ-55)/E(44-55)', xlabel='Wavelength [Å]')
ax[1].set(title='XP extinction curve (optical and infrared)', ylabel='E(λ-55)/E(44-55)', xlabel='Wavelength number [1/μm]')
for axi in ax:
    axi.minorticks_on()
    axi.grid(True, ls=':', color='dimgray', lw=0.2, which='major', zorder=1)
~~~

## 2. 对输入光谱进行消光改正
### 2.1 单一光谱
~~~python
ExtTool = ExtinctionToolkit()
# 获取 BOSZ 光谱
BOSZ = np.load(ExtTool.start_path+'/BOSZ_spectra.npy', allow_pickle=1) 

# 获取波长，已修改为匹配 XP 消光曲线。
λ, _ = ExtTool.ext_curve()

# 随机选择一个作为示例，获取原始光谱及其 Teff
ind = np.random.choice(range(len(BOSZ)))
original_spec = BOSZ[ind]
teff = ExtTool.teff_list[ind] 

# 以 E(B-V) 为 0.1 对原始光谱进行红化。
# 输入的 ebv 转换为 E4455 需要温度信息
reddened_spec = original_spec / ExtTool.deredden(ebv=0.1, Teff=teff, wave=λ)

# 对红化的光谱进行去红化
# （这里 intrinsic_spec 即等于 original_spec，仅展示去红化的用法）
intrinsic_spec = reddened_spec * ExtTool.deredden(ebv=0.1, Teff=teff, wave=λ)

# 绘图
fig, ax = plt.subplots(1, 1, figsize=(8, 4), tight_layout=True, dpi=300)
ax.plot(λ, intrinsic_spec, lw=0.9, label='本征光谱')
ax.plot(λ, reddened_spec, lw=0.9, label='红化光谱')
ax.set(title=f'Teff = {teff} [K]', ylabel='Flux [$cm^{-2} s^{-1} Å^{-1}$]', xlabel='波长 [Å]')
ax.minorticks_on()
ax.legend()
~~~

### 2.2 多光谱的红化/去红化
~~~python
ebvlist = np.random.uniform(low=0, high=2, size=len(BOSZ))
reddened_spec_matrix = 1 / ExtTool.deredden(ebv=ebvlist, Teff=ExtTool.teff_list, wave=λ) * BOSZ
~~~

## 3. 将 E(B-V) 或 E(B-V)_SFD 转换为 E(440-550)
### 3.1 输入 E(B-V)
~~~python
from xp_extinction_toolkit import ExtinctionToolkit
ebv = [0.1, 0.3, 0.5]
Teff = [4000, 6000, 8000]
E4455 = ExtTool.Cal_E4455(ebv=ebv, Teff=Teff)
~~~

### 3.2 输入 E(B-V)_SFD
~~~python
sfdebv = 0.2
E4455 = ExtTool.Cal_E4455(sfdebv=sfdebv)
~~~

## 4. 计算 1000 颗恒星的四个波段的消光
~~~python
import numpy as np
from xp_extinction_toolkit import ExtinctionToolkit

ExtTool = ExtinctionToolkit()

# 检查内置滤光片通带
print(ExtTool.Filters.keys())

# 生成 Teff 和 E(440-550) 的随机值
Band = ['2MASS.Ks', 'PS1.i', 'SDSS.g', 'WISE.W2']
Teff = np.random.uniform(low=4000, high=10000, size=1000)
E4455 = np.random.uniform(low=0, high=2, size=1000)

# 计算定义波段的消光
Extinction = ExtTool.star_ext(Band, E4455=E4455, Teff=Teff)
~~~

## 5. 根据 Teff 和 E(B-V)，计算恒星的红化
### 5.1 使用内置滤光片
本程序包内置了 GCPD_Johnson、GALEX、PS1、SDSS、Gaia3、2MASS 和 WISE 的通带。
~~~python
# 定义 E(B-V)_SFD 和 Teff（这里也可以使用array、tuple等）
sfdebv = [0.3, 0.7, 1, 1.3]
Teff = [6542, 6000, 3900, 7659]

# 定义色指数
Color = ['Johnson.B-Johnson.V', 'GAIA3.G-GAIA3.Grp']

# 计算红化
Reddening = ExtTool.star_reddening(Color, sfdebv=sfdebv, Teff=Teff)
~~~

### 5.2 使用外部输入滤光片
请将 '/path/to/new/filters/' 替换为实际的滤光片存放的绝对路径。
~~~python
ExtTool = ExtinctionToolkit(filter_dir='/path/to/new/filters/')
Band = ['New_Band']
Extinction = ExtTool.star_ext(Band, E4455=0.1, Teff=5000)
~~~

## 6. 计算不同 SED 时的模拟消光
### 6.1 使用 XP 消光曲线
~~~python
ExtTool_XP = ExtinctionToolkit(model='XP')
Band = ['GALEX.NUV', 'PS1.z', 'SDSS.u', 'WISE.W4']

# 计算模拟消光
Model_Ext, _ = ExtTool_XP.model_extinction(Band, cal_lamb=False, energy=False)
~~~

### 6.2 获取有效波长
设置 'cal_lamb' 为 True 以计算有效波长（单位：埃（Å））。

~~~python
# 计算模拟消光和相应的有效波长
Model_Ext, Model_λeff = ExtTool_XP.model_extinction(['GALEX.NUV', 'PS1.z', 'SDSS.u'], cal_lamb=True, energy=False)
~~~

部分项目滤光片类型为光子计数（如 Gaia、GALEX、PAN-STARRS、SDSS），部分为能量计数（如 2MASS、WISE）。使用能量计数类型的滤光片时，将 'energy' 设置为 True 以获取正确的有效波长。

~~~python
Model_Ext, Model_λeff = ExtTool_XP.model_extinction(['2MASS.H', 'WISE.W3'], cal_lamb=True, energy=True)
~~~

### 6.3 使用 [`dust_extinction`](https://dust-extinction.readthedocs.io/en/stable/index.html) 提供的 Rv 依赖消光曲线
~~~python
from dust_extinction.parameter_averages import CCM89, F99, F19, G23

# 使用 F19 消光曲线
ExtTool_F19 = ExtinctionToolkit(model=F19)
Model_Ext, Model_λeff = ExtTool_F19.model_extinction(['PS1.z', 'SDSS.u'], cal_lamb=True, energy=False)

# 使用不同 Rv 值（默认值为 3.1）的 G23 模型
ExtTool_G23 = ExtinctionToolkit(model=G23, Rv=2.5)
Model_Ext, Model_λeff = ExtTool_G23.model_extinction(['2MASS.H', 'WISE.W4'], cal_lamb=True, energy=True)
~~~
