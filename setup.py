from setuptools import setup, find_packages

setup(
    name='XP_Extinction_Toolkit',
    version='1.0',
    packages=find_packages(),
    description='Extinction correction package using extinction curves derived from XP spectra',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Ruoyi Zhang',
    author_email='zry@mail.bnu.edu.cn',
    url='https://github.com/vnohhf/XP_Extinction_Toolkit',
    include_package_data=True,
    package_data={
        'xp_extinction_toolkit': ['filters/*', '*.npy', '*.csv'],
    },
    install_requires=[
        'numpy',
        'pandas',
        'tqdm',
        'astropy',
        'scipy',
    ],
    entry_points={
        'console_scripts': [
            'extinction_toolkit=xp_extinction_toolkit.ExtinctionToolkit:main_function',  # 根据实际情况调整
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
)
