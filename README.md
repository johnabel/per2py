# Code for processing PER2iLUC
This repository contains code for processing SCN color-switching Per2iLuc data in a high-throughput manner, with a manuscript under review. This package uses the scientific Python stack to 

Below, please find instructions for the use of this package, and a step-by-step walkthrough of the analysis it performs.

# Table of Contents


# Installation
We recommend cloning this repository, or downloading a .zip copy of it. The methods for processing default to analysis of the data within the `Demo` directory, and enable the user to direct the analysis at other directories. Several dependencies are necessary

## Dependencies
The scripts in this package are written in Python 2.7. The following dependencies are required for running this code. We recommend using Anaconda to manage the Python environment in which this analysis is done, and using Anaconda's pip install to install the necessary packages. To check for pip, run in a terminal command line interpreter:
```
python -m ensurepip
```

If pip is installed, the associated packages may be installed with:
```
pip install -U [packagename]
```

The packages necessary for running this core are

# Using per2iluc automated data analysis
We have provided two interfaces for running the data analysis tools within this package. A Jupyter Notebook method is provided for a simplified interface.

## Running analysis in a Jupyter/iPython Notebook
The Jupyter Notebook provides a simple interface to the computational tools within this package.

## Running analysis using a terminal and command line
In addition to a Jupyter Notebook, we have provided a command line interface to the analysis tools.

The data may be imported directly from the imageJ output .xls file. If so, set:
```
pull_from_imagej = True
input_folder = [path to folder containing data]
input_ij_file   = [name of imageJ file with no extension]
data_type = [type of data being imported]
```
Note that data_type will be used to name all future outputs of the file. For example, the raw signal and raw XY tracking will be named:
```
raw_signal = output_folder +data_type+'_signal'+output_extension
raw_xy = output_folder +data_type+'_XY'+output_extension
```
If the data is _not_ being taken from an imageJ file but has already been converted into `_signal.csv` and
`_XY.csv` files, one need only adjust `data_type` so that `raw_signal` and `raw_xy` look in the correct
locations. The flag for pulling from imageJ should be set to False.


## Interpreting results
The data produced during this analysis includes detrended signal, detrended and denoised signal, z-scored detrended and denoised signal, Lomb-Scargle periodogram results, a sinusoidal approximation to the data, the phases of the sinusoidal approximation, and parameters describing the oscillation for each trajectory provided. These results are exported as `.csv` files in the following manner.
    
|    Suffix | Data        |
|-----------|-------------|
|`_signal.csv`| Raw bioluminescence.|
|`_XY.csv`| Raw location values.|
| `_signal_detrend.csv`| Data interpolated and detrended via HP filter.|
| `_signal_detrend_denoise.csv`| Detrended, denoised, 12h truncated data. | 
|`_XY_detrend_denoise.csv` | Detrended, denoised, 12h truncated locations. | 
|`_zscore.csv` | z-scored detrended and denoised data |
|`_lombscargle.csv`| Normalized Lomb-Scargle periodogram values. |
|`_cosine.csv`| The sinusoid fit to each cell. Only created if cell is rhythmic.|
|`_cosine_phases.csv`| The sinusoid phases corresponding to the fit.|
|`_oscillatory_params.csv`| Rhythmic, estimated phase, normalized LS circadian peak, period (h), amplitude, decay, $R^2$.|
<br><br>

# Step-by-step details on bioluminescence processing
1. **Import the data.**
    Details on the data import are provided in the instructions section.
<br><br>
    
2. **Detect outliers and interpolate missing intermediate values in the dataset.**
    Outliers > 5 standard deviations from the mean are removed. Note that interpolation is linear and only applied to intermediate points (i.e., no extrapolation is performed).<br><br>

3. **Detrending via Hodrick-Prescott filter.**
    A Hodrick-Prescott filter is applied to detrend the data. We (JH Abel and Y Shan) found that this performs better than polynomial detrending in that it does not overfit the trend and begin to fit the oscillatory components. Parameter selection for the filter is performed as in Ravn, Uhlig 2004, and implemented as in St. John and Doyle 2015.<br><br>
    
4. **Eigendecomposition and reconstruction of the signal to denoise the data.**
    An eigendecomposition of the autocorrelation matrix is performed. Any eigenvectors with eigenvalue >5% of the total sum of all eigenvalues is kept, and the signal is reconstructed from these eigenvectors. The process is theoretically explained here: http://www.fil.ion.ucl.ac.uk/~wpenny/course/subspace.pdf<br>
     **Then truncate first 12 h to remove artifact.**
    All trajectories truncated to start at 12 h. If a cell is first found after 12 h, no time is truncated since it does not exhibit an artifact from plating.<br><br>
    
5. **Lomb-scargle periodogram to determine which cells have oscillation.**
    The Lomb-Scargle periodogram is used to identify statistically significant circadian oscillation. Here,
    we have defined this as having a dominant peak of $P<0.05$ between periods of 18 and 30 h. The method for 
    periodogram calculation and corresponding P-values is from WH Press, Numerical Recipes, 1994, p576. The final periodogram is normalized using the standard normalization in Astropy or Scipy. We applied this to the detrended, denoised, and truncated data so that the noise difference between channels is corrected, and no statistical difference in peak distribution exists within a WT AVPCre SCN between channels (t-test, $P=0.75$). See `lomb_testing.py` for details.<br><br>
    
6. **Fit a damped sinusoid to the data to yield sine fit, phase, amplitude, damping rate, goodness of fit.**
    All rhythmic cells are now fit by a damped sinusoid plus a polynomial to fit any non-oscillatory trend in the data. Only the sinusoid data is returned.<br><br>
    
7. **Export data from the analyses.** 
    The data produced at most steps of this process is then saved in the output folder, as delineated in the . 
    
|    Suffix | Data        |
|-----------|-------------|
|`_signal.csv`| Raw bioluminescence.|
|`_XY.csv`| Raw location values.|
| `_signal_detrend.csv`| Data interpolated and detrended via HP filter.|
| `_signal_detrend_denoise.csv`| Detrended, denoised, 12h truncated data. | 
|`_XY_detrend_denoise.csv` | Detrended, denoised, 12h truncated locations. | 
|`_zscore.csv` | z-scored detrended and denoised data |
|`_lombscargle.csv`| Normalized Lomb-Scargle periodogram values. |
|`_cosine.csv`| The sinusoid fit to each cell. Only created if cell is rhythmic.|
|`_cosine_phases.csv`| The sinusoid phases corresponding to the fit.|
|`_oscillatory_params.csv`| Rhythmic, estimated phase, normalized LS circadian peak, period (h), amplitude, decay, $R^2$.|
<br><br>

**Lastly, the code generates summary plots at random.**
    Summaries of the data processing are performed for three cells selected at random. Subpanels are A-F, left to right, top to bottom. (A) Raw bioluminescence data. (B) Raw bioluminescence with trend (gold). (C) Detrended bioluminescence. (D) Eigendecomposition (gold) with threshold for reconstruction (red). (E) Detrended and denoised reconstructed signal. (F) Lomb-Scargle periodogram with rhythmic test in y-axis. (G) Sinusoid fit to data with $R^2$ value in y-axis label. (H) Detrended, denoised, and sine fit plotted simultaneously.

# License

This package is licensed under the GNU General Public Licesnce vserion 3. Please see LICENSE for more information.

# Authors
This package was written by John H. Abel and Yongli Shan.

# Funding
This package was funded by NIH/NIA F32 AG064886 and NIH/NHLBI T32 HL007901 (JHA).