# EGFAnalysisTimeFreq

Surface Wave EGF/CF Phase Velocity Dipersion Analysis Software.

The software can do individual, semiautomated or automated phase velocity dispersion measurements.

by Huajian Yao, USTC

## Modified: 

Using the Automatic mode, the program can automatically pick the dispersion curves using a simple algorithm and save the dispersion images. The dispersion images can be input to DisperPicker to automatically pick the dispersion curves. Please read Yang et al. (2022) for details.

## Usage:

There is an example in the `Feidong_disp` folder. 

1. Run `EGFAnalysisTimeFreq.m`

2. The parameters I used is saved as `fd_config.png`

3. After selecting the parameters, click the `Start Processing` button to start saving the dispersion files.

For other details, please refer to `EGFAnalysisTimeFreq_Manual.pdf`

The dispersion image file is a 176x49 matrix and saved as a `dat`. The velocity range we set is 0.5-4 km/s and the step is 0.02 km/s. The period range is 0.2-5 s with a step of 0.1 s. 

## References:
Yang, S., Zhang, H., Gu, N., Gao, J., Xu, J., Jin, J., Li, J., and Yao H. (2022). Automatically Extracting Surface Wave Group and Phase Velocity Dispersion Curves from Dispersion Spectrograms Using a Convolutional Neural Network. Seismological Research Letters, doi: https://doi.org/10.1785/0220210280.
