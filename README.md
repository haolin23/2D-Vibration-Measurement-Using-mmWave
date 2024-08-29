<div align="center">

# 2D-Vibration-Measurement-Using-a-Single-mmWave-Radar

</div>

## Overview
- m2Vib is a novel method for measuring 2D vibrations using millimeter waves. This method leverages **a single mmWave radar** and a reflective surface to obtain the 2D vibration observations. The core idea is to decompose the 2D vibration into a set of 1D vibrations in orthogonal directions and leverage the **none-line-of-sight (NLoS) signal** positively. 

## Key Features
- Both mmWave line-of-sight (LoS) and NLoS signals carry the target's vibration signal.
- **A geometric model** is established to unveil the relationship between 2D vibration displacement and mmWave LoS and NLoS signals.
- A new **retroreflector** is designed to enhance LoS and NLoS signals.
- **Advanced signal processing methods** are used to obtain fine-grained geometric parameters.
- m2Vib can achieve sub-mm level 2D vibration measurement.

## Methodology
- **Target Detection and Vibration Signal Extraction**: The mmWave signal is processed by **range-FFT** and **beamforming** to discover the target, then we extract the phase and remove the static component to get the vibration signal.
- **Distance and Angle Estimation**: **Chirp-Z transform** and **Multiple Signal Classification (MUSIC)** methods are used to estimate the distance and angle parameters with high resolutions.
- **2D Vibration Estimation**: Based on the established geometric model, Taylor's expansion is applied to calculate the 2D vibration **along the X-axis and Y-axis**. The displacement and frequency along each axis and the trajectory can then be obtained.

##
- **More details can be found in our paper "m2Vib: 2D Vibration Measurement Using a Single mmWave Radar"**.
- You can get partial mmWave data from [here](https://pan.baidu.com/s/1NusPUI72karg5FebRWfwug?pwd=5zxc).

