# Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging

This repository contains the source code and data for the paper, "Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging," by Daiki Tamada et al.

## About The Project

This work introduces and evaluates a novel MRI method, Phase-Based Diffusion (PBD), for quantitative mapping of the apparent diffusion coefficient (ADC) and T2. Our approach uses an RF-phase modulated gradient echo (GRE) sequence to encode diffusion and T2 information directly into the signal phase.

This technique overcomes key limitations of conventional single-shot echo-planar imaging (SS-EPI), including geometric distortions and limited spatial resolution. The result is high-fidelity, high-resolution ADC and T2 maps from a single, motion-robust acquisition.

This repository will provide the MATLAB code used for:
* The closed-form analytical signal model simulations.
* Monte Carlo noise performance simulations.
* Iterative reconstruction of ADC and T2 maps from acquired data.
* Processing of the phantom and volunteer data presented in the paper.

## License

The primary source code in this repository, including the core PBD algorithms and associated scripts, is licensed under the **ACADEMIC SOFTWARE END USER LICENSE AGREEMENT**. Please refer to the `LICENSE` file in the root of this repository for the full terms and conditions.

## Third-Party Licenses

This repository incorporates certain third-party components, each under its own respective license. These licenses apply only to the specific portions of the code they cover.

* **Statistical Analysis Scripts:** The scripts `regression_line_ci.m` and `icc.m`, located in the `./tools/` directory, are subject to the license terms provided in `./tools/LICENSE_THIRDPARTY/LICENSE_ICC.txt`. This license applies only to these specific statistical code portions.
* **Colormaps:** The colormaps used in this project, located in `./tools/Colormaps/`, are licensed under a CC0 "no rights reserved" license. The full license details can be found in `./tools/LICENSE_THIRDPARTY/LICENSE_COLORMAP.txt`.

## Usage

All scripts in this repository are written in MATLAB.

To generate the figures and/or plots presented in the paper:
* Navigate to the respective figure directory (e.g., `Figure1-3`, `Figure5`, ..., `Figure8`).
* Each directory contains a main script (e.g., `Figure1.m`, `Figure2.m`).
* Simply run the `Figure*.m` script within its directory to generate the corresponding figures or plots. All dependent scripts required for execution are located in the `./tools/` directory.

## Contact

Daiki Tamada - dtamada@wisc.edu

Project Link: [https://github.com/dtamadauw/phase_based_diffusion_code](https://github.com/dtamadauw/phase_based_diffusion_code)






