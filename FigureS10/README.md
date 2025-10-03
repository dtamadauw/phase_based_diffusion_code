
# PBD Reconstruction with B1+ Correction

This folder contains the MATLAB code to perform the joint ADC/T2 mapping from Phase-Based Diffusion (PBD) data, including a voxelwise correction for B1+ field inhomogeneity.

This implementation was developed to validate the B1+ sensitivity of the PBD method, as described in the Supporting Information of the manuscript "Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging."

## Description

The core of the B1+ correction is the extension of the forward signal model to be a function of three variables: **T2**, **ADC**, and the **local flip angle (FA)**. This is achieved by pre-calculating a 3D look-up table (LUT) and using a tri-variate spline interpolant during the iterative reconstruction.

The workflow consists of two main scripts:

1.  `build_LUTs_B1.m`: Generates and saves the 3D look-up table.
2.  `TV_mapping_fast_B1.m`: Performs the iterative, regularized reconstruction using the 3D LUT and an acquired B1+ map.

## Workflow

### 1. Build the 3D Look-Up Table (LUT)

**Script:** `build_LUTs_B1.m`

This script pre-calculates the expected PBD signal phases for a grid of tissue properties (T2, ADC) and sequence parameters (local flip angle).

* **Inputs:**
    * `params`, `ADC_v`, `FA_v`: Sequence parameters (TR, RF phase increment, Gradient Moments, etc.) defined within the script.
* **Process:**
    1.  It iterates through every combination of T2, ADC, and local FA.
    2.  For each combination, it calls the `analytical_SPGR_W_diffusion.m` function to calculate the complex steady-state signals for the low (L) and high (H) gradient moments.
    3.  It calculates the two signal phases (`phase_L` and `phase_H`) that are the outputs of the forward model.
    4.  The results (models) are stored in two 3D matrices, `mdl_theta3D` and `mdl_theta23D`.
* **Output:**
    * `LUTs`: 3D LUTs (ADC*T2*beta).

### 2. Perform B1+-Corrected Reconstruction

**Script:** `TV_mapping_fast_B1.m`

This script loads the acquired PBD data, the B1+ map, and the 3D LUT, then performs the iterative reconstruction to estimate the ADC and T2 maps.

* **Inputs:**
    * `sig`: The four complex-valued images from the PBD acquisition.
    * `beta`: A co-registered B1+ map, where each voxel value represents the local flip angle scaling factor (e.g., 1.0 = nominal FA, 0.9 = 90% of nominal FA).
    * `LUTs`: The pre-calculated 3D look-up table.
    * Regularization parameters (`lambda`, `beta`) and optimization settings defined within the script.
* **Process:**
    1.  **Data Preparation:** Loads the data and calculates the background-corrected input phase maps.
    2.  **Iterative Optimization:** The script iterates through each voxel of the image.
        * For a given voxel, it retrieves the local flip angle from the co-registered `beta`.
        * It then solves the non-linear, regularized inverse problem (Equation 14 from the main text) to find the ADC and T2 values that best match the measured phases.
        * Crucially, during the data fidelity calculation at each iteration, the spline interpolant is evaluated using the **voxel-specific flip angle**. This ensures that the forward model is tailored to the local B1+ environment.
* **Output:**
    * `x1_est`, `x2_est`: The final, B1+-corrected quantitative ADC and T2 maps.

## How to Run

1.  Run `FigureS10.m` to generate B1+-corrected ADC and T2 maps.

