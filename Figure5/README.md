# Figure 5 Simulation Code

This folder contains the MATLAB code used to generate the Monte Carlo simulation results presented in **Figure 5** of the main paper. The simulation compares the quantitative accuracy and noise performance of the proposed Phase-Based Diffusion (PBD) method against conventional reference methods (SS-EPI for ADC, SESE for T2) under time-normalized conditions.

---

## Files

* `Figure5.m`: The main script that orchestrates the entire simulation, from parameter setup to final plotting.
* `pbd_recon.m`: A function that estimates ADC and T2 values from noisy PBD signals using the analytical model.
* `ssepi_recon.m`: A function that estimates ADC from noisy Single-Shot Echo-Planar Imaging (SS-EPI) signals.
* `mese_recon.m`: A function that estimates T2 from noisy Multi-Echo Spin-Echo (MESE) or Single-Echo Spin-Echo (SESE) signals.

---

## How to Run

To execute the simulation and reproduce the plots shown in Figure 5, simply run the main script in MATLAB:

```matlab
Figure5a
Figure5b
```

## Detailed Implementation

* **Key Steps:**

    1.  Parameter Definition: The script begins by defining all necessary parameters:

        * Ground Truth Values: Sets the true ADC (800, 1300, 1800 µm²/s) and T2 (100, 120, 200 ms) values for the simulation.

        * Sequence Parameters: Defines the acquisition parameters (TR, TE, flip angle, bandwidth, etc.) for PBD, SS-EPI, and SESE, mirroring those used in the in vivo experiments.

        * Simulation Controls: Sets the number of Monte Carlo iterations and the range of Signal-to-Noise Ratios (SNRs) to be tested.

    2.  Time Normalization: To ensure a fair comparison, the script calculates an effective number of signal averages (NSA) for each method that would be achievable within a fixed total acquisition time (10 minutes). This normalization accounts for differences in the acquisition time per slice for each sequence.

    3.  Signal Generation: For each set of ground truth parameters, the script generates the corresponding noiseless complex signals for each imaging method (PBD, SS-EPI, SESE) using their respective physical models.

    4.  Monte Carlo Loop: The script iterates through each SNR level and each set of ground truth values. Inside the loop:

        * Noise Addition: Complex Gaussian noise is added to the noiseless signals. The noise variance is carefully scaled based on the receiver bandwidth and the time-normalized NSA for each method.

        * Parameter Estimation: The noisy signals are passed to the respective reconstruction functions (pbd_recon.m, ssepi_recon.m, mese_recon.m) to estimate the ADC and T2 values.

        * Store Results: The estimated values from all iterations are stored for statistical analysis.

    6.  Statistical Analysis & Plotting: After the loops complete, the script calculates the mean, median, and standard deviation of the estimated ADC and T2 distributions for each condition. Finally, it generates the plots presented in Figure 5, showing the relative bias and standard deviation as a function of SNR.

