# Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging

This repository contains the source code and data for the paper, "Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging," by Daiki Tamada et al.

## About The Project

This work introduces and evaluates a novel MRI method, Phase-Based Diffusion (PBD), for quantitative mapping of the apparent diffusion coefficient (ADC) and T2. Our approach uses an RF-phase modulated gradient echo (GRE) sequence to encode diffusion and T2 information directly into the signal phase.

This technique overcomes key limitations of conventional single-shot echo-planar imaging (SS-EPI), including geometric distortions and limited spatial resolution. The result is high-fidelity, high-resolution ADC and T2 maps from a single, motion-robust acquisition.

This repository will provide the MATLAB code used for:
*   The closed-form analytical signal model simulations.
*   Monte Carlo noise performance simulations.
*   Iterative reconstruction of ADC and T2 maps from acquired data.
*   Processing of the volunteer data presented in the paper.

---

## **Important Notice Regarding Code Availability**

**The core technique described in our manuscript is subject to a pending patent filed by our institution.**

Due to these intellectual property constraints, we are currently unable to release the source code publicly. We are working with our institution's patent management organization to obtain the necessary approvals for a public release.

We are fully committed to reproducibility and open science. **We will upload the complete source code and relevant data to this repository as soon as we receive clearance.**

In the meantime, for the purpose of peer review, we are able to provide the source code privately. If you are a reviewer and require access to the code, please contact the corresponding author.

We appreciate your understanding and patience.

---

## Contact

Daiki Tamada - dtamada@wisc.edu

Project Link: [https://github.com/dtamadauw/phase_based_diffusion_code](https://github.com/dtamadauw/phase_based_diffusion_code)

