# Figure S5 Simulation Code

This folder contains the MATLAB code to reproduce the simulation results presented in **Figure S5** of the supporting information. This simulation was designed to validate the choice of regularization weights (λ for Total Variation and β for L2) used for the *in vivo* reconstructions in the paper.

---

## Files

* `FigureS5.m`: The main script to run the simulation and generate the plots shown in Figure S5.
* `TV_mapping_penalty.m`: A function that performs the iterative reconstruction of ADC and T2 maps from noisy Phase-Based Diffusion (PBD) signals, applying TV and L2 regularization.

---

## How to Run

To run the simulation and generate the figure, simply execute the `FigureS5.m` script in MATLAB.


## Detailed Implementation

`FigureS5.m`

This script performs a numerical simulation to evaluate the trade-off between bias and variance in ADC and T2 estimates when using different regularization weights.


* **Key Steps:**

    1.  Phantom Creation: A 2D numerical phantom is created with three distinct tissue compartments, each with different ground truth ADC and T2 values.

    2.  Signal Generation: For each pixel in the phantom, noiseless PBD signals are generated using the analytical model described in the main paper.

    3.  Noise Addition: Gaussian noise is added to the complex signals to achieve a signal-to-noise ratio (SNR) of 5.

    4.  Iterative Reconstruction: The script then iterates through a range of λ (TV) and β (L2) regularization weights. For each pair of weights, the TV_mapping_penalty.m function is called to reconstruct the ADC and T2 maps from the noisy signals.

    5.  Bias and Variance Calculation: The bias and standard deviation of the ADC and T2 estimates are calculated for each tissue compartment and for each set of regularization weights.

    6.  Plotting: Finally, the script generates the plots seen in Figure S5, showing the ADC and T2 bias and standard deviation as a a function of the regularization weights.
 
`TV_mapping_penalty.m`

This function implements the iterative algorithm for TV-regularized T2 and ADC mapping, as described in Algorithm 1 in the supporting information. It takes noisy phase images and regularization parameters as input and iteratively solves the optimization problem to estimate the T2 and ADC maps.

* **Inputs:**

    * y_H, y_L: The phase images for high and low gradient moments.

    * lambda1, lambda2: The TV regularization parameter (λ).

    * beta1, beta2: The L2 regularization parameter (β).

    * LUTs: Lookup table and a spline-fitted function representing the analytical signal model for rapid computation.

* **Output:**

    * x1_est: The estimated T2 map.

    * x2_est: The estimated ADC map.

* **Implementation Steps:**

    1.  Initialize Maps:

        * Initialize the T2 map (x1) and ADC map (x2) with zeros.

    2.  Begin Iterative Loop (for k = 1 to max_iterations):

        * Compute Predicted Phases & Gradients: Use the spline function to compute the predicted phases (f1, f2) based on the current x1 and x2 maps.

        * Calculate Residuals: Compute the difference between the measured phases (y1, y2) and the predicted phases (f1, f2) to get the residuals (r1, r2).

        * Compute Gradients:

          * Calculate the gradient of the data fidelity term (grad_data) from the residuals.

          * Calculate the gradient of the TV penalty term (grad_TV) from the current x1 and x2 maps.

          * Calculate the gradient of the L2 penalty term (grad_L2) from the current x1 and x2 maps.

        * Calculate Total Gradient: Combine the individual gradients, weighted by lambda and beta, to get the total gradient for both x1 and x2.

        * Update Step Size: Use the Barzilai-Borwein method to determine an adaptive step size for the update.

        * Update Maps: Update the T2 and ADC maps by taking a step in the direction of the negative gradient.
          
        * Check for Convergence: Evaluate the relative change in the maps over the last 5 iterations. If the change is less than the tolerance (tol), exit the loop.
      
    3.  End Loop: The loop terminates when convergence is reached or maxIter is exceeded.




