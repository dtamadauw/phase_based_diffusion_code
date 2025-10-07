
#Author: Daiki Tamada
#Affiliation: Department of Radiology, University of Wisconsin-Madison
#Date: 10/1/2025
#Email: dtamada@wisc.edu
#By downloading, installing, or otherwise accessing or using the Software , you (“Recipient”) agree to receive and use the above-
#identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
#obligations and restrictions you are not permitted to download, install,
#execute, access, or use the SOFTWARE:

#1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE (“Provider’s Rights”) to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation (“WARF”). 
#2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider’s Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient’s responsibility. 
#3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider’s proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
#4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
#5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN “AS IS” BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
#6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER’S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient’s use of the SOFTWARE.
#7.	No use of names/trademarks.  Recipient shall not use Provider’s name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
#8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient’s possession, custody, or control.	


# --- 1. SETUP: Import necessary packages ---
using LsqFit
using Plots
using Random
using Statistics
using MAT


# Set a seed for reproducibility of the random noise
Random.seed!(1234)

# --- 2. MODEL AND PARAMETER DEFINITION ---

# Define the T2 decay model function
# t: A vector of echo times (TEs)
# p: A vector of parameters [S0, T2]
# S(t) = S0 * exp(-t / T2)
@. t2_decay_model(t, p) = p[1] * exp(-t / p[2])

# Define Ground Truth Parameters (the "true" values we want to recover)
true_S0 = 1.0  # Initial signal intensity
true_T2 = 80.0   # True T2 relaxation time (e.g., in ms)
T2_true_s_values = [100, 120, 200];

# Define Acquisition Parameters
TEs = [15, 30, 60]; # Vector of echo times
TR = 500;
SNR = 2:0.5:50;
noise = 1 ./ SNR;
# Define Simulation Parameters
num_simulations = 100000  # Number of Monte Carlo iterations

println("--- Simulation Configuration ---")
println("True S0: $true_S0")
println("True T2: $true_T2 ms")
println("Number of Simulations: $num_simulations")
println("--------------------------------\n")


# --- 3. FITTING FUNCTION ---

"""
    fit_t2_signal(model, tes, all_signals, p0)

Performs a non-linear least squares fit on a collection of signals.

# Arguments
- `model`: The function defining the mathematical model (e.g., t2_decay_model).
- `tes`: A vector of echo times.
- `all_signals`: An array of signal data arrays to be fitted.
- `p0`: A vector of initial parameter guesses.

# Returns
- A vector of estimated T2 values. Failed fits are marked as `NaN`.
"""
function fit_t2_signal(model, tes, all_signals, p0, lower_bounds, upper_bounds)
    num_fits = size(all_signals,1)
    estimated_t2s = zeros(Float64, num_fits, size(all_signals,3))
    
    println("Fitting each dataset...")
    for noise_ind in 1:size(all_signals,3)
        for i in 1:size(all_signals,1)
            signal = all_signals[i,:,noise_ind]
            try
                fit_result = curve_fit(model, tes, signal, p0, lower=lower_bounds, upper=upper_bounds)

                estimated_t2 = fit_result.param[2]
                
                if isnan(estimated_t2)
                    println("Warning: Fit resulted in NaN at iteration $i. Skipping.")
                end
                estimated_t2s[i,noise_ind] = estimated_t2
            catch e
                println("Catched Warning: Fit failed at iteration $i. Skipping.")
                estimated_t2s[i,noise_ind] = NaN
            end
        end
    end
    return estimated_t2s
end


# --- 3a. Prepare noisy data ---

noisy_signals = zeros(num_simulations,length(TEs),length(noise),length(T2_true_s_values));
Random.seed!(321);

for t2_idx = 1:length(T2_true_s_values)
    T2_true_s = T2_true_s_values[t2_idx];
    S_MESE = 1.0 * exp(-TEs[1]/T2_true_s) * (1.0-exp(-TR/1000));

    for ii=1:length(noise)
        ideal_signal = t2_decay_model(TEs', [true_S0, T2_true_s])
        noise_array = S_MESE .* noise[ii] .* complex.(randn(num_simulations,length(TEs)),randn(num_simulations,length(TEs)))/sqrt(2);
        noisy_signals[:,:,ii,t2_idx] = abs.(ideal_signal .+ noise_array);
    end
end


# --- 4. MONTE CARLO SIMULATION ---

# --- Step B: Fit all noisy datasets ---
p0 = [1.0, 75.0] # A good initial guess, slightly off from the truth.
lower_bounds = [-1.0, -1000.0]       # S0 >= 0, T2 >= 0
upper_bounds = [1.0, 1000.0]    # S0 has no upper limit, T2 <= 1000

estimated_t2s = zeros(num_simulations,length(noise),length(T2_true_s_values));
for t2_idx = 1:length(T2_true_s_values)
     @time estimated_t2s[:,:,t2_idx] = fit_t2_signal(t2_decay_model, TEs, noisy_signals[:,:,:,t2_idx], p0, lower_bounds, upper_bounds);
end

data_to_save = Dict(
    "MESET2_jul" => estimated_t2s
)

matwrite("MESET2_Julia.mat", data_to_save)
