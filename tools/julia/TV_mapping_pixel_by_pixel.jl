
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



using LinearAlgebra, Statistics

# Include the original vectorized function file to use its helper functions
include("TV_mapping_no_penalty_ori.jl")

"""
    estimate_single_pixel(y1_val, y2_val, LUTs)

Performs the gradient descent optimization for a single pixel.
This is a scalar version of the logic inside `TV_mapping_no_penalty_ori`.
"""
function estimate_single_pixel(y1_val::Float64, y2_val::Float64, LUTs)
    # Get spline models from LUTs
    mdl_theta = LUTs.mdl_theta
    mdl_theta2 = LUTs.mdl_theta2

    # Initial guess for this single pixel (starting at 0)
    x1, x2 = 0.0, 0.0
    x1_old, x2_old = 0.0, 0.0
    grad_x1_old, grad_x2_old = 0.0, 0.0
    
    max_iter = 250
    tol = 1e-8
    
    # Optimization loop for the single pixel
    for iter in 1:max_iter
        # Calculate residuals (all scalar operations)
        f1 = mdl_theta.sp(x1, x2)
        f2 = mdl_theta2.sp(x1, x2)
        r1 = wrapToPi(f1 - y1_val)
        r2 = wrapToPi(f2 - y2_val)

        # Calculate gradients for the pixel
        grad_x1_f1 = mdl_theta.dx1(x1, x2)
        grad_x2_f1 = mdl_theta.dx2(x1, x2)
        grad_x1_f2 = mdl_theta2.dx1(x1, x2)
        grad_x2_f2 = mdl_theta2.dx2(x1, x2)

        grad_x1 = 2 * (grad_x1_f1 * r1 + grad_x1_f2 * r2)
        grad_x2 = 2 * (grad_x2_f1 * r1 + grad_x2_f2 * r2)

        # Barzilai-Borwein step size (scalar version)
        alpha_init = 1e-6 # Fallback
        if iter > 1
            s_vec = [x1 - x1_old; x2 - x2_old]
            y_bb_vec = [grad_x1 - grad_x1_old; grad_x2 - grad_x2_old]
            
            y_bb_vec_norm_sq = dot(y_bb_vec, y_bb_vec)
            if y_bb_vec_norm_sq > 1e-9
                s_dot_y = dot(s_vec, y_bb_vec)
                s_dot_s = dot(s_vec, s_vec)
                alpha_bb1 = clamp(abs(s_dot_y) / y_bb_vec_norm_sq, 1e-9, 1.0)
                alpha_bb2 = clamp(s_dot_s / (abs(s_dot_y) + eps()), 1e-9, 1.0)
                alpha_init = min(alpha_bb1, alpha_bb2)
            end
        end

        # Store old values for next iteration's BB step
        x1_old, x2_old = x1, x2
        grad_x1_old, grad_x2_old = grad_x1, grad_x2
        

        # implementing the full backtracking line search for every pixel, which would be very slow.
        # You can add the line search back here if precision is more critical than speed.
        alpha = alpha_init

        # Update estimates
        x1_new = x1 - alpha * grad_x1
        x2_new = x2 - alpha * grad_x2

        # Convergence check
        if abs(x1_new - x1) < tol && abs(x2_new - x2) < tol
             break
        end
        x1, x2 = x1_new, x2_new
    end
    
    return x1, x2
end


"""
    TV_mapping_pixel_by_pixel(y1, y2, LUTs)

Loops through each pixel of the input images `y1` and `y2` and calls
`estimate_single_pixel` to perform the optimization independently for each.
"""
function TV_mapping_pixel_by_pixel(y1, y2, LUTs)
    rows, cols = size(y1)
    x1_est = zeros(rows, cols)
    x2_est = zeros(rows, cols)

    # Use nested loops to iterate through each pixel
    for j in 1:cols
        for i in 1:rows
            # Get the y-values for the current pixel
            y1_pixel = y1[i, j]
            y2_pixel = y2[i, j]

            # Run the estimation for this single pixel
            x1_pixel, x2_pixel = estimate_single_pixel(y1_pixel, y2_pixel, LUTs)
            
            # Store the result in the output maps
            x1_est[i, j] = x1_pixel
            x2_est[i, j] = x2_pixel
        end
    end

    return x1_est, x2_est, nothing # Return nothing for loss for simplicity
end
