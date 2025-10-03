%TITLE: MATLAB code for "Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging"
%Matlab script to generate figures using the closed-form equation proposed in the paper
%Author: Daiki Tamada
%Affiliation: Department of Radiology, University of Wisconsin-Madison
%Date: 7/1/2025
%Email: dtamada@wisc.edu

%By downloading, installing, or otherwise accessing or using the Software , you (“Recipient”) agree to receive and use the above-
%identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
%obligations and restrictions you are not permitted to download, install,
%execute, access, or use the SOFTWARE:

%1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE (“Provider’s Rights”) to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation (“WARF”). 
%2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider’s Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient’s responsibility. 
%3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider’s proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
%4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
%5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN “AS IS” BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
%6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER’S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient’s use of the SOFTWARE.
%7.	No use of names/trademarks.  Recipient shall not use Provider’s name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
%8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient’s possession, custody, or control.	

function [x1_est, x2_est, loss] = TV_mapping_fast_B1(sig, lambda1, lambda2, beta1, beta2, LUTs, beta)
    % sig: Complex images [rows, cols, 1, 4]
    % lambda1, lambda2, beta1, beta2: Regularization coefficients.
    % LUTs: Struct containing 3D spline models.
    % beta_or_FA: B1+ map (in %) or absolute FA map (in degrees).

    % Extract measured phases for high and low gradient moments from input signals
    y_high = angle(sig(:,:,1,1) .* conj(sig(:,:,1,3))) / 2 + pi / 2;
    y_low = angle(sig(:,:,1,2) .* conj(sig(:,:,1,4))) / 2 + pi / 2;

    % Initialize parameters
    [rows, cols] = size(y_low);
    x1 = 0.2 * ones(rows, cols); % Initial normalized T2 map
    x2 = 0.2 * ones(rows, cols); % Initial normalized ADC map
    max_iter = 150; % Reduced max iterations due to faster convergence
    tol = 1e-6;
    
    % --- B1+ Correction: Compute per-voxel Flip Angle Map ---
    FAq = LUTs.FA0 * (beta / 100.0);

    % Clamp FA values to the LUT range to prevent extrapolation errors
    FAq = min(max(FAq, LUTs.FAmin), LUTs.FAmax);
    
    % --- Mask Creation (same as original) ---
    mask_sig = abs(sig(:,:,1,1)); 
    th_sig = 0.5 * multithresh(mask_sig(:));
    mask_sig = mask_sig > th_sig;
    se = strel('disk', 5);
    mask_sig = imclose(mask_sig, se);
    x1 = x1 .* mask_sig; 
    x2 = x2 .* mask_sig;
    
    mask_weight = abs(sig(:,:,1,1));
    mask_mean = mean(mask_weight(mask_sig > 0));
    mask_weight = mask_weight ./ mask_mean;
    
    % Store necessary items in a params struct for helper functions
    params.mask = mask_weight(:);
    params.FA_vec = FAq(:)'; % Pass flattened FA map

    % --- Optimization Initialization ---
    grad_x1_old = zeros(rows * cols, 1);
    grad_x2_old = zeros(rows * cols, 1);
    x1_old = x1;
    x2_old = x2;
    
    window_size = 5;
    obj_history = zeros(window_size, 1);
    loss.res1 = [];
    loss.res2 = [];

    for iter = 1:max_iter
        % Flatten unknowns for vectorized spline evaluation
        xn1 = x1(:)';
        xn2 = x2(:)';

        % --- Vectorized Spline and Gradient Evaluation ---
        eval_points = [xn1; xn2; params.FA_vec];
        fL = fnval(LUTs.mdl_theta3D.sp, eval_points)';
        fH = fnval(LUTs.mdl_theta23D.sp, eval_points)';
        
        dT2_fL = fnval(LUTs.mdl_theta3D.dx1, eval_points)';
        dD_fL  = fnval(LUTs.mdl_theta3D.dx2, eval_points)';
        dT2_fH = fnval(LUTs.mdl_theta23D.dx1, eval_points)';
        dD_fH  = fnval(LUTs.mdl_theta23D.dx2, eval_points)';

        % Calculate residuals (masked)
        rL = (fL - y_low(:)) .* params.mask;
        rH = (fH - y_high(:)) .* params.mask;
        
        % --- Total Gradient Calculation ---
        grad_data_x1 = 2 * (dT2_fL .* rL + dT2_fH .* rH);
        grad_data_x2 = 2 * (dD_fL .* rL + dD_fH .* rH);
        
        grad_tv_x1 = tv_gradient(x1);
        grad_tv_x2 = tv_gradient(x2);
        
        grad_x1 = grad_data_x1 + lambda1 * grad_tv_x1(:) + (2 * beta1 * x1(:));
        grad_x2 = grad_data_x2 + lambda2 * grad_tv_x2(:) + (2 * beta2 * x2(:));

        % --- Barzilai-Borwein Step Size ---
        if iter > 1
            s_vec = [(x1(:) - x1_old(:)); (x2(:) - x2_old(:))];
            y_bb_vec = [(grad_x1 - grad_x1_old); (grad_x2 - grad_x2_old)];
            
            if norm(y_bb_vec) > 1e-9
                alpha_bb1 = abs(s_vec' * y_bb_vec) / (y_bb_vec' * y_bb_vec);
                alpha_bb2 = (s_vec' * s_vec) / abs(s_vec' * y_bb_vec);
                alpha_init = min(max(1e-9, alpha_bb1), max(1e-9, alpha_bb2));
                alpha_init = min(1.0, alpha_init); % Cap initial step
            else
                alpha_init = 1e-7;
            end
        else
            alpha_init = 1e-7;
        end

        x1_old = x1(:); x2_old = x2(:);
        grad_x1_old = grad_x1; grad_x2_old = grad_x2;

        % --- Backtracking Line Search ---
        alpha = backtracking_line_search(x1, x2, grad_x1, grad_x2, alpha_init, ...
                                         y_low, y_high, lambda1, lambda2, beta1, beta2, LUTs, params);

        % --- Update Estimates ---
        x1 = x1 - alpha * reshape(grad_x1, rows, cols);
        x2 = x2 - alpha * reshape(grad_x2, rows, cols);

        % Enforce positivity and apply mask
        x1 = max(x1, 0) .* mask_sig; 
        x2 = max(x2, 0) .* mask_sig;

        % --- Check Convergence ---
        obj_current = objective_function(x1, x2, y_low, y_high, lambda1, lambda2, beta1, beta2, LUTs, params);
        obj_history = [obj_current; obj_history(1:end-1)];
        loss.res1(iter) = obj_current;
        
        if iter >= window_size
            obj_change = abs(mean(diff(obj_history))) / (mean(obj_history) + eps);
            loss.res2(iter) = obj_change;
            if obj_change < tol
                %fprintf('Converged at iteration %d\n', iter);
                break;
            end
        end
    end

    x1_est = x1;
    x2_est = x2;
end

% --- Helper Functions for Optimization ---

function alpha = backtracking_line_search(x1, x2, grad_x1, grad_x2, alpha_init, ...
    y_low, y_high, lambda1, lambda2, beta1, beta2, LUTs, params)
    
    beta_bt = 0.5; c = 1e-4; alpha = alpha_init;
    f_curr = objective_function(x1, x2, y_low, y_high, lambda1, lambda2, beta1, beta2, LUTs, params);
    dir_deriv = sum(grad_x1(:).^2) + sum(grad_x2(:).^2);
    
    while true
        x1_new = x1(:) - alpha * grad_x1(:);
        x2_new = x2(:) - alpha * grad_x2(:);
        f_new = objective_function(reshape(x1_new, size(x1)), reshape(x2_new, size(x2)), ...
                                   y_low, y_high, lambda1, lambda2, beta1, beta2, LUTs, params);
        if f_new <= f_curr - c * alpha * dir_deriv
            break;
        end
        alpha = beta_bt * alpha;
        if alpha < 1e-10
            break;
        end
    end
end

function f = objective_function(x1, x2, y_low, y_high, lambda1, lambda2, beta1, beta2, LUTs, params)
    % Calculate objective function value with efficient TV term
    xn1 = x1(:)'; xn2 = x2(:)';
    eval_points = [xn1; xn2; params.FA_vec];
    
    fL = fnval(LUTs.mdl_theta3D.sp, eval_points)';
    fH = fnval(LUTs.mdl_theta23D.sp, eval_points)';
    
    rL = (fL - y_low(:)) .* params.mask;
    rH = (fH - y_high(:)) .* params.mask;
    
    % Efficient TV calculation using forward differences
    tv_x1 = tv_value(x1, params.mask);
    tv_x2 = tv_value(x2, params.mask);
    
    f = sum(rL.^2) + sum(rH.^2) + ...
        lambda1 * tv_x1 + lambda2 * tv_x2 + ...
        beta1 * sum((x1(:).*params.mask).^2) + beta2 * sum((x2(:).*params.mask).^2);
end

function grad = tv_gradient(x)
    % Same as original implementation
    epsilon_tv = 1e-8;
    [rows, cols] = size(x);
    grad = zeros(rows, cols);
    
    diff_h = x(:,2:end) - x(:,1:end-1);
    diff_v = x(2:end,:) - x(1:end-1,:);
    
    weight_h = 1.0 ./ sqrt(epsilon_tv + diff_h.^2);
    weight_v = 1.0 ./ sqrt(epsilon_tv + diff_v.^2);
    
    grad(:,1:end-1) = grad(:,1:end-1) - diff_h .* weight_h;
    grad(:,2:end) = grad(:,2:end) + diff_h .* weight_h;
    
    grad(1:end-1,:) = grad(1:end-1,:) - diff_v .* weight_v;
    grad(2:end,:) = grad(2:end,:) + diff_v .* weight_v;
end

function val = tv_value(x, mask)
    % Efficient anisotropic TV calculation for objective function
    epsilon_tv = 1e-8;
    [rows, cols] = size(x);
    mask_2d = reshape(mask, rows, cols);
    
    dx = diff(x, 1, 2);
    dy = diff(x, 1, 1);
    
    % Apply mask to differences to only include valid regions
    mask_dx = mask_2d(:, 1:end-1) > 0 & mask_2d(:, 2:end) > 0;
    mask_dy = mask_2d(1:end-1, :) > 0 & mask_2d(2:end, :) > 0;

    val = sum(sqrt(dx(mask_dx).^2 + epsilon_tv)) + sum(sqrt(dy(mask_dy).^2 + epsilon_tv));
end
