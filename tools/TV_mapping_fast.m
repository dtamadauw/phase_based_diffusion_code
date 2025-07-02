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


function [x1_est, x2_est, loss] = TV_mapping_fast(sig, lambda1, lambda2, beta1, beta2, LUTs)
    %Sig: Complex images
    %lambda1, lambda2, beta1, beta2: Coefficients for regularization.
    %LUTs: contains spline model which is used for optimization.

    %Scale coefficient:
    mdl_theta = LUTs.mdl_theta;
    mdl_theta2 = LUTs.mdl_theta2;
    y1 = angle(sig(:,:,1,2).*conj(sig(:,:,1,4)))/2+pi/2; %phi_H: phase for high gradient moment
    y2 = angle(sig(:,:,1,1).*conj(sig(:,:,1,3)))/2+pi/2; %phi_L: phase for low gradient moment
    % Initialize parameters
    [rows, cols] = size(y1);
    x1 = 0.2 * ones(rows, cols); % Initial T2 map
    x2 = 0.2 * ones(rows, cols); % Initial ADC map Max:3.5000e-09
    max_iter = 250;
    tol = 1e-6;
    alpha = 0.1; % Step size
    alpha2 = 0.9; % Step size

    %Create mask
    %Noisy areas have random phase, which results in difficulty in
    %convergence.
    mask_sig = abs(sig(:,:,1,1)); th_sig = 0.5*multithresh(mask_sig(:));
    mask_sig = mask_sig > th_sig;
    se = strel('disk',5);
    mask_erode = imclose(mask_sig,se);
    mask_sig = mask_erode;

    mask_weight = abs(sig(:,:,1,1)); mask_mean = mean(mask_weight(mask_sig>0));
    mask_weight = mask_weight./mask_mean;


    %Maskout initial x2
    x1 = x1.*mask_sig; x2 = x2.*mask_sig;
    mask = mask_sig(:);

    % Parameters for BB step size
    alpha1_max = 1.0;
    alpha1_min = 1e-8;%1e-3;
    alpha2_max = 1.0;  % Smaller maximum step size for x2
    alpha2_min = 1e-9;%1e-4; % Smaller minimum step size for x2


    % Initialize previous gradients and steps
    grad_x1_old = zeros(rows, cols);
    grad_x2_old = zeros(rows, cols);
    s1_old = zeros(rows, cols);
    s2_old = zeros(rows, cols);
    x2_old = x2;

    params.mdl_theta = mdl_theta;
    params.mdl_theta2 = mdl_theta2;
    params.delta = 0.5;
    params.mask = mask_weight(:);%mask;

    window_size = 5;  % Window size for moving average
    obj_history = zeros(window_size, 1);
    loss.res1 = [];
    loss.res2 = [];
    alpha2_ar = [];

    for iter = 1:max_iter
        xn = [x1(:) x2(:)];

        % Calculate residuals
        f1 = fnval(mdl_theta.sp, xn');
        f2 = fnval(mdl_theta2.sp, xn');
        r1 = f1' - y1(:); r1 = r1.*params.mask;
        r2 = f2' - y2(:); r2 = r2.*params.mask;

        % Calculate gradients of data fidelity terms
        grad_x1_f1=fnval(mdl_theta.dx1, xn'); grad_x1_f1 = grad_x1_f1';
        grad_x2_f1=fnval(mdl_theta.dx2, xn'); grad_x2_f1 = grad_x2_f1';
        grad_x1_f2=fnval(mdl_theta2.dx1, xn'); grad_x1_f2 = grad_x1_f2';
        grad_x2_f2=fnval(mdl_theta2.dx2, xn'); grad_x2_f2 = grad_x2_f2';

        % Calculate TV gradients        
        grad_tv_x1 = tv_gradient(x1);
        grad_tv_x2 = tv_gradient(x2);

        % Calculate L2 gradients
        grad_l2_x1 = 2 * x1(:);  % Gradient of ||x1||^2
        grad_l2_x2 = 2 * x2(:);  % Gradient of ||x2||^2

        % Total gradient
        grad_x1 = 2 * (grad_x1_f1 .* r1 + grad_x1_f2 .* r2) + lambda1 * grad_tv_x1(:) + beta1 * grad_l2_x1;
        grad_x2 = 2 * (grad_x2_f1 .* r1 + grad_x2_f2 .* r2) + lambda2 * grad_tv_x2(:) + beta2 .* grad_l2_x2;

        % Current objective function value
        grad__abs_x1 = gradient(x1); grad__abs_x2 = gradient(x2);
        grad__abs_x1 = grad__abs_x1(:).*params.mask; grad__abs_x2 = grad__abs_x2(:).*params.mask;
        obj_current = sum(r1(:).^2) + sum(r2(:).^2) + ...
                     lambda1 * sum(abs(grad__abs_x1(:))) + lambda2 * sum(abs(grad__abs_x2(:))) + ...
                     beta1 * sum(x1(:).^2) + beta2 * sum(x2(:).^2);
        % Update objective history
        obj_history = [obj_current; obj_history(1:end-1)];


        %Use BB step as initial guess for backtracking ---
        if iter > 1
            % Calculate combined s and y_bb for a single BB step estimate
            s_vec = [ (x1(:) - x1_old(:)); (x2(:) - x2_old(:)) ];
            y_bb_vec = [ (grad_x1(:) - grad_x1_old(:)); (grad_x2(:) - grad_x2_old(:)) ];
            
            if norm(y_bb_vec) > 0
                alpha_bb1 = max(1e-9, min(1.0, abs(s_vec' * y_bb_vec) / (y_bb_vec' * y_bb_vec)));
                alpha_bb2 = max(1e-9, min(1.0, (s_vec' * s_vec) / abs(s_vec' * y_bb_vec)));
                alpha_init = min(alpha_bb1, alpha_bb2); % Use BB as initial guess
            else
                alpha_init = 1e-7; % Fallback initial guess
            end
        else
            alpha_init = 1e-7; % Initial guess for first iteration
        end


        % Store current values and gradients for next iteration
        x1_old = x1;
        x2_old = x2;
        grad_x1_old = grad_x1;
        grad_x2_old = grad_x2;


        % Update with backtracking line search
        alpha = backtracking_line_search(x1, x2, grad_x1(:), grad_x2(:), alpha_init, ...
                                            y1, y2, lambda1, lambda2, beta1, beta2, params);

        % Update estimates
        x1 = x1 - alpha * reshape(grad_x1,rows,cols);
        x2 = x2 - alpha * reshape(grad_x2,rows,cols);

        % Ensure positivity
        x1 = max(x1, 0);
        x2 = max(x2, 0);

        loss.res1(iter) = obj_current;

        % Check convergence using multiple criteria
        if iter >= window_size
            % Objective function relative change
            obj_change = abs(mean(diff(obj_history(1:window_size)))) / (mean(obj_history(1:window_size)) + eps);            
            loss.res2(iter) = obj_change;
            % Check all convergence criteria
            if obj_change < tol
                break;
            end

        end


    end

    x1_est = x1;
    x2_est = x2;
end


function alpha = backtracking_line_search(x1, x2, grad_x1, grad_x2, alpha_init, ...
    y1, y2, lambda1, lambda2, beta1, beta2, params)
    % Backtracking line search parameters
    beta = 0.5;
    c = 1e-4;
    alpha = alpha_init;
    
    % Initial objective value
    f_curr = objective_function(x1, x2, y1, y2, lambda1, lambda2, beta1, beta2, params);
    
    % Directional derivative
    dir_deriv = sum(grad_x1(:).^2) + sum(grad_x2(:).^2);
    
    % Backtracking
    while true
        % Try step
        x1_new = x1(:) - alpha * grad_x1;
        x2_new = x2(:) - alpha * grad_x2;
        
        % Evaluate new objective
        f_new = objective_function(x1_new, x2_new, y1, y2, lambda1, lambda2, beta1, beta2, params);
        
        % Check Armijo condition
        if f_new <= f_curr - c * alpha * dir_deriv
            break;
        end
        
        % Reduce step size
        alpha = beta * alpha;
        
        % Prevent too small step sizes
        if alpha < 1e-10
            break;
        end
    end
end


function f = objective_function(x1, x2, y1, y2, lambda1, lambda2, beta1, beta2, params)
    % Calculate objective function value

    xn = [x1(:) x2(:)];
    f1 = fnval(params.mdl_theta.sp, xn');
    f2 = fnval(params.mdl_theta2.sp, xn');
    r1 = (f1' - y1(:)).*params.mask;
    r2 = (f2' - y2(:)).*params.mask;    
    
    grad__abs_x1 = gradient(x1); grad__abs_x2 = gradient(x2);
    grad__abs_x1 = grad__abs_x1(:).*params.mask; grad__abs_x2 = grad__abs_x2(:).*params.mask;
    f = sum(r1(:).^2) + sum(r2(:).^2) + ...
        lambda1 * sum(abs(grad__abs_x1(:))) + lambda2 * sum(abs(grad__abs_x2(:))) + ...
        beta1 * sum(x1(:).^2) + beta2 * sum(x2(:).^2);
  
end


function grad = tv_gradient(x)
    epsilon_tv = 1e-8;
    % Calculate gradient of TV norm
    [rows, cols] = size(x);
    
    % Initialize gradient array
    grad = zeros(rows, cols);
    
    % Calculate differences with proper handling of boundaries
    % Horizontal differences
    diff_h = zeros(rows, cols-1);
    diff_h = x(:,2:end) - x(:,1:end-1);
    
    % Vertical differences
    diff_v = zeros(rows-1, cols);
    diff_v = x(2:end,:) - x(1:end-1,:);
    
    % Calculate weights for horizontal differences
    weight_h = 1.0 ./ sqrt(epsilon_tv + diff_h.^2);
    
    % Calculate weights for vertical differences
    weight_v = 1.0 ./ sqrt(epsilon_tv + diff_v.^2);
    
    % Contribute to gradient from horizontal differences
    grad(:,1:end-1) = grad(:,1:end-1) - diff_h .* weight_h;
    grad(:,2:end) = grad(:,2:end) + diff_h .* weight_h;
    
    % Contribute to gradient from vertical differences
    grad(1:end-1,:) = grad(1:end-1,:) - diff_v .* weight_v;
    grad(2:end,:) = grad(2:end,:) + diff_v .* weight_v;
end
