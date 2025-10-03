%TITLE: MATLAB code for "Quantitative Diffusion and T2 Mapping Using RF-Modulated Phase-Based Gradient Echo Imaging"
%Matlab script to generate figures using the closed-form equation proposed in the paper
%Author: Daiki Tamada
%Affiliation: Department of Radiology, University of Wisconsin-Madison
%Date: 10/1/2025
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



clear; close all; clc;
addpath('../tools/')
%% --- 1. Simulation Setup ---

TR = 10.9e-3;      % Repetition Time [s]
T1 = 1.0;        % T1 value [s]
T2_true = 120e-3;% True T2 [s]
D_true = 1.3e-9; % True Diffusion (ADC) [m^2/s]

% --- Gradient Parameters ---
Gamma = 4285 * 2* pi; %[rad/Gauss/s]
FOV = 0.256;
numRO = 256;
Resolution = FOV/numRO;%[m/px]
G = 2*pi/(Gamma*Resolution)*(1/TR);
G_2pi = G*Gamma;

G_L = G_2pi;
T_grad = TR;

% --- Define the 4D Parameter Sweep Ranges ---
theta_deg_vec = [1 2];
% Define fixed noise levels corresponding to a reference SNR of 10 and 40

%Calculate base signal for noise magnitude.
[f0_base, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, deg2rad(20), 0.5*deg2rad(20), D_true, G_L, T_grad);

% for a typical signal magnitude of 0.08.
%sigma_n_vec     = [0.08/5, 0.08/20]; 
SNR_n_vec = [5 10];
alpha_deg_vec   = [10:0.5:30];
G_H_ratio_vec   = [2:0.5:10];

% --- Pre-allocate result arrays for all metrics ---
dims = [length(theta_deg_vec), length(SNR_n_vec), length(alpha_deg_vec), length(G_H_ratio_vec)];
std_ADC_realistic = zeros(dims);
std_T2_realistic  = zeros(dims);
std_ADC_const_noise = zeros(dims);
std_T2_const_noise  = zeros(dims);
ratio_ADC         = zeros(dims);
ratio_T2          = zeros(dims);
condition_number  = zeros(dims);
correlation       = zeros(dims);

%% --- 2. Perform the 4D Parameter Sweep ---

fprintf('Starting detailed 4D CRLB analysis...\n');
tic;

for i_theta = 1:length(theta_deg_vec)
    theta_deg = theta_deg_vec(i_theta);
    theta_rad = 0.5*deg2rad(theta_deg);
    
    for i_sigma = 1:length(SNR_n_vec)
       
        ref_SNR = SNR_n_vec(i_sigma);
        fprintf('  Running for theta = %.1f deg (Ref SNR ~%d)\n', theta_deg, round(ref_SNR));
        
        for i_alpha = 1:length(alpha_deg_vec)
            alpha_deg = alpha_deg_vec(i_alpha);
            alpha_rad = deg2rad(alpha_deg);
            
            for i_gh = 1:length(G_H_ratio_vec)
                G_H_ratio = G_H_ratio_vec(i_gh);
                G_H = G_L * G_H_ratio;
                
                % Define a reference phase variance for the constant noise analysis
                
                var_phi_ref_single_acq = 1 / (2 * ref_SNR^2);
                var_phi_ref_final = var_phi_ref_single_acq / 2;
                Sigma_const_noise = [var_phi_ref_final, 0; 0, var_phi_ref_final];

                % --- Jacobian Calculation (same for all analyses) ---
                dD = D_true * 1e-2;
                dT2 = T2_true * 1e-2;
                
                [f0_L_p, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad, theta_rad, D_true + dD, G_L, T_grad);
                [f0_L_m, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad, theta_rad, D_true - dD, G_L, T_grad);
                J_11 = (angle(f0_L_p) - angle(f0_L_m)) / (2 * dD);
                
                [f0_H_p, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad, theta_rad, D_true + dD, G_H, T_grad);
                [f0_H_m, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad, theta_rad, D_true - dD, G_H, T_grad);
                J_21 = (angle(f0_H_p) - angle(f0_H_m)) / (2 * dD);
                
                [f0_L_p, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true + dT2, alpha_rad, theta_rad, D_true, G_L, T_grad);
                [f0_L_m, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true - dT2, alpha_rad, theta_rad, D_true, G_L, T_grad);
                J_12 = (angle(f0_L_p) - angle(f0_L_m)) / (2 * dT2);
                
                [f0_H_p, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true + dT2, alpha_rad, theta_rad, D_true, G_H, T_grad);
                [f0_H_m, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true - dT2, alpha_rad, theta_rad, D_true, G_H, T_grad);
                J_22 = (angle(f0_H_p) - angle(f0_H_m)) / (2 * dT2);

                %[FIM, CRLBt, out] = pbd_phase_fim(TR, T1, T2_true, alpha_rad, theta_rad, D_true, G_L, G_H, T_grad, 0.1);

                J = [J_11, J_12; J_21, J_22];
                
                % --- Analysis 1: Realistic Noise Model (Signal-Dependent) ---
                [S_L_p, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad,  theta_rad, D_true, G_L, T_grad);
                [S_L_m, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad, -theta_rad, D_true, G_L, T_grad);
                [S_H_p, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad,  theta_rad, D_true, G_H, T_grad);
                [S_H_m, ~] = analytical_SPGR_W_diffusion(TR, T1, T2_true, alpha_rad, -theta_rad, D_true, G_H, T_grad);

                sigma_n = abs(S_L_p)/ref_SNR;

                var_phi_L_p = 1 / (2 * (abs(S_L_p) / sigma_n)^2);
                var_phi_L_m = 1 / (2 * (abs(S_L_m) / sigma_n)^2);
                var_phi_H_p = 1 / (2 * (abs(S_H_p) / sigma_n)^2);
                var_phi_H_m = 1 / (2 * (abs(S_H_m) / sigma_n)^2);
                
                var_phi_L_final = (var_phi_L_p + var_phi_L_m) / 4;
                var_phi_H_final = (var_phi_H_p + var_phi_H_m) / 4;
                
                Sigma_realistic = [var_phi_L_final, 0; 0, var_phi_H_final];
                FIM_realistic = J' * (Sigma_realistic \ J);
                
                % --- Analysis 2: Conditioning-Only (Constant Reference Noise) ---
                FIM_const_noise = J' * (Sigma_const_noise \ J);
                
                % --- Calculate all metrics ---
                if rcond(FIM_realistic) > 1e-16 && rcond(FIM_const_noise) > 1e-16
                    CRLB_realistic = inv((FIM_realistic));
                    CRLB_const_noise = inv(FIM_const_noise);
                    
                    % Store realistic precision
                    std_ADC_realistic(i_theta, i_sigma, i_alpha, i_gh) = sqrt(CRLB_realistic(1, 1));
                    std_T2_realistic(i_theta, i_sigma, i_alpha, i_gh)  = sqrt(CRLB_realistic(2, 2));
                    
                    % Store conditioning-only precision
                    std_ADC_const_noise(i_theta, i_sigma, i_alpha, i_gh) = sqrt(CRLB_const_noise(1, 1));
                    std_T2_const_noise(i_theta, i_sigma, i_alpha, i_gh)  = sqrt(CRLB_const_noise(2, 2));
                    
                    % Store meaningful, unitless ratios
                    ratio_ADC(i_theta, i_sigma, i_alpha, i_gh) = std_ADC_realistic(i_theta, i_sigma, i_alpha, i_gh) / std_ADC_const_noise(i_theta, i_sigma, i_alpha, i_gh);
                    ratio_T2(i_theta, i_sigma, i_alpha, i_gh)  = std_T2_realistic(i_theta, i_sigma, i_alpha, i_gh) / std_T2_const_noise(i_theta, i_sigma, i_alpha, i_gh);
                    
                    % Store condition number and correlation
                    condition_number(i_theta, i_sigma, i_alpha, i_gh) = cond(FIM_realistic);
                    correlation(i_theta, i_sigma, i_alpha, i_gh) = -CRLB_realistic(1, 2) / sqrt(CRLB_realistic(1, 1) * CRLB_realistic(2, 2));
                else
                    % Mark ill-conditioned points as Inf
                    [std_ADC_realistic(i_theta, i_sigma, i_alpha, i_gh), std_T2_realistic(i_theta, i_sigma, i_alpha, i_gh), ...
                     std_ADC_const_noise(i_theta, i_sigma, i_alpha, i_gh), std_T2_const_noise(i_theta, i_sigma, i_alpha, i_gh), ...
                     ratio_ADC(i_theta, i_sigma, i_alpha, i_gh), ratio_T2(i_theta, i_sigma, i_alpha, i_gh), ...
                     condition_number(i_theta, i_sigma, i_alpha, i_gh), correlation(i_theta, i_sigma, i_alpha, i_gh)] = deal(inf);
                end
            end
        end
    end
end
toc;
fprintf('Simulation finished.\n\n');

%% --- 3. Visualize the Results ---

% --- FIGURE FOR MAIN PAPER ---
hFig_Main = figure('Name', 'CRLB Analysis: Main Paper Figure', 'Position', [50, 50, 800, 600]);
sgtitle('CRLB Analysis: Decomposing ADC Variance Sources (\theta = 3^{\circ})', 'FontSize', 16, 'FontWeight', 'bold');

% Select data for theta = 2 degrees (index 2)
i_theta_main = 2; 
main_adc_metrics = {std_ADC_realistic*1e12, std_ADC_const_noise*1e12, ratio_ADC};
main_adc_titles = {'Realistic std(ADC) (um^2/ms)', 'Constant std(ADC) (um^2/ms)', 'Ratio (Realistic / Constant Noise)'};
main_adc_cmaps = {'parula', 'parula', 'jet'};

num_rows = length(main_adc_metrics);
num_cols = length(SNR_n_vec);

for i_metric = 1:num_rows
    metric_data = main_adc_metrics{i_metric};
    
    % Find min/max for color scaling based on the selected theta
    slice_for_clim = squeeze(metric_data(i_theta_main, :, :, :));
    data_min = prctile(slice_for_clim(isfinite(slice_for_clim(:))), 2);

    if i_metric < 3
        data_max = 3000;
        data_min = 0;
    else
        data_max = 1.5;
        data_min = 1;
    end


    for i_sigma = 1:length(SNR_n_vec)
        plot_idx = (i_metric - 1) * num_cols + i_sigma;
        ax = subplot(num_rows, num_cols, plot_idx);
        
        data_slice = squeeze(metric_data(i_theta_main, i_sigma, :, :))';
        imagesc(alpha_deg_vec, G_H_ratio_vec, data_slice);
        set(gca, 'YDir', 'normal');
        
        if i_metric == 1
            title_str = sprintf('Ref SNR ~%d', round(SNR_n_vec(i_sigma)));
            title(title_str);
        end
        if i_sigma == 1
            ylabel(main_adc_titles{i_metric}, 'FontWeight', 'bold');
        end
        if i_metric == num_rows
             xlabel('Flip Angle (deg)');
        end
        
        colormap(ax, main_adc_cmaps{i_metric});
        colorbar;
        caxis([data_min, data_max]);
        set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
    end
end


set(gcf,'units','inches','position',[0,0,10,12]);

%%

% --- FIGURE FOR MAIN PAPER (T2) ---
hFig_Main = figure('Name', 'CRLB Analysis: Main Paper Figure', 'Position', [50, 50, 800, 600]);
sgtitle('CRLB Analysis: Decomposing T2 Variance Sources (\theta = 2^{\circ})', 'FontSize', 16, 'FontWeight', 'bold');

% Select data for theta = 2 degrees (index 2)
i_theta_main = 2; 
main_t2_metrics = {std_T2_realistic*1e3, std_T2_const_noise*1e3, ratio_T2};
main_t2_titles = {'Realistic std(T2) (ms)', 'Constant std(T2) (ms)', 'Ratio (Realistic / Constant Noise)'};
main_t2_cmaps = {'parula', 'parula', 'jet'};

num_rows = length(main_t2_metrics);
num_cols = length(SNR_n_vec);

for i_metric = 1:num_rows
    metric_data = main_t2_metrics{i_metric};
    
    % Find min/max for color scaling based on the selected theta
    slice_for_clim = squeeze(metric_data(i_theta_main, :, :, :));

    if i_metric < 3
        data_max = 40;
        data_min = 15;
    else
        data_min = 1.0;
        data_max = 1.5;
    end

    for i_sigma = 1:length(SNR_n_vec)
        plot_idx = (i_metric - 1) * num_cols + i_sigma;
        ax = subplot(num_rows, num_cols, plot_idx);
        
        data_slice = squeeze(metric_data(i_theta_main, i_sigma, :, :))';
        imagesc(alpha_deg_vec, G_H_ratio_vec, data_slice);
        set(gca, 'YDir', 'normal');
        
        if i_metric == 1
            title_str = sprintf('Ref SNR ~%d', round(SNR_n_vec(i_sigma)));
            title(title_str);
        end
        if i_sigma == 1
            ylabel(main_t2_titles{i_metric}, 'FontWeight', 'bold');
        end
        if i_metric == num_rows
             xlabel('Flip Angle (deg)');
        end
        
        colormap(ax, main_t2_cmaps{i_metric});
        colorbar;
        caxis([data_min, data_max]);
        set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
    end
end


set(gcf,'units','inches','position',[0,0,10,12]);
