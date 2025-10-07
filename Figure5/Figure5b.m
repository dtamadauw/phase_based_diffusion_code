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


clear;
close all;
clc;

%% 1. Define Realistic Sequence Parameters

% --- PBD Parameters (from paper Table 1) ---
params_PBD.TR = 10.9e-3; % ms
params_PBD.FA_deg = 20; % degrees
params_PBD.ScanTime = 4*60 + 45; % (4:45)
params_PBD.VoxelVol = 1.5 * 1.5 * 4.0; % mm^3
params_PBD.Bandwidth = 325; % Hz
params_PBD.NSlice = 32;

% --- SESE Parameters (from paper Table 1) ---
params_SESE.TR = 500e-3; % 
params_SESE.TE = [15 30 60] * 1e-3;
params_SESE.ScanTime = 56*3; %
params_SESE.VoxelVol = 1.5 * 1.5 * 4.0; % mm^3%
params_SESE.Bandwidth = 244; % Hz
params_SESE.r_factor = 1; 
% --- MESE Parameters (from paper Table 1) ---
params_MESE.TR = 1200e-3; % 
params_MESE.TE = [8.5 17 25 34 42 51 59 67] * 1e-3;
params_MESE.ScanTime = 3*60 + 25; %
params_MESE.VoxelVol = 1.5 * 1.5 * 4.0;%
params_MESE.Bandwidth = 488; % Hz
params_MESE.r_factor = 1; 


% --- Tissue & Simulation Parameters ---
T1 = 1.0;   % T1 in seconds [s]
T2_values = [0.100, 0.120, 0.200]; % T2 in seconds [s] for different tissues
T2_true_s_values = T2_values*1e+3; 
ADC_true_m2s_values = [800, 1300, 1800] * 1e-12; % ADC in m^2/s
N_trials = 100000;       % Number of Monte Carlo trials
SNR_0_range = [2:0.5:50];   % Range of base SNR to test (SNR for a single, unaveraged acquisition with S0=1)
noise = 1./SNR_0_range;
S0_EPI = 1.0;
TEs = [15 30 60]*1e-3;
T2Decay = @(T2)(exp(-TEs./T2));
%% 2. Calculate the "Fair Comparison" Normalization Factors

% Factor 1: Slice-encoding
F_Slice = 1./sqrt(params_PBD.NSlice);

% Factor 2: Voxel Volume
F_Voxel = params_MESE.VoxelVol / params_PBD.VoxelVol;

% Factor 3: Bandwidth (SNR ~ 1/sqrt(BW))
F_BW_SESE = sqrt(params_PBD.Bandwidth / params_SESE.Bandwidth);
F_BW_MESE = sqrt(params_PBD.Bandwidth / params_MESE.Bandwidth);

% Factor 4: Scan Time (SNR ~ sqrt(Time))
F_Time_MESE = sqrt(params_PBD.ScanTime / params_MESE.ScanTime);
F_Time_SESE = sqrt(params_PBD.ScanTime / params_SESE.ScanTime);


F_ASSET = 1/sqrt(params_MESE.r_factor);

% --- Combined Factor (for single-shot comparison) ---
Total_Correction_Factor_SESE =  F_Time_SESE * F_Voxel * F_BW_SESE * F_Slice * F_ASSET;
Total_Correction_Factor_MESE =  F_Time_MESE * F_Voxel * F_BW_MESE * F_Slice * F_ASSET;


fprintf('Effective PBD STD per unit time (SESE) is %.2f times the EPI SNR.\n\n', Total_Correction_Factor_SESE);
fprintf('Effective PBD STD per unit time (MESE) is %.2f times the EPI SNR.\n\n', Total_Correction_Factor_MESE);




%% 4. Main Monte Carlo Loop
% =========================================================================

command = system('julia --threads auto ../tools/julia/main_T2.jl');
command = system(['julia --threads auto ../tools/julia/main_PBT2.jl ',num2str(Total_Correction_Factor_SESE)]);

load(['./MESET2_Julia.mat']);
load(['./PBT2_Julia_' num2str(Total_Correction_Factor_SESE) '.mat']);
%%

results = struct();


for t2_idx = 1:length(ADC_true_m2s_values)
    T2_true_s = T2_values(t2_idx);
    T2_true_plot = T2_true_s * 1000; % ms   




    results(t2_idx).T2_true = T2_true_plot;
    results(t2_idx).MESE = MESET2_jul(:,:,t2_idx)'; 
    results(t2_idx).PBD = 1000*PBT2_jul(:,:,t2_idx)'; % convert to ms
end




%% 5. Analyze and Plot Results with Improved Legend
% =========================================================================
fprintf('5. Analyzing and plotting results...\n');

% --- Define Plotting Styles ---
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]}; % Blue, Red, Yellow
line_styles = {'-', '--'}; % Solid for PBD, Dashed for MESE
method_names = {'PBD', 'SESE'};

figure('Position', [50 50, 1600, 550], 'Color', 'w');

% --- Subplot 1: Mean T2 (Bias) ---
ax1 = subplot(1,3,1); hold on;
for t2_idx = 1:length(T2_true_s_values)
    plot(SNR_0_range, mean(results(t2_idx).PBD, 2),  line_styles{1}, 'Color', colors{t2_idx}, 'LineWidth', 2.5);
    plot(SNR_0_range, mean(results(t2_idx).MESE, 2), line_styles{2}, 'Color', colors{t2_idx}, 'LineWidth', 2);
    yline(results(t2_idx).T2_true, '-', 'Color', [colors{t2_idx} 0.5]);
end
title('Mean T2');xlim([3 50]);ylim([80 350]);
xlabel('SNR of Reference SESE'); ylabel('T2 (ms)');
grid on; box on;
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);


% --- Subplot 2: Median T2 (Robustness to Bias) ---
subplot(1,3,2); hold on;
for t2_idx = 1:length(T2_true_s_values)
    plot(SNR_0_range, median(results(t2_idx).PBD, 2),  line_styles{1}, 'Color', colors{t2_idx}, 'LineWidth', 2.5);
    plot(SNR_0_range, median(results(t2_idx).MESE, 2), line_styles{2}, 'Color', colors{t2_idx}, 'LineWidth', 2);
    yline(results(t2_idx).T2_true, '-', 'Color', [colors{t2_idx} 0.5]);
end
title('Median T2');xlim([3 50]);ylim([80 350]);
xlabel('SNR of Reference SESE');
grid on; box on;
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

% --- Subplot 3: Standard Deviation (Precision) ---
subplot(1,3,3); hold on;
for t2_idx = 1:length(T2_true_s_values)
    plot(SNR_0_range, std(results(t2_idx).PBD, 0, 2),  line_styles{1}, 'Color', colors{t2_idx}, 'LineWidth', 2.5);
    plot(SNR_0_range, std(results(t2_idx).MESE, 0, 2), line_styles{2}, 'Color', colors{t2_idx}, 'LineWidth', 2);
end
title('Standard Deviation');xlim([3 50]);
xlabel('SNR of Reference SESE'); ylabel('Std. Dev. (ms)');
set(gca, 'YScale', 'log');
grid on; box on;
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);


% --- Create Two Separate, Clean Legends ---
% This is a trick to create two legends on the same subplot.
% First legend for methods (line styles)
hold(ax1, 'on');
h_methods = gobjects(length(method_names), 1);
for i = 1:length(method_names)
    h_methods(i) = plot(ax1, NaN, NaN, line_styles{i}, 'Color', 'k', 'LineWidth', 2);
end
lgd1 = legend(h_methods, method_names, 'Location', 'northwest', 'Box', 'off');
legend box off
lgd1.Title.String = '';
hold(ax1, 'off');


% Second legend for T2 values (colors) using a hidden axes overlay
ax_hidden = axes('Position', get(ax1, 'Position'), 'Visible', 'off');
hold(ax_hidden, 'on');
h_t2 = gobjects(length(T2_true_s_values), 1);
legend_labels_t2 = cell(size(T2_true_s_values));
for i = 1:length(T2_true_s_values)
    h_t2(i) = plot(ax_hidden, NaN, NaN, 's', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    legend_labels_t2{i} = sprintf('T2 = %.0f ms', T2_true_s_values(i));
end
lgd2 = legend(h_t2, legend_labels_t2, 'Location', 'northeast', 'Box', 'off');
lgd2.Title.String = '';
hold(ax_hidden, 'off');

set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);


sgtitle(sprintf('Time-Normalized T2 Estimation Performance (Total Time = %.1f min)', T_total/60), 'FontSize', 16, 'FontWeight', 'bold');

disp('T2 simulation complete.');