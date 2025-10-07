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

% --- SS-EPI Parameters (from paper Table 1) ---
params_EPI.TR = 5000e-3; % ms
params_EPI.TE = 82e-3; % ms
params_EPI.ScanTime_1dir = 1*60 + 0; % 60 seconds 
params_EPI.ScanTime_3dir = 2*60 + 0; % 150 seconds
params_EPI.VoxelVol = 1.5 * 1.5 * 4.0; % mm^3
params_EPI.Bandwidth = 1953; % Hz
params_EPI.b_values = [0, 1000]*1e6; % s/mm^2, matching PBD's effective b-value
params_EPI.nex_b_high = 3; 
params_EPI.r_factor = 2; 


% --- Tissue & Simulation Parameters ---
T1 = 1.0;   % T1 in seconds [s]
T2_values = [0.100, 0.120, 0.200]; % T2 in seconds [s] for different tissues
ADC_true_m2s_values = [800, 1300, 1800] * 1e-12; % ADC in m^2/s
N_trials = 100000;       % Number of Monte Carlo trials
SNR_0_range = [2:0.5:50];   % Range of base SNR to test (SNR for a single, unaveraged acquisition with S0=1)
S0_EPI = 1.0;

%% 2. Calculate the "Fair Comparison" Normalization Factors

% Factor 1: Slice-encoding
F_Slice = 1./sqrt(params_PBD.NSlice);

% Factor 2: Voxel Volume
F_Voxel = params_EPI.VoxelVol / params_PBD.VoxelVol;

% Factor 3: Bandwidth (SNR ~ 1/sqrt(BW))
F_BW = sqrt(params_PBD.Bandwidth / params_EPI.Bandwidth);

% Factor 4: Scan Time (SNR ~ sqrt(Time))
F_Time_1dir = sqrt(params_PBD.ScanTime / params_EPI.ScanTime_1dir);
F_Time_3dir = sqrt(params_PBD.ScanTime/ params_EPI.ScanTime_3dir);

F_ASSET = 1/sqrt(params_EPI.r_factor);

% --- Combined Factor (for single-shot comparison) ---
Total_Correction_Factor =  F_Time_1dir * F_Voxel * F_BW * F_Slice * F_ASSET;

fprintf('--- Fair Comparison Factors (EPI relative to PBD) ---\n');
fprintf('Encoding Factor (F_Slice):         %.2f\n', F_Slice);
fprintf('Voxel Volume Ratio (F_Voxel):    %.2f\n', F_Voxel);
fprintf('Bandwidth Factor (F_BW):         %.2f\n', F_BW);
fprintf('Scan Time Factor 1Dir (F_Time):       %.2f\n', F_Time_1dir);
fprintf('Scan Time Factor 3 Dir (F_Time):       %.2f\n', F_Time_3dir);
fprintf('PI Factor 3 Dir (F_Time):       %.2f\n', F_ASSET);
fprintf('---------------------------------------------------\n');
fprintf('Effective PBD per unit time (1dir) is %.2f times the EPI SNR.\n\n', Total_Correction_Factor);




%% 4. Main Monte Carlo Loop
%Run Julia code for faster processing

command = system(['julia --threads auto ../tools/julia/main_PBD.jl ',num2str(Total_Correction_Factor)]);
load(['./PBD_Julia_' num2str(Total_Correction_Factor) '.mat']);

%%
% =========================================================================
results = struct();
fprintf('4. Running Monte Carlo\n');

for adc_idx = 1:length(ADC_true_m2s_values)
    ADC_true_m2s = ADC_true_m2s_values(adc_idx);
    T2 = T2_values(adc_idx);
    ADC_true_plot = ADC_true_m2s * 1e12;
   
    % --- Generate Noiseless Signals ---
    S_EPI_b_low = S0_EPI * exp(-params_EPI.TE/T2) * (1.0-exp(-params_EPI.TR/1000e-3));
    S_EPI_b_high = S0_EPI * exp(-params_EPI.TE/T2) * exp(-params_EPI.b_values(2) * ADC_true_m2s) * (1.0-exp(-params_EPI.TR/1000e-3));


    mag_low = zeros(length(SNR_0_range),N_trials);
    mag_high_1dir = zeros(length(SNR_0_range),N_trials);
    mag_low_3dir = zeros(length(SNR_0_range),N_trials);
    mag_high_3dir = zeros(length(SNR_0_range),N_trials);
    noise_epi_b_low = zeros(length(SNR_0_range),N_trials);
    noise_epi_b_high = zeros(length(SNR_0_range),N_trials);

    for i=1:length(SNR_0_range)
        noise_std_epi = abs(S_EPI_b_low)/SNR_0_range(i);
        noise_epi_b_low(i,:) = noise_std_epi.*(randn(1,N_trials) + 1i*randn(1,N_trials))./ sqrt(2);
        noise_epi_b_high(i,:) = noise_std_epi.*(randn(1,N_trials) + 1i*randn(1,N_trials))./ sqrt(2);
    end

    mag_low = abs(S_EPI_b_low + noise_epi_b_low);
    mag_high_1dir = abs(S_EPI_b_high + noise_epi_b_high/sqrt(params_EPI.nex_b_high));
    mag_high_3dir = abs(S_EPI_b_high + noise_epi_b_high/sqrt((3*params_EPI.nex_b_high)));


    est_adc = @(mag_low, mag_high, b_val) -log(max(1e-9, mag_high) ./ max(1e-9, mag_low)) ./ b_val;

    ADC_estimates_EPI_1dir = est_adc(mag_low, mag_high_1dir, params_EPI.b_values(2));
    ADC_estimates_EPI_3dir = est_adc(mag_low, mag_high_3dir, params_EPI.b_values(2));


    results(adc_idx).ADC_true = ADC_true_plot;
    results(adc_idx).EPI_3dir = ADC_estimates_EPI_3dir * 1e12;
    results(adc_idx).EPI_1dir = ADC_estimates_EPI_1dir * 1e12;
    results(adc_idx).PBD = PBD_jul(:,:,adc_idx)'  * 1e12;
end


for adc_idx = 1:length(ADC_true_m2s_values)

    std_pbd_1 = std(results(adc_idx).PBD, 0, 2);
    std_epi_1 = std(results(adc_idx).EPI_1dir,0,2);
    std_epi_3 = std(results(adc_idx).EPI_3dir,0,2);

    fprintf('ADC: %f, std(PBD/EPI_1dir): %f, std(PBD/EPI_3dir): %f \n',...
        ADC_true_m2s_values(adc_idx)*1e+12, mean(std_pbd_1./std_epi_1), mean(std_pbd_1./std_epi_3));

end




%% 5. Analyze and Plot Results with Improved Legend
% =========================================================================
fprintf('5. Analyzing and plotting results...\n');

% --- Define Plotting Styles ---
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]}; % Blue, Red, Yellow for ADC values
line_styles = {'-', '--', ':'}; % Solid for PBD, Dashed for EPI 3-dir, Dotted for EPI 1-dir
method_names = {'PBD', 'EPI 3-dir', 'EPI 1-dir'};

figure('Position', [50 50 1600 550], 'Color', 'w');

% --- Subplot 1: Mean ADC (Bias) ---
subplot(1,3,1); hold on;
for adc_idx = 1:length(ADC_true_m2s_values)
    plot(SNR_0_range, mean(results(adc_idx).PBD, 2),      line_styles{1}, 'Color', colors{adc_idx}, 'LineWidth', 2.5);
    plot(SNR_0_range, mean(results(adc_idx).EPI_3dir, 2), line_styles{2}, 'Color', colors{adc_idx}, 'LineWidth', 2);
    plot(SNR_0_range, mean(results(adc_idx).EPI_1dir, 2), line_styles{3}, 'Color', colors{adc_idx}, 'LineWidth', 1.5);
    yline(results(adc_idx).ADC_true, '-', 'Color', [colors{adc_idx} 0.5]);
end
title('Mean ADC');xlim([3 50]);ylim([700 2100]);
xlabel('SNR of Reference SS-EPI'); ylabel('ADC (\mum^2/s)');
grid on; box on;
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

% --- Subplot 2: Median ADC (Robustness to Bias) ---
subplot(1,3,2); hold on;
for adc_idx = 1:length(ADC_true_m2s_values)
    plot(SNR_0_range, median(results(adc_idx).PBD, 2),      line_styles{1}, 'Color', colors{adc_idx}, 'LineWidth', 2.5);
    plot(SNR_0_range, median(results(adc_idx).EPI_3dir, 2), line_styles{2}, 'Color', colors{adc_idx}, 'LineWidth', 2);
    plot(SNR_0_range, median(results(adc_idx).EPI_1dir, 2), line_styles{3}, 'Color', colors{adc_idx}, 'LineWidth', 1.5);
    yline(results(adc_idx).ADC_true, '-', 'Color', [colors{adc_idx} 0.5]);
end
title('Median ADC');xlim([3 50]);ylim([700 2100]);
xlabel('SNR of Reference SS-EPI');
grid on; box on;
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

% --- Subplot 3: Standard Deviation (Precision) ---
subplot(1,3,3); hold on;
for adc_idx = 1:length(ADC_true_m2s_values)
    plot(SNR_0_range, std(results(adc_idx).PBD, 0, 2),      line_styles{1}, 'Color', colors{adc_idx}, 'LineWidth', 2.5);
    plot(SNR_0_range, std(results(adc_idx).EPI_3dir, 0, 2), line_styles{2}, 'Color', colors{adc_idx}, 'LineWidth', 2);
    plot(SNR_0_range, std(results(adc_idx).EPI_1dir, 0, 2), line_styles{3}, 'Color', colors{adc_idx}, 'LineWidth', 1.5);
end
title('Standard Deviation'); xlim([3 50]);
xlabel('SNR of Reference SS-EPI'); ylabel('Std. Dev. (\mum^2/s)');
set(gca, 'YScale', 'log');
grid on; box on;
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
% Create invisible dummy plots to generate the legend handles
h_dummy = gobjects(length(method_names) + length(ADC_true_m2s_values), 1);
ax_legend = subplot(1,3,1); % Attach legend to the first subplot
hold(ax_legend, 'on');
for i = 1:length(method_names)
    h_dummy(i) = plot(ax_legend, NaN, NaN, line_styles{i}, 'Color', 'k', 'LineWidth', 2);
end
legend_labels_methods = method_names;

for i = 1:length(ADC_true_m2s_values)
    idx = length(method_names) + i;
    h_dummy(idx) = plot(ax_legend, NaN, NaN, 's', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', colors{i}, 'MarkerSize', 10);
    legend_labels_adc{i} = sprintf('ADC = %.0f', ADC_true_m2s_values(i)*1e12);
end
hold(ax_legend, 'off');

% Display the combined legend
lgd = legend(h_dummy, [legend_labels_methods, legend_labels_adc]);
lgd.NumColumns = 2;
lgd.Location = 'northeast';
legend box off
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

disp('Simulation complete.');
