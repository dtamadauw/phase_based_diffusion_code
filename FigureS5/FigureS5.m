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


addpath('../tools/');


%%
%Monte-Carlo simulation for phase-based ADC mapping
%using PBD equation


TR = 10.9e-3;
numRO = 256;
Gamma = 4285 * 2* pi; %[rad/Gauss/s]
FOV = 0.256;
BW = 4*pi*2*(125)*1e+3/numRO;%[Rad/px]
Resolution = FOV/numRO;%[m/px]
Tsamp = 1/BW;
G = 2*pi/(Gamma*Resolution)*(1/TR);
G_2pi = G*Gamma;
SNR = [5];
noise = 1./SNR;
nex = 100*100;

disp(['Gradiet Moment (2*pi): ' num2str(G) ' (Gauss/m)']);
rng(321);%Fixed seed

lambdas = [0 1e-5 0.5e-4 1e-4 0.5e-3 1e-3 0.5e-2 1e-2 1e-1];
betas =  [0 1e-5 0.5e-4 1e-4 0.5e-3 1e-3 0.5e-2 1e-2 1e-1];


%T1 dependency
dphi = 2;
FA = 20;
C = (pi/180)*dphi/2;
alpha = FA*(pi/180);

T1 = 1000*1e-3;
%https://www.sciencedirect.com/science/article/pii/S1076633218301880
T2 = [100 120 200]*1e-3;
%Ds = [500 1000 1500 2000]*1e-12;%[100:100:2000]*1e-12;
Ds = [800 1300 1800]*1e-12;

G1 = G_2pi;%4.6367e+05;
G2 = 7*G_2pi;%4.3054e+06;

y1p = zeros(nex, length(noise), length(Ds));
y2p = zeros(nex, length(noise), length(Ds));
y1n = zeros(nex, length(noise), length(Ds));
y2n = zeros(nex, length(noise), length(Ds));

ys = zeros(1,4);
noise_dit = zeros(nex,length(noise),4);

for kk=1:length(Ds)

    [ys(1),~,~,~] = analytical_SPGR_W_diffusion(TR, T1, T2(kk), alpha, C, Ds(kk), G1, TR);
    [ys(2),~,~,~] = analytical_SPGR_W_diffusion(TR, T1, T2(kk), alpha, C, Ds(kk), G2, TR);
    [ys(3),~,~,~] = analytical_SPGR_W_diffusion(TR, T1, T2(kk), alpha, -C, Ds(kk), G1, TR);
    [ys(4),~,~,~] = analytical_SPGR_W_diffusion(TR, T1, T2(kk), alpha, -C, Ds(kk), G2, TR);
    ys = ys./abs(ys(1));


    for ii=1:length(noise)
        n_s = zeros(nex,4);
        for jj=1:10
            rng(42+jj+kk);%Reduce seed-dependency
            ind_seed = ((jj-1)*nex/10 + 1):(jj*nex/10);
            n_s(ind_seed,:) = abs(ys(1))*complex(normrnd(0,noise(ii),nex/10,4),normrnd(0,noise(ii),nex/10,4))/sqrt(2);
        end
        noise_dit(:,ii,:) = n_s;
        y1p(:,ii,kk) = ys(1)+n_s(:,1); y2p(:,ii,kk) = ys(2)+n_s(:,2);
        y1n(:,ii,kk) = ys(3)+n_s(:,3); y2n(:,ii,kk) = ys(4)+n_s(:,4);
    end

end

F0 = angle((y1n).*conj((y1p)))/2+pi/2;
F1 = angle((y2n).*conj((y2p)))/2+pi/2;

PH0 = angle(y1n.*conj(y1p))/2+pi/2;
PH1 = angle(y2n.*conj(y2p))/2+pi/2;


params.TR = TR;
params.G1 = TR*1e+4*G1/Gamma;
params.G2 = TR*1e+4*G2/Gamma;
params.te1 = 0;
params.FA = FA;
params.dphi = dphi;
params.opuser8 = 0;

[LUTs] = build_LUTs_ROA(params);
T2_PBD_lambda=[];D_PBD_lambda=[];
T2_PBD_beta=[];D_PBD_beta=[];
PH0_2D = reshape(PH0,sqrt(size(PH0,1)),sqrt(size(PH0,1)), size(PH0,2)*size(PH0,3));
PH1_2D = reshape(PH1,sqrt(size(PH1,1)),sqrt(size(PH1,1)), size(PH1,2)*size(PH1,3));
for jj=1:length(lambdas)
    for ii = 1:size(PH0_2D,3)
        [T2_PBD_lambda(:,:,ii,jj), D_PBD_lambda(:,:,ii,jj), loss] = TV_mapping_penalty(PH0_2D(:,:,ii), PH1_2D(:,:,ii), lambdas(jj), lambdas(jj), 0, 0, LUTs);
    end
end

for jj=1:length(betas)
    for ii = 1:size(PH0_2D,3)
        [T2_PBD_beta(:,:,ii,jj), D_PBD_beta(:,:,ii,jj), loss] = TV_mapping_penalty(PH0_2D(:,:,ii), PH1_2D(:,:,ii), 0, 0, betas(jj), betas(jj), LUTs);
    end
end

T2_PBD_lambda = reshape(T2_PBD_lambda,nex,size(PH0_2D,3),length(betas));
D_PBD_lambda = reshape(D_PBD_lambda,nex,size(PH0_2D,3),length(betas));
T2_PBD_beta = reshape(T2_PBD_beta,nex,size(PH0_2D,3),length(betas));
D_PBD_beta = reshape(D_PBD_beta,nex,size(PH0_2D,3),length(betas));


T2_PBD_lambda_u = T2_PBD_lambda*LUTs.x1*1e+3;
D_PBD_lambda_u = D_PBD_lambda*LUTs.x2*1e+12;
T2_PBD_beta_u = T2_PBD_beta*LUTs.x1*1e+3;
D_PBD_beta_u = D_PBD_beta*LUTs.x2*1e+12;





%% --- 3. Create the 2x3 Figure ---
Ds_true = Ds*1e+12;

T2s_true = T2*1e+3;


% --- Define plot properties for consistency ---
colors = {'#0072BD', '#D95319', '#77AC30'}; % Blue, Orange, Green
line_styles = {'-', ':', '-'};
line_width = 4;
font_size = 12;
SNRt = 1./noise;
% --- Top Row: ADC Results ---

figure('Position', [50 50, 120, 550], 'Color', 'w');
% Panel A: Mean ADC
ax1 = subplot(1,2,1); 
hold on;
Ds_unit = Ds.*1e+12;
temp = squeeze(mean(D_PBD_lambda_u(:,:,:),1))';temp = 100.*temp./Ds_unit-100;
temp_median = squeeze(median(D_PBD_lambda_u(:,:,:),1))';temp_median = 100.*temp_median./Ds_unit-100;
for i = 1:length(Ds)
    plot(lambdas,temp(:,i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
    plot(lambdas,temp_median(:,i), 'Color', colors{i}, 'LineStyle', line_styles{2}, 'LineWidth', line_width);
end
yline(0, '--', 'Color', 'k'); % Zero bias line
hold off;
%title('Mean/Median Bias', 'FontSize', font_size+2, 'FontWeight', 'bold');
xlabel('\lambda', 'FontSize', font_size, 'FontWeight', 'bold');
ylabel('ADC Bias (%)', 'FontSize', font_size, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);set(gca, 'XScale', 'log');


% Panel C: Standard Deviation of ADC
ax2 = subplot(1,2,2); 
hold on;
Ds_unit = Ds.*1e+12;
temp = squeeze(std(D_PBD_lambda_u(:,:,:),1))';
for i = 1:length(Ds)
    plot(lambdas,temp(:,i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
end
hold off;
ylabel('Standard Deviation (\mumm^2/s)', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('\lambda', 'FontSize', font_size, 'FontWeight', 'bold');
%title('Standard Deviation', 'FontSize', font_size+2, 'FontWeight', 'bold');
grid on;
box on;
legend({'ADC=800\mumm^2/s', 'ADC=1300\mumm^2/s', 'ADC=1800\mumm^2/s'});legend boxoff;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,15,5]);set(gca, 'XScale', 'log');


% This is a trick to create two legends on the same subplot.
% First legend for methods (line styles)
method_names = {'Mean', 'Median'};

hold(ax1, 'on');
h_methods = gobjects(length(method_names), 1);
for i = 1:length(method_names)
    h_methods(i) = plot(ax1, NaN, NaN, line_styles{i}, 'Color', 'k', 'LineWidth', 2);
end
lgd1 = legend(h_methods, method_names, 'Location', 'northwest', 'Box', 'off');
legend box off
lgd1.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax1, 'off');


% Second legend for T2 values (colors) using a hidden axes overlay
ax_hidden = axes('Position', get(ax1, 'Position'), 'Visible', 'off');
hold(ax_hidden, 'on');
h_t2 = gobjects(length(T2s_true), 1);
legend_labels_t2 = cell(size(T2s_true));
for i = 1:length(T2s_true)
    h_t2(i) = plot(ax_hidden, NaN, NaN, 's', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    legend_labels_t2{i} = sprintf('ADC = %.0f mumm^2/s', Ds(i)*1e+12);
end
lgd2 = legend(h_t2, legend_labels_t2, 'Location', 'northeast', 'Box', 'off');
lgd2.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax_hidden, 'off');

set(gcf,'units','inches','position',[0,0,16,5]);set(gca, 'XScale', 'log');


%%
%% --- 3. Create the 2x3 Figure ---
%Beta
Ds_true = Ds*1e+12;

T2s_true = T2*1e+3;


% --- Define plot properties for consistency ---
colors = {'#0072BD', '#D95319', '#77AC30'}; % Blue, Orange, Green
line_styles = {'-', ':', '-'};
line_width = 4;
font_size = 12;
SNRt = 1./noise;
% --- Top Row: ADC Results ---

figure('Position', [50 50, 120, 550], 'Color', 'w');
% Panel A: Mean ADC
ax1 = subplot(1,2,1); 
hold on;
Ds_unit = Ds.*1e+12;
temp = squeeze(mean(D_PBD_beta_u(:,:,:),1))';temp = 100.*temp./Ds_unit-100;
temp_median = squeeze(median(D_PBD_beta_u(:,:,:),1))';temp_median = 100.*temp_median./Ds_unit-100;
for i = 1:length(Ds)
    plot(betas,temp(:,i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
    plot(betas,temp_median(:,i), 'Color', colors{i}, 'LineStyle', line_styles{2}, 'LineWidth', line_width);
end
yline(0, '--', 'Color', 'k'); % Zero bias line
hold off;
%title('Mean/Median Bias', 'FontSize', font_size+2, 'FontWeight', 'bold');
xlabel('\beta', 'FontSize', font_size, 'FontWeight', 'bold');
ylabel('ADC Bias (%)', 'FontSize', font_size, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);set(gca, 'XScale', 'log');


% Panel C: Standard Deviation of ADC
ax2 = subplot(1,2,2); 
hold on;
Ds_unit = Ds.*1e+12;
temp = squeeze(std(D_PBD_beta_u(:,:,:),1))';
for i = 1:length(Ds)
    plot(betas,temp(:,i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
end
hold off;
ylabel('Standard Deviation (\mumm^2/s)', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('\beta', 'FontSize', font_size, 'FontWeight', 'bold');
%title('Standard Deviation', 'FontSize', font_size+2, 'FontWeight', 'bold');
grid on;
box on;
legend({'ADC=800\mumm^2/s', 'ADC=1300\mumm^2/s', 'ADC=1800\mumm^2/s'});legend boxoff;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,15,5]);set(gca, 'XScale', 'log');

method_names = {'Mean', 'Median'};

hold(ax1, 'on');
h_methods = gobjects(length(method_names), 1);
for i = 1:length(method_names)
    h_methods(i) = plot(ax1, NaN, NaN, line_styles{i}, 'Color', 'k', 'LineWidth', 2);
end
lgd1 = legend(h_methods, method_names, 'Location', 'northwest', 'Box', 'off');
legend box off
lgd1.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax1, 'off');


% Second legend for T2 values (colors) using a hidden axes overlay
ax_hidden = axes('Position', get(ax1, 'Position'), 'Visible', 'off');
hold(ax_hidden, 'on');
h_t2 = gobjects(length(T2s_true), 1);
legend_labels_t2 = cell(size(T2s_true));
for i = 1:length(T2s_true)
    h_t2(i) = plot(ax_hidden, NaN, NaN, 's', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    legend_labels_t2{i} = sprintf('ADC = %.0f mumm^2/s', Ds(i)*1e+12);
end
lgd2 = legend(h_t2, legend_labels_t2, 'Location', 'northeast', 'Box', 'off');
lgd2.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax_hidden, 'off');

set(gcf,'units','inches','position',[0,0,16,5]);set(gca, 'XScale', 'log');















%%
% --- Bottom Row: T2 Results ---

% Panel D: Mean T2
figure
ax1 = subplot(1,2,1); hold on;
T2s_unit = T2s_true;
temp = squeeze(mean(T2_PBD_lambda_u(:,:,:),1))';temp = 100.*temp./T2s_unit-100;
temp_median = squeeze(median(T2_PBD_lambda_u(:,:,:),1))';temp_median = 100.*temp_median./T2s_unit-100;
for i = 1:length(T2)
    plot(lambdas, temp(:, i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
    plot(lambdas, temp_median(:, i), 'Color', colors{i}, 'LineStyle', line_styles{2}, 'LineWidth', line_width);
end
yline(0, '--', 'Color', 'k'); % Zero bias line
hold off;
ylabel('T2 Bias (%)', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('\lambda', 'FontSize', font_size, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);set(gca, 'XScale', 'log');

% Panel F: Standard Deviation of T2
ax2 = subplot(1,2,2); hold on;
temp = squeeze(std(T2_PBD_lambda_u(:,:,:),1))';
for i = 1:length(T2)
    plot(lambdas, temp(:, i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
end
hold off;
ylabel('Standard Deviation (ms)', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('\lambda', 'FontSize', font_size, 'FontWeight', 'bold');
grid on;
box on;
legend({'T2=100ms', 'T2=120ms', 'T2=200ms'});legend boxoff;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,16,5]);set(gca, 'XScale', 'log');


% This is a trick to create two legends on the same subplot.
% First legend for methods (line styles)
method_names = {'Mean', 'Median'};

hold(ax1, 'on');
h_methods = gobjects(length(method_names), 1);
for i = 1:length(method_names)
    h_methods(i) = plot(ax1, NaN, NaN, line_styles{i}, 'Color', 'k', 'LineWidth', 2);
end
lgd1 = legend(h_methods, method_names, 'Location', 'northwest', 'Box', 'off');
legend box off
lgd1.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax1, 'off');


% Second legend for T2 values (colors) using a hidden axes overlay
ax_hidden = axes('Position', get(ax1, 'Position'), 'Visible', 'off');
hold(ax_hidden, 'on');
h_t2 = gobjects(length(T2s_true), 1);
legend_labels_t2 = cell(size(T2s_true));
for i = 1:length(T2s_true)
    h_t2(i) = plot(ax_hidden, NaN, NaN, 's', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    legend_labels_t2{i} = sprintf('T2 = %.0f ms', T2s_true(i));
end
lgd2 = legend(h_t2, legend_labels_t2, 'Location', 'northeast', 'Box', 'off');
lgd2.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax_hidden, 'off');

set(gcf,'units','inches','position',[0,0,16,5]);set(gca, 'XScale', 'log');


%%
% --- Bottom Row: T2 Results ---

% Panel D: Mean T2, Beta
figure
ax1 = subplot(1,2,1); hold on;
T2s_unit = T2s_true;
temp = squeeze(mean(T2_PBD_beta_u(:,:,:),1))';temp = 100.*temp./T2s_unit-100;
temp_median = squeeze(median(T2_PBD_beta_u(:,:,:),1))';temp_median = 100.*temp_median./T2s_unit-100;
for i = 1:length(T2)
    plot(lambdas, temp(:, i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
    plot(lambdas, temp_median(:, i), 'Color', colors{i}, 'LineStyle', line_styles{2}, 'LineWidth', line_width);
end
yline(0, '--', 'Color', 'k'); % Zero bias line
hold off;
ylabel('T2 Bias (%)', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('\beta', 'FontSize', font_size, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);set(gca, 'XScale', 'log');

% Panel F: Standard Deviation of T2
ax2 = subplot(1,2,2); hold on;
temp = squeeze(std(T2_PBD_beta_u(:,:,:),1))';
for i = 1:length(T2)
    plot(lambdas, temp(:, i), 'Color', colors{i}, 'LineStyle', line_styles{1}, 'LineWidth', line_width);
end
hold off;
ylabel('Standard Deviation (ms)', 'FontSize', font_size, 'FontWeight', 'bold');
xlabel('\beta', 'FontSize', font_size, 'FontWeight', 'bold');
grid on;
box on;
legend({'T2=100ms', 'T2=120ms', 'T2=200ms'});legend boxoff;
set(gca, 'FontSize', font_size);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,16,5]);set(gca, 'XScale', 'log');

% This is a trick to create two legends on the same subplot.
% First legend for methods (line styles)
method_names = {'Mean', 'Median'};

hold(ax1, 'on');
h_methods = gobjects(length(method_names), 1);
for i = 1:length(method_names)
    h_methods(i) = plot(ax1, NaN, NaN, line_styles{i}, 'Color', 'k', 'LineWidth', 2);
end
lgd1 = legend(h_methods, method_names, 'Location', 'northwest', 'Box', 'off');
legend box off
lgd1.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax1, 'off');


% Second legend for T2 values (colors) using a hidden axes overlay
ax_hidden = axes('Position', get(ax1, 'Position'), 'Visible', 'off');
hold(ax_hidden, 'on');
h_t2 = gobjects(length(T2s_true), 1);
legend_labels_t2 = cell(size(T2s_true));
for i = 1:length(T2s_true)
    h_t2(i) = plot(ax_hidden, NaN, NaN, 's', 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    legend_labels_t2{i} = sprintf('T2 = %.0f ms', T2s_true(i));
end
lgd2 = legend(h_t2, legend_labels_t2, 'Location', 'northeast', 'Box', 'off');
lgd2.Title.String = '';
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
hold(ax_hidden, 'off');

set(gcf,'units','inches','position',[0,0,16,5]);set(gca, 'XScale', 'log');

