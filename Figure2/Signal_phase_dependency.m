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

%%
addpath('../tools/');

%Setting up pparameters
TR = 10e-3;
numRO = 256;
Gamma = 4285 * 2* pi;%[rad/Gauss/s]
FOV = 0.256;
BW = 4*pi*2*(125)*1e+3/numRO;%[Rad/px]
Resolution = FOV/numRO;%[m/px]
Tsamp = 1/BW;
G = 2*pi/(Gamma*Resolution)*(1/TR);
G_2pi = G*Gamma;

disp(['Gradiet Moment (2*pi): ' num2str(G) ' (Gauss/m)']);


%%
%T2 dependency

dphis = [1 2 4];
alpha = 20*(pi/180);

T1s = [1000 700  1300]*1e-3;
T2s = [1:10:300]*1e-3%
Ds = [1000*1e-12];%


Phase_D1 = zeros(length(T2s),length(dphis),length(T1s));
Phase_D2 = zeros(length(T2s),length(dphis),length(T1s));
G1 = 4.6367e+05;
G2 = 4.3054e+06;

for kk = 1:length(T1s)
    for ii = 1:length(dphis)
        for jj=1:length(T2s)
    
                dphi = dphis(ii);
                C = (pi/180)*dphi/2;
                T2 = T2s(jj);
                T1 = T1s(kk);
                D = Ds(1);
    
                [y1, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C, D, G1, TR);
                [y2, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C, D, G2, TR);
    
                Phase_D1(jj,ii,kk) = angle(y1) + pi/2;
                Phase_D2(jj,ii,kk) = angle(y2) + pi/2;
    
        end
    end
end


% --- Plotting and Dual Legend Creation ---

cmap = colororder();
LineStyle_lst = {'-', '--', ':'};
h_for_theta_legend = gobjects(length(dphis), 1); % Handles for the theta legend

figure;
hold on;

% Plot all the lines
for kk=1:length(T1s)
    % Plot lines for each dphi value, cycling through colors
    p1 = plot(T2s*1000, Phase_D2(:,1,kk)*(180/pi),'LineWidth',2,'Color',cmap(1,:),'LineStyle',LineStyle_lst{kk});
    p2 = plot(T2s*1000, Phase_D2(:,2,kk)*(180/pi),'LineWidth',2,'Color',cmap(2,:),'LineStyle',LineStyle_lst{kk});
    p3 = plot(T2s*1000, Phase_D2(:,3,kk)*(180/pi),'LineWidth',2,'Color',cmap(3,:),'LineStyle',LineStyle_lst{kk});
    set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
    % On the first loop (kk=1), store the handles of the three colored lines.
    % These will serve as the representative items for the 'theta' legend.
    if kk == 1
        h_for_theta_legend = [p1, p2, p3];
    end
end
grid on;
% --- Finalize Plot Formatting ---
xlabel(gca, 'T2 (ms)');
ylabel(gca, 'Signal Phase (\circ)');
%title(gca, 'Moment = 14\pi');
ylim(gca, [0 50]);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);
box on


% --- Create the two legends ---

% 1. Create the first legend for Theta (which corresponds to line color)
legend1 = legend(h_for_theta_legend, {'\theta = 1\circ', '\theta = 2\circ', '\theta = 4\circ'});
set(legend1, 'Location', 'northeast', 'Box', 'off', 'FontSize', 18);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

% 2. Create a new, invisible axes that sits on top of the current one.
% This new axes will hold our second legend.
ax_new = axes('Position', get(gca,'Position'), 'Visible', 'off');

% 3. Create "dummy" line plots on the new axes for the T1 legend.
% These plots are not visible because their data is NaN, but they hold the
% line style information we need for the legend.
h_for_T1_legend = gobjects(length(T1s), 1);
for kk = 1:length(T1s)
    h_for_T1_legend(kk) = line(NaN, NaN, 'Parent', ax_new, ...
                               'LineStyle', LineStyle_lst{kk}, ...
                               'Color', 'k', ...
                               'LineWidth', 2);
end
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

% Create legend strings dynamically from the T1s values
T1_legend_strs = cellfun(@(x) sprintf('T1 = %03d ms', x), num2cell(int16(T1s*1000)), 'UniformOutput', false);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

% 4. Create the second legend for T1 (line styles) on the invisible axes
legend2 = legend(h_for_T1_legend, T1_legend_strs);
set(legend2, 'Location', 'northwest', 'Box', 'off', 'FontSize', 18);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

hold off;


% legend('boxoff');
% xlabel('T2 (ms)'); ylabel('Signal Phase (\circ)');title('Moment = 14\pi');ylim([0 50]);
% set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','bold','LineWidth',2);
% set(gcf,'units','inches','position',[0,0,5,5])


%%
%ADC dependency
dphis = [1 2 4];
alpha = 20*(pi/180);

T1s = [1000 700 1300]*1e-3;
T2 = 100*1e-3
Ds = [100:100:2000]*1e-12;


Phase_D1 = zeros(length(Ds),length(dphis));
Phase_D2 = zeros(length(Ds),length(dphis));
G1 = 4.6367e+05;
G2 = 4.3054e+06;

for kk = 1:length(T1s)
    for ii = 1:length(dphis)
        for jj=1:length(Ds)
    
                C = (pi/180)*dphis(ii)/2;
                D = Ds(jj);
                T1 = T1s(kk);
    
                [y1, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C, D, G1, TR);
                [y2, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C, D, G2, TR);
    
                Phase_D1(jj,ii,kk) = angle(y1) + pi/2;
                Phase_D2(jj,ii,kk) = angle(y2) + pi/2;
    
        end
    end
end


cmap = colororder();
LineStyle_lst = {'-', '--', ':'};

figure;
for kk = 1:length(T1s)
plot(Ds*1e+12, Phase_D2(:,1,kk)*(180/pi),'LineWidth',2,'Color',cmap(1,:),'LineStyle',LineStyle_lst{kk});hold on;
plot(Ds*1e+12, Phase_D2(:,2,kk)*(180/pi),'LineWidth',2,'Color',cmap(2,:),'LineStyle',LineStyle_lst{kk});
plot(Ds*1e+12, Phase_D2(:,3,kk)*(180/pi),'LineWidth',2,'Color',cmap(3,:),'LineStyle',LineStyle_lst{kk});
end
grid on;

legend({'\theta = 1\circ', '\theta = 2\circ', '\theta = 4\circ'});legend('boxoff');
xlabel('ADC (\mu mm^2/s)'); ylabel('Signal Phase (\circ)');ylim([0 50]);xlim([100 2000]);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5])



%%
%FA dependency
dphis = [1 2 4];
FAs = [0:40];


T1s = [1000 700 1300]*1e-3;
T2 = 100*1e-3;
D = 1000*1e-12;


Phase_D1 = zeros(length(Ds),length(dphis));
Phase_D2 = zeros(length(Ds),length(dphis));
G1 = 4.6367e+05;
G2 = 4.3054e+06;

for kk = 1:length(T1s)
    for ii = 1:length(dphis)
        for jj=1:length(FAs)
    
                alpha = FAs(jj)*(pi/180);
                T1 = T1s(kk);
                C = (pi/180)*dphis(ii)/2;
    
                [y1, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C, D, G1, TR);
                [y2, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C, D, G2, TR);
    
                Phase_D1(jj,ii,kk) = angle(y1) + pi/2;
                Phase_D2(jj,ii,kk) = angle(y2) + pi/2;
    
        end
    end
end

cmap = colororder();
LineStyle_lst = {'-', '--', ':'};

figure;
for kk = 1:length(T1s)
plot(FAs, Phase_D2(:,1,kk)*(180/pi),'LineWidth',2,'Color',cmap(1,:),'LineStyle',LineStyle_lst{kk});hold on;
plot(FAs, Phase_D2(:,2,kk)*(180/pi),'LineWidth',2,'Color',cmap(2,:),'LineStyle',LineStyle_lst{kk});
plot(FAs, Phase_D2(:,3,kk)*(180/pi),'LineWidth',2,'Color',cmap(3,:),'LineStyle',LineStyle_lst{kk});
end
grid on;

legend({'\theta = 1\circ', '\theta = 2\circ', '\theta = 4\circ'}); legend('boxoff');
xlabel('FA (\circ)'); ylabel('Signal Phase (\circ)');ylim([0 50]);xlim([0 40]);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5])



