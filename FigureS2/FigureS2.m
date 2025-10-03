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
%% Demo code to plot the magnitude and phase of the steady-state of SPGRE with diffusion
%% Takes long to complete whole calculation!

addpath('../tools/');

Gamma = 4285 * 2* pi; %[rad/Gauss]
TR = 10e-3; % TR [s]
T1 = 1000e-3; % T1 [s]
T2 = 100*1e-3; % T2 [s]
flip    = [20];  % Flip angles in degree
inc     = [0.5 1 2 4];  % RF phase increments in degree
Nex     = 600;      % Number of excitations, about 6*TR to reach steady-state
Nf      = 181;      % The number of isochromats.
TG = TR; %Duration for gradient. We use simple model that gradient is applied during TR.

%Diffusion for simulation
Ds = [800:200:3000]*1e-12;

f0_Bloch_D = zeros(length(Ds), length(flip), length(inc));
f0_Analytical_D = zeros(length(Ds), length(flip), length(inc));


%Gradient strength
%This is about 14*pi Moment
G_ave = 7*23.4852;%[Gauss/m]
G = G_ave*Gamma; %[rad/m]



%Calculate signal evolution
for jj=1:length(flip)
    for kk=1:length(inc)
        parfor ii=1:length(Ds)
            D = Ds(ii);
            [Msig,Rfph]=SPGRE_Bloch_diffusion(T1,T2,TR,flip(jj),inc(kk),Nex,Nf,D,TG,G_ave);
            f0_Bloch_D(ii,jj,kk) = mean(Msig((end-10):end)); %Obtain steady-state signal
        end
    end
end

Da = [800:20:4000]*1e-12;
%Calculate analytical signal evolution
for jj=1:length(flip)
    for kk=1:length(inc)
        alpha = flip(jj)*(pi/180); 
        C0 = inc(kk)*(pi/180)*(1/2);
        for ii=1:length(Da)
            D = Da(ii);
            f0_Analytical_D(ii,jj,kk) = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C0, D, G, TG);
        end
    end
end


%%

color_array = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
%Plot magnitude against diffusion coefficients
figure;
subplot(1,2,1);hold on;
legend_text = {};
for ii=1:length(inc)
    plot(Ds*1e+12, abs(f0_Bloch_D(:,1,ii)),'o', 'MarkerEdgeColor', color_array{ii}, 'MarkerSize', 12); 
    legend_text{1,ii} = strcat('\theta=', num2str(inc(ii)), '\circ');
end
for ii=1:length(inc)
    plot(Da*1e+12, abs(f0_Analytical_D(:,1,ii)),'-', 'Color', color_array{ii}, 'LineWidth',2); 
    legend_text{1,length(inc)+ii} = strcat('\theta=', num2str(inc(ii)), '\circ');
end
xlabel('Diffusion coefficient (\mu mm^2/s)'); ylabel('Magnitude (a.u.)');
legend(legend_text); xlim([800 3000]);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);

%Plot phase against diffusion coefficients
subplot(1,2,2);hold on;
legend_text = {};
for ii=1:length(inc)
    plot(Ds*1e+12, (180/pi)*angle(f0_Bloch_D(:,1,ii))+90,'o', 'MarkerEdgeColor', color_array{ii}, 'MarkerSize', 12); 
    legend_text{1,ii} = strcat('\theta=', num2str(inc(ii)), '\circ (Bloch Equation)');
end
for ii=1:length(inc)
    plot(Da*1e+12, (180/pi)*angle(f0_Analytical_D(:,1,ii))+90,'-', 'Color', color_array{ii}, 'LineWidth',2); 
    legend_text{1,length(inc)+ii} = strcat('\theta=', num2str(inc(ii)), '\circ (Analytical)');
end
xlabel('Diffusion coefficient (\mu mm^2/s)'); ylabel('Phase (\circ)');
legend(legend_text); xlim([800 3000]);
set(gca, 'fontname', 'Arial', 'FontSize',18,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,6])