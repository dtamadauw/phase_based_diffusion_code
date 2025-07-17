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


function [LUTs] = build_LUTs_ROA(params, varargin)

if length(varargin) > 0
    T1 = varargin{1}*1e-3;
else
    T1 = 1000*1e-3;
end



T2s = [20:20:300]*1e-3;
Ds = [0:250:3000]*1e-12;
T2s_re = [20:1:300]*1e-3;
Ds_re = [0:10:3000]*1e-12;
dphis = [0];


normalize_sig = true;
if normalize_sig
    x1_norm = max(T2s);
    x2_norm = max(Ds);
else
    x1_norm = 1;
    x2_norm = 1;
end


area = [params.G1 params.G2];
G_ave = 100*area/(params.TR*1e+6); 
Gamma = 4285 * 2* pi; %[rad/Gauss]
G = G_ave*Gamma; %[rad/m]


TE = params.te1*1e-6;
TR = params.TR;
FAs1 = (params.FA/90)*ones(1,500);
alpha = [params.FA*(pi/180)];

LUTs_mag = zeros(length(Ds), length(T2s), length(dphis));
LUTs_magSm = zeros(length(Ds), length(T2s), length(dphis));
LUTs_DESS = zeros(length(Ds), length(T2s), length(dphis));
LUTs_DIFF = zeros(length(Ds), length(T2s), length(dphis));
LUTs_dS = zeros(length(Ds), length(T2s), length(dphis));
LUTs_eta = zeros(length(Ds), length(T2s), length(dphis));
LUTs_eta2 = zeros(length(Ds), length(T2s), length(dphis));
LUTs_theta = zeros(length(Ds), length(T2s), length(dphis));
LUTs_theta2 = zeros(length(Ds), length(T2s), length(dphis));
LUTs_epsilon_eta0 = zeros(length(Ds), length(T2s), length(dphis));
LUTs_epsilon_eta1 = zeros(length(Ds), length(T2s), length(dphis));
LUTs_index = zeros(length(Ds), length(T2s), length(dphis),3);
LUTs_dict = zeros(length(Ds), length(T2s), length(dphis),4);

D_ind = 1;


disp(['T1: ' num2str(T1)]);


if params.opuser8 == 0
    isROA1 = 0;
    isROA2 = 0;
elseif params.opuser8 == 1
    isROA1 = 2;
    isROA2 = 2;
elseif params.opuser8 == 2
    isROA1 = 2;
    isROA2 = 2;
elseif params.opuser8 == 3
    isROA1 = 4;
    isROA2 = 4;
elseif params.opuser8 == 4
    isROA1 = 4;
    isROA2 = 4;
elseif params.opuser8 == 5
    isROA1 = 1;
    isROA2 = 1;
end


parfor D_ind = 1:length(Ds)
    
    LUTs_mag_temp = zeros(1,length(T2s), length(dphis));
    LUTs_magSm_temp = zeros(1,length(T2s), length(dphis));
    LUTs_DESS_temp = zeros(1,length(T2s), length(dphis));
    LUTs_DIFF_temp = zeros(1,length(T2s), length(dphis));
    LUTs_dS_temp = zeros(1,length(T2s), length(dphis));
    LUTs_eta_temp = zeros(1,length(T2s), length(dphis));
    LUTs_eta2_temp = zeros(1,length(T2s), length(dphis));
    LUTs_theta_temp = zeros(1,length(T2s), length(dphis));
    LUTs_theta2_temp = zeros(1,length(T2s), length(dphis));
    LUTs_epsilon_eta0_temp = zeros(1,length(T2s), length(dphis));
    LUTs_epsilon_eta1_temp = zeros(1,length(T2s), length(dphis));
    LUTs_index_temp = zeros(1,length(T2s), length(dphis),3);
    LUTs_dict_temp = zeros(1,length(T2s), length(dphis),4);

    dphi_ind = 1;
    for dphi = dphis   

        C0 = params.dphi*(pi/180)*(1/2);
        C = [(dphi+params.dphi) -1.0*(dphi+params.dphi)] * (pi/180)*(1/2);
    
        T2_ind = 1;
        for T2 = T2s
            [y0, ys0, epsilon_eta0] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha(1), C0, Ds(D_ind), G(1), TR);
            [y1, ys1, epsilon_eta1] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha(1), C(1), Ds(D_ind), G(2), TR);
            [y1r, ys1r] = analytical_SPGR_W_diffusion(TR, T1, T2, alpha(1), C(2), Ds(D_ind), G(2), TR);
    
            y_G1 = y0;%mean((y0((end-10):end)));%
            ys_G1 = ys0;%mean((ys0((end-10):end)));%
            y_G2 = y1;%mean((y1((end-10):end)));%
            ys_G2 = ys1;%mean((ys1((end-10):end)));%
            yr_G2 = y1r;%mean((y1((end-10):end)));%
            ysr_G2 = ys1r;%mean((ys1((end-10):end)));%

            real_G2 = abs(y_G2-yr_G2)/2;
            imag_G2 = abs(y_G2+yr_G2)/2;

            LUTs_mag_temp(1,T2_ind,dphi_ind) = imag(y_G2);
            LUTs_magSm_temp(1,T2_ind,dphi_ind) = imag(y_G1);
            LUTs_DESS_temp(1,T2_ind,dphi_ind) = abs(ys_G1)/abs(y_G1);
            LUTs_DIFF_temp(1,T2_ind,dphi_ind) = abs(ys_G2)/abs(y_G2);
            LUTs_dS_temp(1,T2_ind,dphi_ind) = (real(y_G2)/real(y_G1)) * (imag(y_G2)/imag(y_G1)); %1 - (ys_G2*y_G1/(ys_G1*y_G2));
            LUTs_eta_temp(1,T2_ind,dphi_ind) = abs(real_G2/real(y_G1));%(real(y_G2)/real(y_G1))./ (abs(y_G2)/abs(y_G1));%
            LUTs_eta2_temp(1,T2_ind,dphi_ind) = abs(imag_G2/imag(y_G1));%%1 - (real(y_G2)*abs(y_G1)/(real(y_G1)*abs(y_G2)));%%imag(y_G2)/imag(y_G1);%imag(y_G2)/imag(y_G1);%1 - (real(y_G2)*abs(y_G1)/(real(y_G1)*abs(y_G2)));%
            LUTs_theta_temp(1,T2_ind,dphi_ind) = abs(angle(y_G1));
            LUTs_theta2_temp(1,T2_ind,dphi_ind) = abs(angle(y_G2.*conj(yr_G2))/2-pi/2);
            LUTs_epsilon_eta0_temp(1,T2_ind,dphi_ind) = epsilon_eta0;
            LUTs_epsilon_eta1_temp(1,T2_ind,dphi_ind) = epsilon_eta1;
            LUTs_index_temp(1,T2_ind,dphi_ind,1) = Ds(D_ind)/x2_norm;
            LUTs_index_temp(1,T2_ind,dphi_ind,2) = T2/x1_norm;
            LUTs_index_temp(1,T2_ind,dphi_ind,3) = dphi;
            
            y1s = complex(imag(y_G2), real(y_G2));
            y0s = complex(imag(y_G1), real(y_G1));
            norm_y = [y1s y0s conj(y1s) conj(y0s)];% norm_y = norm_y./norm(norm_y);
            LUTs_dict_temp(1,T2_ind,dphi_ind,:) = norm_y;
    
            T2_ind = T2_ind + 1;
        end

        dphi_ind = dphi_ind + 1;
    end

    LUTs_mag(D_ind, :, :) = LUTs_mag_temp;
    LUTs_magSm(D_ind, :, :) = LUTs_magSm_temp;
    LUTs_DESS(D_ind, :, :) = LUTs_DESS_temp;
    LUTs_DIFF(D_ind, :, :) = LUTs_DIFF_temp;
    LUTs_dS(D_ind, :, :) = LUTs_dS_temp;
    LUTs_eta(D_ind, :, :) = LUTs_eta_temp;
    LUTs_eta2(D_ind, :, :) = LUTs_eta2_temp;
    LUTs_theta(D_ind, :, :) = LUTs_theta_temp;
    LUTs_index(D_ind, :, :, :) = LUTs_index_temp;
    LUTs_theta2(D_ind, :, :, :) = LUTs_theta2_temp;
    LUTs_epsilon_eta0(D_ind, :, :, :) = LUTs_epsilon_eta0_temp;
    LUTs_epsilon_eta1(D_ind, :, :, :) = LUTs_epsilon_eta1_temp;
    LUTs_dict(D_ind, :, :, :) = LUTs_dict_temp;
end


if normalize_sig
    [xt2,xd] = meshgrid(T2s,Ds);
    xx = [xt2(:) xd(:)*(1e+8)];
    mdl_theta = build_models(T2s./x1_norm, Ds./x2_norm, LUTs_theta');
    mdl_theta2 = build_models(T2s./x1_norm, Ds./x2_norm, LUTs_theta2');
else
    [xt2,xd] = meshgrid(T2s,Ds);
    xx = [xt2(:) xd(:)*(1e+8)];
    mdl_theta = build_models(T2s, Ds*(1e+8), LUTs_theta');
    mdl_theta2 = build_models(T2s, Ds*(1e+8), LUTs_theta2');
end




LUTs_mag_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_magSm_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_DESS_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_DIFF_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_dS_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_eta_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_eta2_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_theta_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_theta2_re = zeros(length(Ds_re), length(T2s_re), length(dphis));
LUTs_index_re = zeros(length(Ds_re), length(T2s_re), length(dphis),3);
LUTs_dict_re = zeros(length(Ds_re), length(T2s_re), length(dphis),4);

[xx,yy] = meshgrid(T2s,Ds);
[xxi,yyi] = meshgrid(T2s_re,Ds_re);
for dd = 1:length(dphis)
    LUTs_mag_re(:,:,dd) = interp2(xx, yy, LUTs_mag(:,:,dd), xxi, yyi);
    LUTs_magSm_re(:,:,dd) = interp2(xx, yy, LUTs_magSm(:,:,dd), xxi, yyi);
    LUTs_DESS_re(:,:,dd) = interp2(xx, yy, LUTs_DESS(:,:,dd), xxi, yyi);
    LUTs_DIFF_re(:,:,dd) = interp2(xx, yy, LUTs_DIFF(:,:,dd), xxi, yyi);
    LUTs_dS_re(:,:,dd) = interp2(xx, yy, LUTs_dS(:,:,dd), xxi, yyi);
    LUTs_eta_re(:,:,dd) = interp2(xx, yy, LUTs_eta(:,:,dd), xxi, yyi);
    LUTs_eta2_re(:,:,dd) = interp2(xx, yy, LUTs_eta2(:,:,dd), xxi, yyi);
    LUTs_theta_re(:,:,dd) = interp2(xx, yy, LUTs_theta(:,:,dd), xxi, yyi);
    LUTs_theta2_re(:,:,dd) = interp2(xx, yy, LUTs_theta2(:,:,dd), xxi, yyi);
    
    for D_ind = 1:length(Ds_re)
        T2_ind = 1;
        for T2 = T2s_re
            LUTs_index_re(D_ind, T2_ind, dd, 1) = Ds_re(D_ind)/x2_norm;
            LUTs_index_re(D_ind, T2_ind, dd, 2) = T2/x1_norm;
            LUTs_index_re(D_ind, T2_ind, dd, 3) = dphis(dd);
            T2_ind = T2_ind + 1;
        end
    end

end


LUTs.mag = LUTs_mag_re;
LUTs.magSm = LUTs_magSm_re;
LUTs.DESS = LUTs_DESS_re;
LUTs.DIFF = LUTs_DIFF_re;
LUTs.dS = LUTs_dS_re;
LUTs.eta = LUTs_eta_re;
LUTs.eta2 = LUTs_eta2_re;
LUTs.theta = LUTs_theta_re;
LUTs.theta2 = LUTs_theta2_re;
LUTs.index = LUTs_index_re;
LUTs.T2s = T2s_re;
LUTs.Ds = Ds_re;
LUTs.dphis = dphis;
LUTs.dict = LUTs_dict_re;
LUTs.mdl_theta = mdl_theta;
LUTs.mdl_theta2 = mdl_theta2;
LUTs.x1 = x1_norm;
LUTs.x2 = x2_norm;
