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


function [Msig,Rfph]=SPGRE_Bloch_diffusion(T1,T2,TR,flip,inc,Nex,Nf,D,TG,G)
%[Msig,Rfph]=SPGRE_Bloch_diffusion(T1,T2,TR,flip,inc,Nex,Nf,D,TG,G)
%T1: T1 [s]
%T2: T2 [s]
%TR: Repetition TIme [s]
%flip: Flip angle [degree]
%inc: RF Phase increment [degree]
%Nex: The Number of TRs
%Nf: The Number of isochromats, larger number is better
%D: Diffusion coefficient [m^2/s]
%TG: Duration time for gradient [s]
%G: Gradient strength [G/m]
%
%This code is implemented base on the Yarnykh's paper
%Yarnykh, Vasily L. "Optimal radiofrequency and gradient spoiling for improved accuracy of T1 and B1 measurements using fast steady‐state techniques." Magnetic Resonance in Medicine: 63.6 (2010): 1610-1626.
    
M0=zeros(3,Nf);
M=zeros(3,Nf);
Msig = zeros(1,Nex);

%Relax matrix for T1 and T2
E1 = exp(-TR/T1);	
E2 = exp(-TR/T2);
Bfp = [0 0 1-E1]';

gamma = 2*pi*42.58e+2;%Gyromagnetic ratio [rad/G]
Vox = 0.001;% We assumed voxel size of 1 mm

%Equilibrium magnetization
%We assume isochromats are alligned along gradient direction from 1 to Nf
M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);

%Array for RF phase increment
Rfph(1) = 0;
Rfinc = inc;

%Calculate gradient area
Ag = G*TG*gamma;


%Define gaussian distribution to simulate attenuation due to random walk
P = @(phi) Vox*(1/Nf)*Ag*(1/(2*Ag*sqrt(pi*D*TR)))*exp( -1*((phi).*(phi))/(4*D*TR*Ag*Ag) );

%Relax operator with diffusion
Ed = exp(-1*Ag*Ag*TG*D/12);
E = @(phi) [Ed*E2*cos(phi/2) Ed*E2*sin(phi/2) 0;-Ed*E2*sin(phi/2) Ed*E2*cos(phi/2) 0;0 0 E1];

%Define range to consider diffusion effect.
%In this case, we used 5 sigma of the distribution
Nmax = ceil(5*Ag*sqrt(2*D*TR)/(Vox*(1/Nf)*Ag)); 

for n=1:Nex

    %Excitation with RF phase increment
    R = throt(flip*pi/180,Rfph(n)*pi/180);
    M0 = R*M;

    %Acquisition
    Msig(n) = mean( complex(M0(1,:),M0(2,:)), 2) * exp(-i*Rfph(n)*pi/180);


    %Diffusion term
    for ii=1:Nf
	%Calculate phase of an isochromat at ii th position
        phi_l = Vox*((ii-floor(Nf/2)-1)/Nf)*Ag;
        M(:,ii) = 0;
        q_range = (ii-Nmax):(ii+Nmax);
        q_range(q_range<1)=q_range(q_range<1)+Nf;
        q_range(q_range>Nf) = q_range(q_range>Nf)-Nf;
        
        Parr=[];phiarr=[];
        for jj=1:(2*Nmax+1)
            q = q_range(jj);
            phi_q = Vox*(( (ii+jj-Nmax-1) -floor(Nf/2)-1)/Nf)*Ag;
            M(:,ii) = M(:,ii) + P(Vox*Ag*(jj-Nmax-1)/Nf)*E(phi_l+phi_q)*M0(:,q);
                        
        end
        M(:,ii) = M(:,ii) + Bfp;
    end
    
    Rfph(n+1) = Rfph(n)+Rfinc;
    Rfinc = (n+1)*inc;
end;


%%
function Rth=throt(phi,theta)

Rz = zrot(-theta);
Rx = xrot(phi);
Rth = inv(Rz)*Rx*Rz;

function Rz=zrot(phi)

Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];

function Rx=xrot(phi)

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
