
#Author: Daiki Tamada
#Affiliation: Department of Radiology, University of Wisconsin-Madison
#Date: 10/1/2025
#Email: dtamada@wisc.edu
#By downloading, installing, or otherwise accessing or using the Software , you (“Recipient”) agree to receive and use the above-
#identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
#obligations and restrictions you are not permitted to download, install,
#execute, access, or use the SOFTWARE:

#1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE (“Provider’s Rights”) to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation (“WARF”). 
#2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider’s Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient’s responsibility. 
#3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider’s proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
#4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
#5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN “AS IS” BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
#6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER’S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient’s use of the SOFTWARE.
#7.	No use of names/trademarks.  Recipient shall not use Provider’s name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
#8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient’s possession, custody, or control.	




using Printf, Statistics
using Plots, Measures # Added for plotting
using Random
using MAT

# Include the translated function modules
include("build_LUTs_ROA.jl")
include("TV_mapping_no_penalty_ori.jl")
include("TV_mapping_pixel_by_pixel.jl") # <-- ADD THIS INCLUDE

"""
    plot_luts(LUTs)

Generates and displays heatmap plots for LUTs.theta and LUTs.theta2.
"""
function plot_luts(LUTs)
    println("\nGenerating LUT plots...")
    
    # Convert units for better plot readability
    # T2s from s to ms
    T2s_ms = LUTs.T2s .* 1000
    # Ds from m²/s to μm²/ms
    Ds_um2_ms = LUTs.Ds .* 1e9

    # CORRECTED: Select the first slice from the 3rd (FA) dimension to get a 2D matrix
    theta_slice = LUTs.theta[:, :, 1]
    theta2_slice = LUTs.theta2[:, :, 1]

    # Plot for theta (phi_H)
    p1 = heatmap(
        T2s_ms,
        Ds_um2_ms,
        theta_slice', # Transpose the 2D slice
        xlabel="T2 (ms)",
        ylabel="D (μm²/ms)",
        title="LUT for θ (High Gradient Moment Phase)",
        colorbar_title="\nPhase (radians)",
        right_margin=10mm
    )

    # Plot for theta2 (phi_L)
    p2 = heatmap(
        T2s_ms,
        Ds_um2_ms,
        theta2_slice', # Transpose the 2D slice
        xlabel="T2 (ms)",
        ylabel="D (μm²/ms)",
        title="LUT for θ2 (Low Gradient Moment Phase)",
        colorbar_title="\nPhase (radians)",
        right_margin=10mm
    )

    # Combine plots and display
    combined_plot = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    display(combined_plot)
    println("Plots displayed. Close the plot window to continue the script.")
end


"""
Main function to demonstrate the usage of the translated Julia code.
It defines parameters, builds LUTs, creates synthetic test data,
runs the mapping algorithm, and prints the results.
"""
function main(Total_Correction_Factor)
    println("Starting Julia translation of TV mapping...")


    TR = 10.9e-3;
    numRO = 256;
    Gamma = 4285 * 2* pi; #[rad/Gauss/s]
    FOV = 0.256;
    BW = 4*pi*2*(125)*1e+3/numRO;#[Rad/px]
    Resolution = FOV/numRO;#[m/px]
    Tsamp = 1/BW;
    G = 2*pi/(Gamma*Resolution)*(1/TR);
    G_2pi = G*Gamma;
    SNR = 2:0.5:50;
    noise = 1 ./ SNR;
    nex = 100000;
    FA = 20;

    #EPI parameters
    S0_EPI = 1.0;
    TR_EPI = 5000e-3;
    TE_EPI = 82e-3;

    
    @printf("Total_Correction_Factor: %.4f\n", Total_Correction_Factor)

    @printf("Gradiet Moment (2*pi): %.4f (Gauss/m)\n", G)
    Random.seed!(321)


    # 1. Define sequence parameters (as a struct)
    params = Params(
        2.0,      # dphi
        2.3337e+03,      # G1
        1.6336e+04, # G2
        TR,    # TR
        0.0,   # te1
        FA,     # FA
        0         # opuser8
    )

    dphi = 2;
    FA = 20;
    C = (pi/180)*dphi/2;
    alpha = FA*(pi/180);
    T1 = 1000*1e-3;
    #https://www.sciencedirect.com/science/article/pii/S1076633218301880
    T2 = [100, 120, 200]*1e-3;
    Ds = [800, 1300, 1800]*1e-12;

    G1 = G_2pi
    G2 = 7*G_2pi


    #1a. Generatae clean data
    y1p = zeros(Complex{Float64}, nex, length(noise), length(Ds));
    y2p = zeros(Complex{Float64}, nex, length(noise), length(Ds));
    y1n = zeros(Complex{Float64}, nex, length(noise), length(Ds));
    y2n = zeros(Complex{Float64}, nex, length(noise), length(Ds));

    ys = zeros(Complex{Float64}, 1,4);
    noise_dit = zeros(Complex{Float64},nex,length(noise),4);


    for kk=1:length(Ds)

        S_EPI_b_low = S0_EPI * exp(-TE_EPI/T2[kk]) * (1.0-exp(-TR_EPI/1000e-3));
        noise_level = Total_Correction_Factor*S_EPI_b_low;

        ys[1],~,~,~ = analytical_SPGR_W_diffusion(TR, T1, T2[kk], alpha, C, Ds[kk], G1, TR);
        ys[2],~,~,~ = analytical_SPGR_W_diffusion(TR, T1, T2[kk], alpha, C, Ds[kk], G2, TR);
        ys[3],~,~,~ = analytical_SPGR_W_diffusion(TR, T1, T2[kk], alpha, -C, Ds[kk], G1, TR);
        ys[4],~,~,~ = analytical_SPGR_W_diffusion(TR, T1, T2[kk], alpha, -C, Ds[kk], G2, TR);

        for ii=1:length(noise)
            n_s = zeros(Complex{Float64},nex,4);
            for jj=1:10
                Random.seed!(42+2*jj+2*kk);#Reduce seed-dependency
                sf = Integer(nex/10);
                ind_seed = ((jj-1)*sf + 1):(jj*sf);
                n_s[ind_seed,:] = noise_level*complex.(noise[ii].*randn(sf,4),noise[ii].*randn(sf,4))/sqrt(2);
            end
            noise_dit[:,ii,:] = n_s;
            y1p[:,ii,kk] = ys[1].+n_s[:,1]; y2p[:,ii,kk] = ys[2].+n_s[:,2];
            y1n[:,ii,kk] = ys[3].+n_s[:,3]; y2n[:,ii,kk] = ys[4].+n_s[:,4];
        end
    end
 
    F0 = angle.(y1n .* conj(y1p)) ./ 2 .+ pi/2;
    F1 = angle.(y2n .* conj(y2p)) ./ 2 .+ pi/2;
    PH0 = angle.(y1n.*conj(y1p))./2 .+pi/2;
    PH1 = angle.(y2n.*conj(y2p))./2 .+pi/2;


    # 2. Build the Look-Up Tables (LUTs) and spline models
    println("\nBuilding LUTs... (This may take a while)")
    @time LUTs = build_LUTs_ROA(params, T1=1.0)
    println("LUTs and spline models built successfully.")



    T2_PBD = reshape(PH0,size(PH0,1), size(PH0,2), size(PH0,3));
    D_PBD = reshape(PH1,size(PH1,1), size(PH1,2), size(PH1,3));
    for kk = 1:size(PH0,3)
        @time T2_PBD[:,:,kk], D_PBD[:,:,kk] = TV_mapping_pixel_by_pixel(PH0[:,:,kk], PH1[:,:,kk], LUTs);
    end

    T2_PBD = T2_PBD .* LUTs.x1;
    D_PBD = D_PBD .* LUTs.x2;

    data_to_save = Dict(
        "PBD_jul" => D_PBD,
        "PBT2_jul" => T2_PBD
    )


    filename = s = "PBD_Julia_$Total_Correction_Factor.mat"
    matwrite(filename, data_to_save)

end

println("Arguments received:")
for arg in ARGS
    println("- ", arg)
end

noise_factor = parse(Float64, ARGS[1])

main(noise_factor)

