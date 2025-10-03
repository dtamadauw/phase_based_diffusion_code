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

