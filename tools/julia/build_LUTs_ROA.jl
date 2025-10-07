
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




using Interpolations, LinearAlgebra, Printf
using .Threads # Using Julia's built-in multi-threading

include("analytical_SPGR_W_diffusion.jl")
include("build_models.jl")

# Define a struct for sequence parameters to make passing them easier
struct Params
    dphi::Float64
    G1::Float64
    G2::Float64
    TR::Float64
    te1::Float64
    FA::Float64
    opuser8::Int
end

# A mutable struct to hold the (potentially large) Look-Up Tables
mutable struct LUTsData
    theta
    theta2
    T2s
    Ds
    FAs
    mdl_theta
    mdl_theta2
    x1 # Normalization factor for T2 (x1_norm)
    x2 # Normalization factor for D (x2_norm)
    
    # Other LUTs from original file (optional, can be added if needed)
    mag; magSm; DESS; DIFF; dS; eta; eta2; index; dict;
    
    # Constructor for empty initialization
    LUTsData() = new()
end

"""
    build_LUTs_ROA(params; T1=1000e-3, FAs_in=nothing)

Julia translation of `build_LUTs_ROA.m`. It constructs look-up tables (LUTs) for MR parameters and then builds spline models for fast mapping.
It uses multi-threading to accelerate the main loop.
"""
function build_LUTs_ROA(params::Params; T1=1000e-3, FAs_in=nothing)
    
    T2s = (20:20:300) .* 1e-3
    Ds = (250:250:3000) .* 1e-12
    T2s_re = (20:1:300) .* 1e-3
    Ds_re = (250:10:3000) .* 1e-12

    FAs = isnothing(FAs_in) ? [params.FA] : FAs_in
    
    normalize_sig = true
    x1_norm, x2_norm = if normalize_sig
        maximum(T2s), maximum(Ds)
    else
        1.0, 1.0
    end

    area = [params.G1, params.G2]
    G_ave = 100 .* area ./ (params.TR * 1e6)
    Gamma = 4285 * 2 * pi # [rad/Gauss]
    G = G_ave .* Gamma # [rad/m]

    TR = params.TR

    # Initialize result arrays
    nDs, nT2s, nFAs = length(Ds), length(T2s), length(FAs)
    LUTs_theta = zeros(nDs, nT2s, nFAs)
    LUTs_theta2 = zeros(nDs, nT2s, nFAs)
    # Note: Only theta and theta2 are required for the TV_mapping function.
    # Other LUTs can be added here if needed for other purposes.

    @printf("Building LUTs for T1 = %.2f ms using %d threads...\n", T1 * 1000, nthreads())

    # The MATLAB code uses parfor; here we use Julia's multi-threading.
    # This loop is thread-safe as each iteration writes to a unique `D_ind`.
    @threads for D_ind in 1:nDs
        for (fa_ind, FA) in enumerate(FAs)
            C0 = params.dphi * (pi/180) * 0.5
            C = [params.dphi, -1.0 * params.dphi] .* (pi/180) * 0.5
            alpha = FA * (pi/180)

            for (T2_ind, T2) in enumerate(T2s)
                y0, _, _, _ = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C0, Ds[D_ind], G[1], TR)
                y1, _, _, _ = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C[1], Ds[D_ind], G[2], TR)
                y1r, _, _, _ = analytical_SPGR_W_diffusion(TR, T1, T2, alpha, C[2], Ds[D_ind], G[2], TR)
    
                LUTs_theta[D_ind, T2_ind, fa_ind] = abs(angle(y0))
                LUTs_theta2[D_ind, T2_ind, fa_ind] = abs(angle(y1 * conj(y1r)) / 2 - pi/2)
            end
        end
    end
    
    println("Interpolating LUTs to finer grid...")
    nDs_re, nT2s_re = length(Ds_re), length(T2s_re)
    LUTs_theta_re = zeros(nDs_re, nT2s_re, nFAs)
    LUTs_theta2_re = zeros(nDs_re, nT2s_re, nFAs)

    knots = (Ds, T2s) # Define the grid for the original LUTs
    for dd in 1:nFAs
        # Create interpolants and evaluate on the finer grid (Ds_re, T2s_re)
        itp_theta = cubic_spline_interpolation(knots, LUTs_theta[:,:,dd], extrapolation_bc=Line())
        LUTs_theta_re[:,:,dd] = [itp_theta(d, t2) for d in Ds_re, t2 in T2s_re]
        
        itp_theta2 = cubic_spline_interpolation(knots, LUTs_theta2[:,:,dd], extrapolation_bc=Line())
        LUTs_theta2_re[:,:,dd] = [itp_theta2(d, t2) for d in Ds_re, t2 in T2s_re]
    end

    println("Building spline models...")
    mdl_theta, mdl_theta2 = if isnothing(FAs_in) && nFAs == 1
        data_theta = normalize_sig ? permutedims(LUTs_theta[:,:,1]) : permutedims(LUTs_theta[:,:,1])
        data_theta2 = normalize_sig ? permutedims(LUTs_theta2[:,:,1]) : permutedims(LUTs_theta2[:,:,1])
        
        knots_x1 = normalize_sig ? T2s ./ x1_norm : T2s
        knots_x2 = normalize_sig ? Ds ./ x2_norm : Ds

        # Data needs to be permuted because `build_models` expects dimensions (x1, x2)
        # and our LUTs have dimensions (Ds, T2s) i.e. (x2, x1)
        build_models(knots_x1, knots_x2, data_theta), build_models(knots_x1, knots_x2, data_theta2)
    else
        nothing, nothing # Models are not built if multiple FAs are processed
    end
    
    # Package results into the struct
    luts = LUTsData()
    luts.theta = LUTs_theta_re
    luts.theta2 = LUTs_theta2_re
    luts.T2s = T2s_re
    luts.Ds = Ds_re
    luts.FAs = FAs
    luts.mdl_theta = mdl_theta
    luts.mdl_theta2 = mdl_theta2
    luts.x1 = x1_norm
    luts.x2 = x2_norm
    
    return luts
end

