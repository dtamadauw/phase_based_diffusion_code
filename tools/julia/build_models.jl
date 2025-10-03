"""
Author: Daiki Tamada
Affiliation: Department of Radiology, University of Wisconsin-Madison
Date: 10/1/2025
Email: dtamada@wisc.edu
By downloading, installing, or otherwise accessing or using the Software , you (“Recipient”) agree to receive and use the above-
identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
obligations and restrictions you are not permitted to download, install,
execute, access, or use the SOFTWARE:

1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE (“Provider’s Rights”) to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation (“WARF”). 
2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider’s Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient’s responsibility. 
3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider’s proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN “AS IS” BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER’S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient’s use of the SOFTWARE.
7.	No use of names/trademarks.  Recipient shall not use Provider’s name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient’s possession, custody, or control.	

"""

using Interpolations

# Define a struct to hold the spline models, similar to MATLAB's struct
struct SplineModel
    sp::AbstractInterpolation
    dx1::AbstractInterpolation
    dx2::AbstractInterpolation
end

"""
    build_models(x1, x2, y_data)

Julia translation of `build_models.m`.
Creates a 2D cubic spline and its partial derivatives using the `Interpolations.jl` package.

# Arguments
- `x1`, `x2`: The grid vectors for the two dimensions.
- `y_data`: The 2D data array for which to create the spline. 
            Dimensions must be `(length(x1), length(x2))`.

# Returns
- A `SplineModel` struct containing the spline interpolant and its partial derivatives.
"""
function build_models(x1, x2, y_data)
    knots = (x1, x2)
    
    # Create the 2D cubic spline interpolation object.
    # We use `Line()` for extrapolation to mimic MATLAB's default behavior.
    itp = cubic_spline_interpolation(knots, y_data, extrapolation_bc=Line())
    
    # The modern approach is to compute the gradients at each grid point
    # and then create new interpolants for the derivative fields.
    
    # 1. Calculate the gradient (a tuple of derivatives) at each grid point.
    #    We use `Ref(itp)` to prevent broadcasting over the interpolant itself.
    #    `x2'` makes it a row vector, enabling broadcasting across the 2D grid.
    grads = Interpolations.gradient.(Ref(itp), x1, x2')
    
    # 2. Separate the gradient components into their own matrices.
    #    `getindex.(grads, 1)` pulls out the first element (dx1) from each tuple.
    dx1_data = getindex.(grads, 1)
    #    `getindex.(grads, 2)` pulls out the second element (dx2).
    dx2_data = getindex.(grads, 2)
    
    # 3. Create new spline interpolants for each partial derivative.
    itp_dx1 = cubic_spline_interpolation(knots, dx1_data, extrapolation_bc=Line())
    itp_dx2 = cubic_spline_interpolation(knots, dx2_data, extrapolation_bc=Line())
    
    # Return the models in a struct
    return SplineModel(itp, itp_dx1, itp_dx2)
end


