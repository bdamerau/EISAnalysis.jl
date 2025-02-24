"""
    evaluate_Z(f_r,f_c)

Computes elements for Z_matrices. The function based on impedance of an RC Circuit

# Arguments
- `f_r`: The Z_matrix row frequency
- `f_c`: The Z_matrix column frequency
"""
function evaluate_Z(f_r,f_c)
    var = f_r/f_c 
    return 1/(1+im*var)
end

"""
    regularizer(p,λ)

Computes the regularization of parameters. Used in regularization method of computing DRT.

See `drt_Z_regular`,`drt_Zimag_regular`, and `drt_Zreal_regular`

#Arguments
- `p`: The parameter list 
- `λ`: Regularization hyperparameter
"""
function regularizer(p,λ)

    return λ*norm(p)^2
end

"""
    build_Z_matrices(ω_in,ω_out)

Builds the matrices (real + imaginary) through which impedance is calculated.

See `compute_drt`, `optimize_lambda`, and `tune_τ`

#Attributes
- `ω_in`: The input EIS frequencies
- `ω_out`: he desired output frequencies for computing DRT
"""
function build_Z_matrices(ω_in,ω_out)
    Z_real,Z_imag = zeros(length(ω_in),length(ω_out)), zeros(length(ω_in),length(ω_out))

    for i in eachindex(ω_in), j in eachindex(ω_out)
        Z_element = evaluate_Z(ω_in[i],ω_out[j])
        Z_real[i,j] = real(Z_element)
        Z_imag[i,j] = imag(Z_element)
    end
    
    return Z_real,Z_imag
end

"""
    drt_Zreal_regular(X_r,Y_r,λ,p_reg)

See `drt_Z_regular`. Outputs the real component.

Used for crossvalidation in `optimize_lambda`.

# Attributes
- `X_r`: The real Z matrix for computing impedance
- `Y_r`: The real component of input impedance data
- `λ`: The regularization hyperparameter
- `p_reg`: The parameter list
"""
function drt_Zreal_regular(X_r,Y_r,λ,p_reg)
    R0,R_drt... = p_reg
    x = R0 .+ (X_r)*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A_r = Y_r + sqrt.(abs2.(x-Y_r) .+ Γ/N)
    return A_r
end

"""
    drt_Zimag_regular(X_i,Y_i,λ,p_reg)

See `drt_Z_regular`. Outputs the imaginary component.

The only difference with drt_Zreal_regular() is ignoring the series resistance. Used for crossvalidation in `optimize_lambda`.

# Attributes
- `X_i`: The imaginary Z matrix for computing impedance
- `Y_i`: The imaginary component of input impedance data
- `λ`: The regularization hyperparameter
- `p_reg`: The parameter list
"""
function drt_Zimag_regular(X_i,Y_i,λ,p_reg)
    R_drt = p_reg[2:end]
    x = X_i*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A_i = Y_i + sqrt.(abs2.(x-Y_i) .+ Γ/N)
    return A_i
end

"""
    drt_Z_regular(X_r,X_i,Y_r,Y_i,λ,p_reg)
    
See `drt_Z`. Incorporates regularization in least-squares minimization.

Since LsqFit minimizes Σ|Xᵢ-Yᵢ|², this function generates A such that
    |Aᵢ-Yᵢ|² = |Xᵢ-Yᵢ|² + (λ/N)|p|²
    Aᵢ = Yᵢ+√(|Xᵢ-Yᵢ|²+(λ/N)|p|²)
where N is the length of the vectors. Since the Real component carries 
the extra parameter (R0), two different functionsare defined for the real 
and imaginary components.Used for regularization method of computing DRT.

# Attributes
- `X_r`: The real Z matrix for computing impedance
- `X_i`: The imaginary Z matrix for computing impedance
- `Y_r`: The real component of input impedance data
- `Y_i`: The imaginary component of input impedance data
- `λ`: The regularization hyperparameter
- `p_reg`: The parameter list
"""
function drt_Z_regular(X_r,X_i,Y_r,Y_i,λ,p_reg)
    R0,R_drt... = p_reg
    x = R0 .+ (X_r+im*X_i)*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A = Y_r + im*Y_i + sqrt.(abs2.(x-Y_r-im*Y_i) .+ Γ/N)
    return A
end

"""
    drt_Z(X_r,X_i,p)
Uses Z matrices to compute impedance for least-squares fit.

# Attributes
- `X_r`: The real Z matrix for computing impedance
- `X_i`: The imaginary Z matrix for computing impedance
- `p`: The parameter list

"""
function drt_Z(X_r,X_i,p)

    R0,R_drt... = p
    return R0 .+ (X_r+im*X_i)*R_drt
end


"""
    compute_drt(ω_exp,Z_exp;ppd,showplot,rtol,regularization)

Main function for computing the DRT from input EIS data.

# Attributes
- `ω_exp`: Input EIS frequency
- `Z_exp`: Input EIS Impedance
- `ppd=7`: Points-per-decade in output τ range for computing DRT
- `showplot=false`: Whether or not to plot DRT results
- `rtol=1e-03`: The desired relative tolerance. Function outputs a warning if the fit result is above this
- `regularization=false`: Whether or no to employ the regularization method

# Examples
```julia
using EISAnalysis, Statistics
eval(initialize())
ω_exp, Z_exp = (r/q).ω, (r/q).Z #replace this with real data
fit = compute_drt(ω_exp,Z_exp;showplot = false,regularization = true)
```
"""
function compute_drt(ω_exp,Z_exp;ppd = 7,showplot = true,rtol = 1e-03,regularization = false)
    τ = logrange(0.1/maximum(ω_exp),10/minimum(ω_exp),floor(Int,log10(100*maximum(ω_exp)/minimum(ω_exp)))*ppd)
    # τ = tune_τ(ω_exp,Z_exp;ppd=ppd)
    ω = 1 ./τ
    n = length(ω) + 1
    dlnτ = log(τ[end]/τ[end-1])

    #cutting down data if it's too dense
    while length(ω_exp)>=length(ω)
        ω_exp = ω_exp[1:2:end]
        Z_exp = Z_exp[1:2:end]
    end

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω)

    if regularization
        λ = optimize_lambda(ω_exp,Z_exp,τ)
        function drt_fit_regular(ω,p)
            Z_drt = drt_Z_regular(Z_real,Z_imag,real(Z_exp),imag(Z_exp),λ,p)
            return vcat(real(Z_drt),imag(Z_drt))
        end 
        fit_funct = drt_fit_regular
        p0 = fill(0.5,n)
    else
        function drt_fit(ω,p)
            Z_drt = drt_Z(Z_real,Z_imag,p)
            return vcat(real(Z_drt),imag(Z_drt))
        end 
        fit_funct = drt_fit
        p0 = abs.(rand(n))
    end 

        fit = curve_fit(fit_funct, ω_exp, vcat(real(Z_exp),imag(Z_exp)), p0;lower = zeros(n),autodiff=:forwarddiff)
        p = fit.param
        Z_fit = drt_Z(Z_real,Z_imag,p)

        loss = mean(abs2.((Z_fit.-Z_exp)./Z_exp))
        println("rerror = $loss")
        if loss > rtol
            println("WARNING: error is above specified tolerance")
        end

    γ_fit = p[2:end]/dlnτ
    if showplot
        Z_matrices = build_Z_matrices(ω,ω)
        Z_expanded = drt_Z(Z_matrices[1],Z_matrices[2],p)
        plt = plot_drt(Z_exp,Z_fit,Z_expanded,τ,γ_fit)
        display(plt)
    end

    return Dict([
        "Z"=>Z_fit
        "R0"=>p[1]
        "drt"=>γ_fit
        "τ"=>τ
    ])
end