function optimize_lambda(ω_exp,Z_exp,τ)
"""
    Description
    -----------
    Main function for finding optimal λ hyperparameter for regularization.

    Parameters
    -----------
    ω_exp   - Input EIS frequency
    Z_exp   - Input EIS Impedance
    τ       - Desired relaxation times for computing DRT: generated in compute_drt
"""
    lambda_values = vcat(logrange(1e-06,1e01,8))
    lambda_crossval = Vector(undef,length(lambda_values))
    ω = 1 ./τ

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω)

    # Initialize the parameters
    n = length(ω) + 1
    
    #finding λ that minimizes error via crossvalidation
    for i in eachindex(lambda_values)
        λ = lambda_values[i]
        p0 = fill(0.05,n)

        function drt_fit_real(Z_real,Z_exp,p;λ_reg)
            Z_drt_real = drt_Zreal_regular(Z_real,real(Z_exp),λ_reg,p)
            return Z_drt_real
        end 
        function drt_fit_imag(Z_imag,Z_exp,p;λ_reg)
            Z_drt_imag = drt_Zimag_regular(Z_imag,imag(Z_exp),λ_reg,p)
            return Z_drt_imag
        end 

        fit_funct_real = (ω,p) -> drt_fit_real(Z_real,Z_exp,p;λ_reg = λ)
        fit_funct_imag = (ω,p) -> drt_fit_imag(Z_imag,Z_exp,p;λ_reg = λ)

        fit_real = curve_fit(fit_funct_real, ω_exp, real(Z_exp), p0;lower = zeros(n),autodiff=:forwarddiff)
        fit_imag = curve_fit(fit_funct_imag, ω_exp, imag(Z_exp), p0;lower = zeros(n),autodiff=:forwarddiff)

        p_real = fit_real.param
        p_imag = fit_imag.param

        crossval_real = norm(p_imag[1] .+ Z_real*p_imag[2:end] - real(Z_exp))^2
        crossval_imag = norm(Z_imag*p_real[2:end] - imag(Z_exp))^2
        lambda_crossval[i] = crossval_real+crossval_imag
    end


    i_min = argmin(lambda_crossval)[1]
    println("Regularization")
    println("--------------")
    println("λ = $(lambda_values[i_min])")

    return lambda_values[i_min]
end

function tune_τ(ω_exp,Z_exp;ppd=ppd,tol = 1e-03)
"""
    Description
    -----------
    Function for finding optimal τ range for computing DRT.
    Currently incomplete. Would be used in compute_drt

    Parameters
    -----------
    ω_exp   - Input EIS frequency
    Z_exp   - Input EIS Impedance
    ppd     - Points-per-decade in output τ range for computing DRT
    tol     - Tolerance for finding τ bounds where DRT impedance is sufficiently small
"""
    # fit = compute_drt(ω_exp,Z_exp)
    τ_init = logrange(0.1/maximum(ω_exp),10/minimum(ω_exp),floor(Int,log10(100*maximum(ω_exp)/minimum(ω_exp)))*ppd)
    ω_init = 1 ./τ_init
    n = length(ω_init) + 1
    dlnτ = log(τ_init[end]/τ_init[end-1])

    #cutting down data if it's too dense
    while length(ω_exp)>=length(ω_init)
        ω_exp = ω_exp[1:2:end]
        Z_exp = Z_exp[1:2:end]
    end

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω_init)

    function drt_fit(ω,p)
        Z_drt = drt_Z(Z_real,Z_imag,p)
        return vcat(real(Z_drt),imag(Z_drt))
    end 
    fit_funct = drt_fit

    p0 = abs.(rand(n))
    fit = curve_fit(fit_funct, ω_exp, vcat(real(Z_exp),imag(Z_exp)), p0;lower = zeros(n),autodiff=:forwarddiff)
    p = fit.param
    γ_init = p[2:end]/dlnτ
    γ_max = maximum(γ_init)
    peaks = findmaxima(γ_init)
    peaks = peakheights(peaks,min = tol*γ_max)
    τ_pks,γ_pks = τ_init[peaks.indices],peaks.heights
    Z_im = ω-> sum([γ_pks[i]*τ_pks[i]*10^ω/(1+τ_pks[i]^2*10^2ω) 
                    for i in eachindex(γ_pks)]) - tol*γ_max
    min,max = find_zeros(Z_im,-10,10)
    τ_min,τ_max = 10^-max ,10^-min ##the 2π was random and I don't think correct
    τ_tuned = logrange(τ_min,τ_max,floor(Int,log10(τ_max/τ_min))*ppd)
    return τ_tuned
end