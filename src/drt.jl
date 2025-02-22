#calculate impedance for a given frequency and timescale (not including resistance)
function evaluate_Z(f_r,f_c)
    var = f_r/f_c 
    return 1/(1+im*var)
end

#regularization function
function regularizer(p,λ)
    return λ*norm(p)^2
end
#building the matrix X above
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
    Since LsqFit minimizes Σ|Xᵢ-Yᵢ|², this function generates A such that
                |Aᵢ-Yᵢ|² = |Xᵢ-Yᵢ|² + (λ/N)|p|²,    where N is the length of the vectors
                      Aᵢ = Yᵢ+√(|Xᵢ-Yᵢ|²+(λ/N)|p|²)
    Since the Real component carries the extra parameter (R0), two different functions
    are defined for the real and imaginary components
"""
function drt_Zreal_regular(X_r,Y_r,λ,p_reg)
    R0,R_drt... = p_reg
    x = R0 .+ (X_r)*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A_r = Y_r + sqrt.(abs2.(x-Y_r) .+ Γ/N)
    return A_r
end
function drt_Zimag_regular(X_i,Y_i,λ,p_reg)
    R_drt = p_reg[2:end]
    x = X_i*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A_i = Y_i + sqrt.(abs2.(x-Y_i) .+ Γ/N)
    return A_i
end
function drt_Z_regular(X_r,X_i,Y_r,Y_i,λ,p_reg)
    R0,R_drt... = p_reg
    x = R0 .+ (X_r+im*X_i)*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A = Y_r + im*Y_i + sqrt.(abs2.(x-Y_r-im*Y_i) .+ Γ/N)
    return A
end

#calculate the impedance, given the matrices and the parameters
function drt_Z(X_r,X_i,p)
    R0,R_drt... = p
    return R0 .+ (X_r+im*X_i)*R_drt
end


#computing DRT
function compute_drt(ω_exp,Z_exp;ppd = 7,showplot = true,rtol = 1e-03,regularization = false)

    # τ = logrange(0.1/maximum(ω_exp),10/minimum(ω_exp),floor(Int,log10(100*maximum(ω_exp)/minimum(ω_exp)))*ppd)
    τ = tune_τ(ω_exp,Z_exp;ppd=ppd)
    ω = 1 ./τ
    n = length(ω) + 1
    dlnτ = log(τ[end]/τ[end-1])

    #cutting down data if it's too dense
    while length(ω_exp)>=length(ω)
        ω_exp = ω_exp[1:2:end]
        Z_exp = Z_exp[1:2:end]
    end

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω)

    ####################
    if regularization
        λ = optimize_lambda(ω_exp,Z_exp,τ)
        function drt_fit_regular(ω,p)
            Z_drt = drt_Z_regular(Z_real,Z_imag,real(Z_exp),imag(Z_exp),λ,p)
            return vcat(real(Z_drt),imag(Z_drt))
        end 
        fit_funct = drt_fit_regular
        p0 = fill(0.05,n)
    else
        function drt_fit(ω,p)
            Z_drt = drt_Z(Z_real,Z_imag,p)
            return vcat(real(Z_drt),imag(Z_drt))
        end 
        fit_funct = drt_fit
        p0 = abs.(rand(n))
    end


        # Initialize the parameters
        
        
        # fit_funct = drt_fit       

        fit = curve_fit(fit_funct, ω_exp, vcat(real(Z_exp),imag(Z_exp)), p0;lower = zeros(n),autodiff=:forwarddiff)
        p = fit.param
        Z_fit = drt_Z(Z_real,Z_imag,p)

        loss = mean(abs2.((Z_fit.-Z_exp)./Z_exp))
        println("rtol = $loss")
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



function optimize_lambda(ω_exp,Z_exp,τ)

    # lambda_values = vcat(0,logrange(1e-06,1e01,8))
    lambda_values = vcat(logrange(1e-06,1e01,8))
    lambda_crossval = Vector(undef,length(lambda_values))
    ω = 1 ./τ

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω)

    # function drt_fit(ω,p)
    #     Z_drt = drt_Z(Z_real,Z_imag,p)
    #     return vcat(real(Z_drt),imag(Z_drt))
    # end 
    function drt_fit_real(ω,p;λ_reg)
        Z_drt_real = drt_Zreal_regular(Z_real,real(Z_exp),λ_reg,p)
        return Z_drt_real
    end 
    function drt_fit_imag(ω,p;λ_reg)
        Z_drt_imag = drt_Zimag_regular(Z_imag,imag(Z_exp),λ_reg,p)
        return Z_drt_imag
    end 


    # Initialize the parameters
    n = length(ω) + 1
    
    #finding λ that minimizes error
    errors,ps = [],[]
    for i in eachindex(lambda_values)
        # if λ==0 #potentially delete
        #     p0 = abs.(rand(n))
        #     fit_funct = drt_fit
        # else
        λ = lambda_values[i]
        p0 = fill(0.05,n)
        #fit_funct = drt_fitregular(ω,p;λ_reg = λ)
        # end
        fit_funct_real(ω,p) = drt_fit_real(ω,p;λ_reg = λ)
        fit_funct_imag(ω,p) = drt_fit_imag(ω,p;λ_reg = λ)

        fit_real = curve_fit(fit_funct_real, ω_exp, real(Z_exp), p0;lower = zeros(n),autodiff=:forwarddiff)
        fit_imag = curve_fit(fit_funct_imag, ω_exp, imag(Z_exp), p0;lower = zeros(n),autodiff=:forwarddiff)

        p_real = fit_real.param
        p_imag = fit_imag.param

        # lambda_discrepencies[i] = norm(p_real-p_imag)
        crossval_real = norm(p_imag[1] .+ Z_real*p_imag[2:end] - real(Z_exp))^2
        crossval_imag = norm(Z_imag*p_real[2:end] - imag(Z_exp))^2
        lambda_crossval[i] = crossval_real+crossval_imag
        # p = fit.param
        # push!(ps,p)
        # Z_fit = drt_Z(Z_real,Z_imag,p)
        # loss = mean(abs2.((Z_fit.-Z_exp)./Z_exp))
        # push!(errors,loss)
    end


    i_min = argmin(lambda_crossval)[1]
    # println(lambda_crossval)#testing
    println("Regularization")
    println("--------------")
    println("λ = $(lambda_values[i_min])")

    return lambda_values[i_min]
    # return lambda_values[2]##testing

    # λ_hat = lambda_values[argmin(errors)[1]]
    # p_hat = ps[argmin(errors)[1]]
    # if λ_hat ==0
    #     println("Regularizaiton Not Used")
    # else
    #     println("Regularization")
    #     println("--------------")
    #     println("λ = $λ_hat")
    # end
    # return p_hat,minimum(errors)
end

function plot_drt(Z_exp,Z_fit,Z_expanded,τ,γ)
    fitplt = scatter(Z_exp,label = "data")
    scatter!(fitplt,Z_fit,markersize = 3,label = "fit")

    drtplt = plot(τ,γ,lw=3)

    #calcualate expanded Z 
    expandedfitplt=scatter(Z_expanded, color=palette(:default)[2])
    R_drt = γ*log(τ[end]/τ[end-1])
    rcs =  [ @. real(Z_expanded[i]) - 0.5R_drt[i]*(cos(0:π/30:π)+ im*sin(0:π/30:π)) for i in eachindex(τ)]
    for rc in rcs
        plot!(expandedfitplt,rc,c=:purple,ls=:dash,lw=2)
    end

    #formatting the figures
    plot!(fitplt,yflip=true,aspect_ratio=:equal,legend = :topleft,ylabel = "-Im(Z) / Ω",xlabel = "Re(Z) / Ω",title = "Fit")
    plot!(drtplt,ylabel = "γ / Ω",xlabel = "τ / s",xaxis=:log,title = "DRT",legend = false,lw=3)
    plot!(expandedfitplt,yflip=true,aspect_ratio=:equal,legend = false,ylabel = "-Im(Z) / Ω",xlabel = "Re(Z) / Ω",title = "Expanded Fit")
    l = @layout [
        a b; c
    ]
    fullplt = plot(fitplt,drtplt,expandedfitplt,layout = l)
    return fullplt
end

function tune_τ(ω_exp,Z_exp;ppd=ppd,tol = 1e-03)
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
    τ_pks,γ_pks = peaks.indices,peaks.heights
    Z_im = ω-> sum([γ_pks[i]*τ_pks[i]*10^ω/(1+τ_pks[i]^2*10^2ω) 
                    for i in eachindex(γ_pks)]) - tol*γ_max
    min,max = find_zeros(Z_im,-10,10)
    τ_min,τ_max = 10^-max,10^-min
    τ_tuned = logrange(τ_min,τ_max,floor(Int,log10(τ_max/τ_min))*ppd)
    return τ_tuned
end