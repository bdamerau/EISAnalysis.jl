using Statistics,Plots,LsqFit,LinearAlgebra

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

#trying to add regularization
function drt_Zregular(X_r,X_i,Y_r,Y_i,λ,p_reg)
    R0,R_drt... = p_reg
    x = R0 .+ (X_r+im*X_i)*R_drt
    N = length(x)
    Γ = regularizer(p_reg,λ)
    A = sqrt.(abs2.(x-Y_r-im*Y_i) .+ Γ/N)
    A += Y_r + im*Y_i
    return A
end

#calculate the impedance, given the matrices and the parameters
function drt_Z(X_r,X_i,p)
    R0,R_drt... = p
    return R0 .+ (X_r+im*X_i)*R_drt
end


#computing DRT
function compute_drt(ω_exp,Z_exp;ppd = 7,showplot = true,rtol = 1e-03,regularization = false)

    τ = logrange(0.1/maximum(ω_exp),10/minimum(ω_exp),floor(Int,log10(100*maximum(ω_exp)/minimum(ω_exp)))*ppd)
    ω = 1 ./τ

    #cutting down data if it's too dense
    while length(ω_exp)>=length(ω)
        ω_exp = ω_exp[1:2:end]
        Z_exp = Z_exp[1:2:end]
    end

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω)

    ####################
    if regularization
        p,loss = optimize_lambda(ω_exp,Z_exp,τ)
        Z_fit = drt_Z(Z_real,Z_imag,p)
        
        println("rtol = $loss")
        if loss > rtol
            println("WARNING: error is above specified tolerance")
        end
    else
        function drt_fit(ω,p)
            Z_drt = drt_Z(Z_real,Z_imag,p)
            return vcat(real(Z_drt),imag(Z_drt))
        end 
        function drt_fitregular(ω,p)
            Z_drt = drt_Zregular(Z_real,Z_imag,real(Z_exp),imag(Z_exp),λ,p)
            return vcat(real(Z_drt),imag(Z_drt))
        end 


        # Initialize the parameters
        n = length(ω) + 1
        p0 = abs.(rand(n))
        fit_funct = drt_fit

        fit = curve_fit(fit_funct, ω_exp, vcat(real(Z_exp),imag(Z_exp)), p0;lower = zeros(n),autodiff=:forwarddiff)
        p = fit.param
        Z_fit = drt_Z(Z_real,Z_imag,p)

        loss = mean(abs2.((Z_fit.-Z_exp)./Z_exp))
        println("rtol = $loss")
        if loss > rtol
            println("WARNING: error is above specified tolerance")
        end
    end


    if showplot
        Z_matrices = build_Z_matrices(ω,ω)
        Z_expanded = drt_Z(Z_matrices[1],Z_matrices[2],p)
        plt = plot_drt(Z_exp,Z_fit,Z_expanded,τ,p)
        display(plt)
    end

    return Dict([
        "Z"=>Z_fit
        "R0"=>p[1]
        "drt"=>p[2:end]
        "τ"=>τ
    ])
end



function optimize_lambda(ω_exp,Z_exp,τ)

    lambda_values = vcat(0,logrange(1e-06,1e01,8))
    ω = 1 ./τ

    Z_real,Z_imag = build_Z_matrices(ω_exp,ω)

    function drt_fit(ω,p)
        Z_drt = drt_Z(Z_real,Z_imag,p)
        return vcat(real(Z_drt),imag(Z_drt))
    end 
    function drt_fitregular(ω,p;λ_reg)
        Z_drt = drt_Zregular(Z_real,Z_imag,real(Z_exp),imag(Z_exp),λ_reg,p)
        return vcat(real(Z_drt),imag(Z_drt))
    end 


    # Initialize the parameters
    n = length(ω) + 1
    
    #finding λ that minimizes error
    errors,ps = [],[]
    for λ in lambda_values
        if λ==0
            p0 = abs.(rand(n))
            fit_funct = drt_fit
        else
            p0 = fill(0.05,n)
            fit_funct(ω,p) = drt_fitregular(ω,p;λ_reg = λ)
        end

        fit = curve_fit(fit_funct, ω_exp, vcat(real(Z_exp),imag(Z_exp)), p0;lower = zeros(n),autodiff=:forwarddiff)

        p = fit.param
        push!(ps,p)
        Z_fit = drt_Z(Z_real,Z_imag,p)
        loss = mean(abs2.((Z_fit.-Z_exp)./Z_exp))
        push!(errors,loss)
    end

    λ_hat = lambda_values[argmin(errors)[1]]
    p_hat = ps[argmin(errors)[1]]
    if λ_hat ==0
        println("Regularizaiton Not Used")
    else
        println("Regularization")
        println("--------------")
        println("λ = $λ_hat")
    end
    return p_hat,minimum(errors)
end

function plot_drt(Z_exp,Z_fit,Z_expanded,τ,p)
    R0, R_drt... = p
    fitplt = scatter(Z_exp,label = "data")
    scatter!(fitplt,Z_fit,markersize = 3,label = "fit")
    plot!(fitplt,yflip=true,aspect_ratio=:equal,legend = :topleft,ylabel = "-Im(Z) / Ω",xlabel = "Re(Z) / Ω",title = "Fit")

    drtplt = plot(τ,R_drt,lw=3)
    plot!(drtplt,ylabel = "R / Ω",xlabel = "τ / s",xaxis=:log,title = "DRT",legend = false,lw=3)

    #calcualate expanded Z 
    expandedfitplt=scatter(Z_expanded, color=palette(:default)[2])
    rcs =  [ @. real(Z_expanded[i]) - 0.5R_drt[i]*(cos(0:π/30:π)+ im*sin(0:π/30:π)) for i in eachindex(τ)]
    for rc in rcs
        plot!(expandedfitplt,rc,c=:purple,ls=:dash,lw=2)
    end
    plot!(expandedfitplt,yflip=true,aspect_ratio=:equal,legend = false,ylabel = "-Im(Z) / Ω",xlabel = "Re(Z) / Ω",title = "Expanded Fit")

    l = @layout [
        a b; c
    ]
    fullplt = plot(fitplt,drtplt,expandedfitplt,layout = l)
    return fullplt
end