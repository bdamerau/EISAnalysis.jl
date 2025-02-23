function plot_Nyquist(a::Circuit...)
"""
    Description
    -----------
    Parameters
    -----------
"""
    plt = plot()
    for circuit in a
        scatter!(plt,circuit.Z)
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return plt
end
function plot_Nyquist!(plt,a::Circuit...)
"""
    Description
    -----------
    Parameters
    -----------
"""
    for circuit in a
        scatter!(plt,circuit.Z)
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    display(plt)
    return nothing
end

function plot_drt(Z_exp,Z_fit,Z_expanded,τ,γ)
"""
    Description
    -----------
    Parameters
    -----------
"""
    fitplt = scatter(Z_exp,label = "data")
    scatter!(fitplt,Z_fit,markersize = 3,label = "fit")

    γ_pks = findmaxima(γ) |> peakproms!(min = maximum(γ)/20) |> peakwidths!()
    drtplt = plotpeaks(τ, γ; peaks=γ_pks.indices, prominences=true, widths=true,lw=2,ms = 2.5)

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

function print_circuit(circuit)
"""
    Description
    -----------
    Prints the elements of a circuit along with its parameters

    Parameters
    -----------
    circuit::Circuit    - The circuit being printed
"""
    for element in circuit.elements
        println(element)
    end
    return nothing
end