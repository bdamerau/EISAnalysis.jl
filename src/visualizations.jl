function plot_Nyquist(a::Circuit...)
    plt = plot()
    for circuit in a
        scatter!(plt,circuit.Z)
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return plt
end
function plot_Nyquist!(plt,a::Circuit...)
    for circuit in a
        scatter!(plt,circuit.Z)
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    display(plt)
    return nothing
end