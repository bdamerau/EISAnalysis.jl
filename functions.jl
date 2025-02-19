using Statistics,Plots,LsqFit
include("library.jl")

#private functions
function build_fit(model)
    function circuit_fit(omega ,p)
        circuit = model(p)
        return vcat(real(circuit.Z),imag(circuit.Z))
    end
    return circuit_fit
end

#public functions

function plot_Nyquist(a::CircuitElement...)
    plt = plot()
    for circuit in a
        scatter!(plt,circuit.Z)
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return plt
end
function plot_Nyquist!(plt,a::CircuitElement...)
    for circuit in a
        scatter!(plt,circuit.Z)
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    display(plt)
    return nothing
end

function fitted_circuit(model, ω_data, Z_data, p0,bounds=[])
    fit_function = build_fit(model)
    if isempty(bounds)
        fit = curve_fit(fit_function, ω_data, vcat(real(Z_data),imag(Z_data)), p0;
        lower = zeros(length(p0)))
    else
        lower_b::Vector{Float64} = [b[1] for b in bounds]
        upper_b::Vector{Float64} = [b[2] for b in bounds]
        fit = curve_fit(fit_function, ω_data, vcat(real(Z_data),imag(Z_data)), p0;
        lower = lower_b,upper=upper_b)
    end
    circuit = model(fit.param)
    Z = circuit.Z
    plt = scatter(Z_data,label = "data")
    scatter!(plt,Z,markersize=3,label = "fit")
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    display(plt)
    println("Parameters")
    println("__________")
    [println(param) for param in fit.param]
    return Dict([
        "Z"=>Z
        "parameters"=>fit.param
    ])
end