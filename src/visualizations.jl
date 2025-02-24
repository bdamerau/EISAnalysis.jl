"""
    plot_Nyquist(a;label)

Creates a Nyquist plot

# Arguments
- `a::Circuit`: The circuits to add to the Nyquist plot
- `label`: kwarg label for plot

# Examples
```julia
julia> using EISAnalysis
julia> eval(initialize());
julia> circuit1 = r-r/c;
julia> circuit2 = r-r/q;
julia> plot_Nyquist(circuit1,circuit2;label= ["R-C Circuit","R-CPE Circuit"]);
```
"""
function plot_Nyquist(a::Circuit...;label::Vector{String} = fill("",length(a)))
    if length(label)!=length(a)
        println("Warning: label must be the same length as data arguments")
        label = fill("",length(a))
    end

    plt = plot()
    for i in eachindex(a)
        scatter!(plt,a[i].Z;label = label[i])
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return plt
end
function plot_Nyquist(a::Circuit;label::String = "")
    plt = scatter(a.Z;label = label)
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return plt
end

"""
    plot_Nyquist!(plt,a;label)

Adds circuits to pre-existing Nyquist plot

# Arguments
- `plt`: The input Nyquist plot to manipulate
- `a::Circuit`: The circuits to add to the Nyquist plot
- `label`: kwarg label for plot

# Examples
```julia
julia> using EISAnalysis
julia> eval(initialize());
julia> circuit1 = r-r/c;
julia> circuit2 = r-r/q;
julia> plot_Nyquist(circuit1,circuit2;label= ["R-C Circuit","R-CPE Circuit"]);
```
"""
function plot_Nyquist!(plt,a::Circuit...;label::Vector{String} = fill("",length(a)))
    if length(label)!=length(a)
        println("Warning: label must be the same length as data arguments")
        label = fill("",length(a))
    end
    for i in eachindex(a)
        scatter!(plt,a[i].Z,label = label[i])
    end
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return nothing
end
function plot_Nyquist!(plt,a::Circuit;label::String = "")
    scatter!(plt,a.Z,label = label)
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    return nothing
end


function plot_drt(Z_exp,Z_fit,Z_expanded,τ,γ)
"""
Called in `compute_drt`

Creates a combined plot of 
1. A fit of the DRT results to input EIS data
2. The DRT plotted with peaks and peakwidths
3. An expanded Nyquist plot of DRT data to aid in DRT interpretation

# Arguments
- `Z_exp`: The input impedance data being fitted
- `Z_fit`: The DRT fit result, matching the frequencies of Z_exp
- `Z_expanded`: The expanded DRT fit using the full range of τ
- `τ`: The range of timescales over which the DRT is calculated
- `γ`: DRT Results

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

"""
    print_circuit(circuit)
    
Prints the elements of a circuit along with its parameters

# Arguments
- `circuit::Circuit`: The circuit being printed

# Examples
```julia
using EISAnalysis
eval(initialize())
randles_circuit = 0.23r-r/0.2q
print_circuit(randles_circuit)

# output 

0.23r
1.0r
0.2 * q ^ 0.8
```
"""
function print_circuit(circuit)
    for element in circuit.elements
        println(element)
    end
    return nothing
end