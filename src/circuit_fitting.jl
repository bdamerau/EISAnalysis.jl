
#shaping the parameters to map onto circuit elements 
function unflatten_parameters(pflat,tuples)
    #initialize
    newp = [pflat[i] for i in 1:tuples[1]-1]
    for i in eachindex(tuples[1:end-1])
        #first add tuple 
        newp = vcat(newp,(pflat[tuples[i]],pflat[tuples[i]+1]))
        #then go until next tuple
        newp = vcat(newp,[pflat[j] for j in tuples[i]+2:tuples[i+1]-1])
    end
    #end
    newp = vcat(newp,(pflat[tuples[end]],pflat[tuples[end]+1]))
    return vcat(newp,[pflat[j] for j in tuples[end]+2:length(pflat)])
end

function circuit_fit(circuit, ω_data,Z_data)

    #getting initial parameters and the parameter shape
    tuples = findall(x->typeof(x)!=Float64,get_params(circuit))
    tuples += collect(0:length(tuples)-1)
    unflatten = x->unflatten_parameters(x,tuples)

    p = get_params(circuit)
    p0 = collect(Iterators.flatten(p))
    if p0 == p
        unflatten = x->x
    end
    function fitting_function(circuit,ω,p)
        p_shaped = unflatten(p)
        Z = (set_params(circuit,p_shaped)~ω).Z
        return vcat(real(Z),imag(Z))
    end
    fit_funct = (ω,p)->fitting_function(circuit,ω,p)
    fit = curve_fit(fit_funct, ω_data, vcat(real(Z_data),imag(Z_data)), p0;
        lower = zeros(length(p0)))
    p_fit = unflatten(round.(fit.param,sigdigits = 5))
    circuit_fit = set_params(circuit,p_fit)
    println("Parameters")
    println("__________")
    print(circuit_fit.elements)
    circuit_fit = circuit_fit ~ ω_data
    Z_fit = circuit_fit.Z
    plt = scatter(Z_data,label = "data")
    scatter!(plt,Z_fit,markersize=3,label = "fit")
    plot!(plt,aspect_ratio=:equal,yflip=true,ylabel = "-Im Z / Ω", xlabel = "Re Z / Ω",legend = :topleft)
    display(plt)
    return circuit_fit
end

