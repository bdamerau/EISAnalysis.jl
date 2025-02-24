import Base.(-)
import Base.(*)
import Base.(^)
import Base.(/)
import Base.(~)

"""
    -(a::Union{Circuit,CircuitElement},b::Union{Circuit,CircuitElement})

Holds the inputs in series and generates a Circuit.

Operates over Circuits and CircuitElements.

# Arguments
- `a::Union{Circuit,CircuitElement}`: circuit or circuit element to hold in series with `b`
- `b::Union{Circuit,CircuitElement}`: circuit or circuit element to hold in series with `a`

#Examples
```julia
julia> using EISAnalysis
julia> eval(initialize());
julia> circuit1 = r-c;
julia> circuit2 = circuit1-c;
julia> circuit2.Z == (r-c-c).Z == (r-0.5c).Z
true
```
"""
function -(a::CircuitElement,b::Circuit)
    return Circuit(a.ω, a.Z+b.Z, vcat(get_symbol(a),b.elements),vcat(-,b.operators),
    vcat(maximum(b.order)+1,b.order),vcat(maximum(b.subcircuits),b.subcircuits))
end
function -(a::Circuit,b::CircuitElement)
    return Circuit(a.ω, a.Z+b.Z, vcat(a.elements,get_symbol(b)),vcat(a.operators,-),
    vcat(a.order,maximum(a.order)+1),vcat(a.subcircuits,maximum(a.subcircuits)))
end
function -(a::CircuitElement,b::CircuitElement)
    return Circuit(a.ω, a.Z+b.Z, [get_symbol(a),get_symbol(b)],[-],[1],[1,1])
end
function -(a::Circuit,b::Circuit)
    return Circuit(a.ω, a.Z+b.Z, vcat(a.elements,b.elements),vcat(a.operators,-,b.operators),
    vcat(a.order,maximum(a.order)+maximum(b.order)+1,b.order .+ maximum(a.order)),vcat(a.subcircuits,b.subcircuits .+ maximum(a.subcircuits)))
end

"""
    Base./(a::Union{Circuit,CircuitElement},b::Union{Circuit,CircuitElement})

Holds the inputs in parallel and generates a Circuit.

Operates over Circuits and CircuitElements.

# Arguments
- `a::Union{Circuit,CircuitElement}`: circuit or circuit element to hold in parallel with `b`
- `b::Union{Circuit,CircuitElement}`: circuit or circuit element to hold in parallel with `a`

#Examples
```julia
julia> using EISAnalysis
julia> eval(initialize());
julia> circuit1 = r/c;
julia> circuit1.elements
2-element Vector{Expr}:
 :(1.0r)
 :(1.0c)
julia> circuit2 = circuit1/c;
julia> circuit2.Z == ((r/c)/c).Z
true
```
"""
function /(a::CircuitElement,b::Circuit)
    return Circuit(a.ω, a.Z.*b.Z ./(a.Z+b.Z), vcat(get_symbol(a),b.elements),vcat(/,b.operators),
    vcat(maximum(b.order)+1,b.order),vcat(maximum(b.subcircuits),b.subcircuits))
end
function /(a::Circuit,b::CircuitElement)
    return Circuit(a.ω, a.Z.*b.Z ./(a.Z+b.Z), vcat(a.elements,get_symbol(b)),vcat(a.operators,/),
    vcat(a.order,maximum(a.order)+1),vcat(a.subcircuits,maximum(a.subcircuits)))
end
function /(a::CircuitElement,b::CircuitElement)
    return Circuit(a.ω, a.Z.*b.Z ./(a.Z+b.Z), [get_symbol(a),get_symbol(b)],[/],[1],[1,1])
end
function /(a::Circuit,b::Circuit)
    return Circuit(a.ω, a.Z.*b.Z ./(a.Z+b.Z), vcat(a.elements,b.elements),vcat(a.operators,/,b.operators),
    vcat(a.order,maximum(a.order)+maximum(b.order)+1,b.order .+ maximum(a.order)),vcat(a.subcircuits,b.subcircuits .+ maximum(a.subcircuits)))
end


"""
    Base.*(a::Real,b::CircuitElement)

Mutates the impedance parameter of CircuitElements.

# Arguments
- `a::Real`: Impedance parameter value
- `b::CircuiElement`: Circuit element 

#Examples
```julia
julia> using EISAnalysis
julia> eval(initialize());
julia> r.R
1.0
julia> twor = 2r;
julia> twor.R
2.0
```
"""
function *(a::Real, b::Resistor)
    return Resistor(b.ω,a*b.R)
end
function *(a::Real, b::Capacitor)
    return Capacitor(b.ω,a*b.C)
end
function *(a::Real, b::CPE)
    return CPE(b.ω,a*b.Q,b.n)
end
function *(a::Real, b::Inductor)
    return Inductor(b.ω,a*b.L)
end
function *(a::Real, b::Warburg)
    return Warburg(b.ω,b.type,a*b.A,b.B)
end

"""
    Base.^(a::CircuitElement,b::Real)

Mutates the exponent parameter of CPE's and Warburgs.

# Arguments
- `a::CircuiElement`: circuit or circuit element to hold in series with `a`
- `b::Real`: Exponential parameter value

#Examples
```julia
julia> using EISAnalysis
julia> eval(initialize());
julia> circuit = q-wo; print_circuit(circuit)
1.0 * q ^ 0.8
1.0 * wo ^ 1.0
julia> circuit2 = q^0.6-wo^5; print_circuit(circuit2)
1.0 * q ^ 0.6
1.0 * wo ^ 5.0
```
"""
function ^(b::CPE,a::Real)
    return CPE(b.ω,b.Q,a)
end
function ^(b::Warburg,a::Real)
    return Warburg(b.ω,b.type,b.A,a)
end

"""
    Base.~(a::Union{CircuitElement,Circuit},ω::Vector{Float64})
Maps the frequencies over which impedance is calculated to the desired frequency.

Operates over Circuits and CircuitElements.

# Arguments
-`a`: Circuit element or circuit
-`ω`: Frequencies to map to

# Examples
```julia
julia> ω = [0.1,1,10]
3-element Vector{Float64}:
  0.1
  1.0
 10.0
julia> using EISAnalysis
julia> eval(initialize());
julia> circuit = r/c ~ω; println(circuit.ω,circuit.Z)
Real[0.1, 1.0, 10.0]ComplexF64[0.9900990099009901 - 0.09900990099009901im, 0.5 - 0.5im, 0.009900990099009903 - 0.09900990099009901im]
```
"""
function ~(a::Resistor,ω::Vector{Float64})
    return Resistor(ω,a.R)
end
function ~(a::Capacitor,ω::Vector{Float64})
    return Capacitor(ω,a.C)
end
function ~(a::CPE,ω::Vector{Float64})
    return CPE(ω,a.Q,a.n)
end
function ~(a::Inductor,ω::Vector{Float64})
    return Inductor(ω,a.L)
end
function ~(a::Warburg,ω::Vector{Float64})
    return Warburg(ω,a.type,a.A,a.B)
end
function ~(a::Circuit,ω::Vector{Float64})
    b = Circuit(a.ω,a.Z,Vector(undef,length(a.elements)),a.operators,a.order,a.subcircuits)
    for i in eachindex(a.elements)
        b.elements[i] = :($(a.elements[i]) ~ $ω)
    end
    return rebuild(b)
end