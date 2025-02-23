"""
Overloaded Operators: -, /, *, ^, and ~

Function /
    Description
    -----------
    Holds the inputs in parallel and generates a Circuit.
    Operates over Circuits and CircuitElements

Function -
    Description
    -----------
    Holds the inputs in series and generates a Circuit.
    Operates over Circuits and CircuitElements

Function *
    Description
    -----------
    Mutates the impedance parameter of CircuitElements

Function ^
    Description
    -----------
    Mutates the exponent parameter of CPE's and Warburgs

Function ~
    Description
    -----------
    Maps the frequencies over which impedance is calculated to the desired frequency.
    Operates over Circuits and CircuitElements.
"""

import Base.(-)
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

import Base.(/)
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


import Base.*
function *(a::Real, r::Resistor)
    return Resistor(r.ω,a*r.R)
end
function *(a::Real, c::Capacitor)
    return Capacitor(c.ω,a*c.C)
end
function *(a::Real, q::CPE)
    return CPE(q.ω,a*q.Q,q.n)
end
function *(a::Real, l::Inductor)
    return Inductor(l.ω,a*l.L)
end
function *(a::Real, w::Warburg)
    return Warburg(w.ω,w.type,a*w.A,w.B)
end

import Base.(^)
function ^(q::CPE,a::Real)
    return CPE(q.ω,q.Q,a)
end
function ^(w::Warburg,a::Real)
    return w = Warburg(w.ω,w.type,w.A,a)
end

##mapping circuit element to desired ω
import Base.(~)
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