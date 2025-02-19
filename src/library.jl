ω_default = collect(logrange(1e05,1e-03,7*Int(log10(1e05/1e-03))+1))

abstract type CircuitElement end

mutable struct Resistor <: CircuitElement
    R       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
end
Resistor(ω=ω_default,R=1.0) = Resistor(R, ω, R*ones(length(ω)))

mutable struct Capacitor <: CircuitElement
    C       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
end
Capacitor(ω=ω_default,C=1.0) = Capacitor(C, ω, -im./(C*ω))

mutable struct Inductor <: CircuitElement
    L       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
end
Inductor(ω=ω_default,L=1.0) = Inductor(L,ω,im*L*ω)

mutable struct CPE <: CircuitElement
    Q       ::Float64
    n       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
end
CPE(ω=ω_default,Q=1.0,n=0.8) = CPE(Q,n,ω,(im*Q*ω).^(-n))

mutable struct Warburg <: CircuitElement
    type    ::String
    A       ::Float64
    B       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
end
Warburg(ω=ω_default,type="short",A=1.0,B=1.0) = begin
    if type=="short"
        Warburg(type,A,B,ω,A*im^(-0.5)*ω.^(-0.5) .* tanh.(B*im^0.5*(ω).^0.5))
    elseif type=="open"
        Warburg(type,A,B,ω,A*im^(-0.5)*ω.^(-0.5) .* coth.(B*im^0.5*(ω).^0.5))
    else
        println("Please specify either \"short\" or \"open\"")
        return nothing
    end
end
Warburg(type::String) = Warburg(ω_default,type)

mutable struct Circuit <: CircuitElement
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
end


#operator functions
import Base.(-)
function -(a::CircuitElement,b::CircuitElement)
    return Circuit(a.ω,  a.Z+b.Z)
end

import Base.(/)
function /(a::CircuitElement,b::CircuitElement)
    return Circuit(a.ω, (a.Z .*b.Z)./(a.Z+b.Z))
end


import Base.*
function *(a::Real, circuit::Circuit)
    return Circuit(circuit.ω,a*circuit.Z)
end
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
    return Warburg(w.ω,w.type,w.A,a)
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