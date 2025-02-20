ω_default = collect(logrange(1e05,1e-03,7*Int(log10(1e05/1e-03))+1))

abstract type CircuitElement end

mutable struct Resistor <: CircuitElement
    R       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
    symb    
end
Resistor(ω=ω_default,R=1.0,symb = :r) = Resistor(R, ω, R*ones(length(ω)),symb)

mutable struct Capacitor <: CircuitElement
    C       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
    symb    
end
Capacitor(ω=ω_default,C=1.0,symb = :c) = Capacitor(C, ω, -im./(C*ω),symb)

mutable struct Inductor <: CircuitElement
    L       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
    symb    
end
Inductor(ω=ω_default,L=1.0,symb = :l) = Inductor(L,ω,im*L*ω,symb)

mutable struct CPE <: CircuitElement
    Q       ::Float64
    n       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
    symb    
end
CPE(ω=ω_default,Q=1.0,n=0.8,symb = :q) = CPE(Q,n,ω,(im*Q*ω).^(-n),symb)

mutable struct Warburg <: CircuitElement
    type    ::String
    A       ::Float64
    B       ::Float64
    ω       ::Vector{Real}
    Z       ::Vector{ComplexF64}
    symb    
end
Warburg(ω=ω_default,type="short",A=1.0,B=1.0) = begin
    if type=="short"
        Warburg(type,A,B,ω,A*im^(-0.5)*ω.^(-0.5) .* tanh.(B*im^0.5*(ω).^0.5),:ws)
    elseif type=="open"
        Warburg(type,A,B,ω,A*im^(-0.5)*ω.^(-0.5) .* coth.(B*im^0.5*(ω).^0.5),:wo)
    else
        println("Please specify either \"short\" or \"open\"")
        return nothing
    end
end
Warburg(type::String) = Warburg(ω_default,type)

mutable struct Circuit <: CircuitElement
    ω           ::Vector{Real}
    Z           ::Vector{ComplexF64}
    elements    ::Vector
    operators   ::Vector
    order       ::Vector
    subcircuits ::Vector
end
# Circuit(ω,Z) = Circuit(ω,Z,[],[])

#function for updating elements and operators of a Circuit
function update_circuit(a::CircuitElement,b::CircuitElement,op)
    if a isa Circuit && b isa Circuit
        elements = vcat(a.elements,b.elements)
        operators = vcat(a.operators,op,b.operators)
        order = vcat(a.order,maximum(a.order)+maximum(b.order)+1,b.order .+ maximum(a.order))
        subcircuits = vcat(a.subcircuits,b.subcircuits .+ maximum(a.subcircuits))
    elseif a isa Circuit && !(b isa Circuit)
        elements = vcat(a.elements,b.symb)
        operators = vcat(a.operators,op)
        order = vcat(a.order,maximum(a.order)+1)
        subcircuits = vcat(a.subcircuits,maximum(a.subcircuits))
    elseif !(a isa Circuit) && b isa Circuit
        elements = vcat(a.symb,b.elements)
        operators = vcat(op,b.operators)
        order = vcat(maximum(b.order)+1,b.order)
        subcircuits = vcat(maximum(b.subcircuits),b.subcircuits)
    elseif !(a isa Circuit) && !(b isa Circuit)
        elements = [a.symb,b.symb]
        operators = [op]
        order = [1]
        subcircuits = [1,1]
    end
    return elements,operators,order,subcircuits
end

#operator functions
import Base.(-)
function -(a::CircuitElement,b::CircuitElement)
    elements,operators,order,subcircuits = update_circuit(a,b,-)
    return Circuit(a.ω, a.Z+b.Z, elements,operators,order,subcircuits)
end

import Base.(/)
function /(a::CircuitElement,b::CircuitElement)
    elements,operators,order,subcircuits = update_circuit(a,b,/)
    return Circuit(a.ω, a.Z.*b.Z ./(a.Z+b.Z), elements,operators,order,subcircuits)
end


import Base.*
function *(a::Real, circuit::Circuit)
    return Circuit(circuit.ω,a*circuit.Z,circuit.elements,
        circuit.operators,circuit.order,circuit.subcircuits) ##this should be either adjusted or deleted
end
function *(a::Real, r::Resistor)
    return Resistor(r.ω,a*r.R,:($a*r))
end
function *(a::Real, c::Capacitor)
    return Capacitor(c.ω,a*c.C,:($a*c))
end
function *(a::Real, q::CPE)
    if q.n == 0.8
        return CPE(q.ω,a*q.Q,q.n,:($a*q))
    else
        return CPE(q.ω,a*q.Q,q.n,:($a*q^$(q.n)))
    end
end
function *(a::Real, l::Inductor)
    return Inductor(l.ω,a*l.L,:($a*l))
end
function *(a::Real, w::Warburg)
    if w.B == 1.0
        w = Warburg(w.ω,w.type,a*w.A,w.B)
        w.symb = :($a*$(w.symb))
    else
        w =  Warburg(w.ω,w.type,a*w.A,w.B)
        w.symb = :($a*$(w.symb)^$(w.B))
    end
    return w
end

import Base.(^)
function ^(q::CPE,a::Real)
    if q.Q == 1
        return CPE(q.ω,q.Q,a,:(q^$a))
    else
        return CPE(q.ω,q.Q,a,:($(q.Q)*q^$a))
    end
end
function ^(w::Warburg,a::Real)
    if w.A == 1.0
        w = Warburg(w.ω,w.type,w.A,a)
        w.symb = :($(w.symb)^$a)
    else
        w = Warburg(w.ω,w.type,w.A,a)
        w.symb = :($(w.A)*$(w.symb)^$(a))
    end
    return w
end

##mapping circuit element to desired ω
import Base.(~)
function ~(a::Resistor,ω::Vector{Float64})
    return Resistor(ω,a.R,a.symb)
end
function ~(a::Capacitor,ω::Vector{Float64})
    return Capacitor(ω,a.C,a.symb)
end
function ~(a::CPE,ω::Vector{Float64})
    return CPE(ω,a.Q,a.n,a.symb)
end
function ~(a::Inductor,ω::Vector{Float64})
    return Inductor(ω,a.L,a.symb)
end
function ~(a::Warburg,ω::Vector{Float64})
    w = Warburg(ω,a.type,a.A,a.B)
    w.symb = a.symb
    return w
end

###this function needs work. It technically works, but I think I'll have to adjust 
###the functions above
function ~(a::Circuit,ω::Vector{Float64})
    for i in eachindex((a.elements))
        element = eval(a.elements[i]) ~ ω
        # element.symb = :($element )
        # a.elements[i] = element.symb
        a.elements[i] = element
        println(a.elements[i])
    end
    return rebuild(a)
end