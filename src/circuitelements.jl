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