# include("library.jl")
ws,wo = Warburg("short"), Warburg("open")
r,c,q,l = Resistor(),Capacitor(),CPE(),Inductor()

mutable struct Circuit
    ω           ::Vector{Real}
    Z           ::Vector{ComplexF64}
    elements    ::Vector
    operators   ::Vector
    order       ::Vector
    subcircuits ::Vector
end
Circuit() = Circuit([],[],[],[],[],[])

get_params(a::Resistor) = a.R
get_params(a::Capacitor) = a.C
get_params(a::Inductor) = a.L
get_params(a::CPE) = (a.Q,a.n)
get_params(a::Warburg) = (a.A,a.B)
get_params(a::Circuit) = [get_params(eval(element)) for element in a.elements]

set_params(a::Resistor,p) = p*r
set_params(a::Capacitor,p) = p*c
set_params(a::Inductor,p) = p*l
set_params(a::CPE,p) = p[1]*q^p[2]
set_params(a::Warburg,p) = (a.type=="short") ? p[1]*ws^p[2] : p[1]*wo^p[2]

#this function is sloppy
function set_params(a::Circuit,p) 
    b = Circuit(a.ω,a.Z,Vector(undef,length(a.elements)),a.operators,a.order,a.subcircuits)
    for i in eachindex(a.elements)
        element = set_params(eval(a.elements[i]),p[i])
        b.elements[i] = get_symbol(element)
    end
    b = rebuild(b)
    return b
end

#these functions are to generate symbols from Circuit elements
function get_symbol(a::Resistor)
    return :($(a.R)*r)
end
function get_symbol(a::Capacitor)
    return :($(a.C)*c)
end
function get_symbol(a::Inductor)
    return :($(a.L)*l)
end
function get_symbol(a::CPE)
    return :($(a.Q)*q^$(a.n))
end
function get_symbol(a::Warburg)
    if a.type =="short"
        return :($(a.A)*ws^$(a.B))
    elseif a.type =="open"
        return :($(a.A)*wo^$(a.B))
    end
end

#this code should be optimized. I think it's too slow and clunky right now
function get_subcircuit(subelements,suboperators,suborder)
    suborder = suborder[2:end]
    suboperators = suboperators[2:end]
    #sorted_index
    si = sortperm(suborder)
    #initialize subcircuit
    if suboperators[si[1]]==-
        # subcircuit = eval(:( $(subelements[si[1]]) - $(subelements[si[1]+1]) ) )
        subcircuit = :( $(subelements[si[1]]) - $(subelements[si[1]+1]) ) 
    elseif suboperators[si[1]]==/
        # subcircuit = eval(:( $(subelements[si[1]]) / $(subelements[si[1]+1]) ) )  
        subcircuit = :( $(subelements[si[1]]) / $(subelements[si[1]+1]) ) 
    end
    #rest of the subcircuit
    for i in si[2:end]
        if suboperators[i]==-
            ##A QUICK FIX BUT NOT ELEGANT
            if i==1
                # subcircuit = eval(:( $subcircuit - $(subelements[i]) ) )
                subcircuit = :( $subcircuit - $(subelements[i]) ) 
            else
                # subcircuit = eval(:( $subcircuit - $(subelements[i+1]) ) )
                subcircuit = :( $subcircuit - $(subelements[i+1]) ) 
            end
        elseif suboperators[i]==/
            # subcircuit = eval(:( $subcircuit / $(subelements[i+1]) ) )  
            subcircuit = :( $subcircuit / $(subelements[i+1]) )
        end
    end
    # return subcircuit
    return eval(subcircuit)
end

function rebuild(circuit)
    fullcircuit = undef
    operators = vcat(undef,circuit.operators)
    order = vcat(maximum(circuit.order),circuit.order)
    for i in 1:maximum(circuit.subcircuits)
        sub_i = findall(isequal(i),circuit.subcircuits)
        subelements = circuit.elements[sub_i]
        suboperators = operators[sub_i]
        suborder = order[sub_i]
        subcircuit = get_subcircuit(subelements,suboperators,suborder)
        if i ==1
            fullcircuit = subcircuit
        else
            fullcircuit = eval(:( $fullcircuit - $subcircuit) )
        end
    end
    fullcircuit.elements = circuit.elements
    fullcircuit.operators = circuit.operators
    fullcircuit.order = circuit.order
    fullcircuit.subcircuits = circuit.subcircuits
    return fullcircuit
end