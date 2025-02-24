# const ws,wo = Warburg("short"), Warburg("open")
# const r,c,q,l = Resistor(),Capacitor(),CPE(),Inductor()

"""
    initialize()

Generates all the circuit elements for ease of use in building circuits
"""
function initialize()
    eval(quote
    ws,wo = Warburg("short"), Warburg("open")
    r,c,q,l = Resistor(),Capacitor(),CPE(),Inductor()  
    end)
    print("Circuit elements r, c, q, l, wo, and ws generated")
    return nothing
end
"""
    Circuit

Struct for storing information about circuits

# Attributes
- `ω`: EIS frequencies
- `Z`: EIS Impedances
- `elements`: List of circuit elements. Also stores paramter information.
- `operators`: List of operators (-,/) between elements
- `order`: Order of operations
- `subcircuits`: List of subcircuits
"""
mutable struct Circuit
    ω           ::Vector{Real}
    Z           ::Vector{ComplexF64}
    elements    ::Vector
    operators   ::Vector
    order       ::Vector
    subcircuits ::Vector
end

"""
    get_params(a::Union{CircuitElement,Circuit})

Gets the parameters for elements in a circuit.

# Attributes
- `a::Union{CircuitElement,Circuit}`: The circuit or circuit element

# Examples
```julia
julia> randles_circuit = 0.23r-(r-0.025wo^80)/0.2q;
julia> p = get_params(randles_circuit)
4-element Vector{Any}:
 0.23
 1.0
  (0.025, 80.0)
  (0.2, 0.8)
```
"""
get_params(a::Resistor) = a.R
get_params(a::Capacitor) = a.C
get_params(a::Inductor) = a.L
get_params(a::CPE) = (a.Q,a.n)
get_params(a::Warburg) = (a.A,a.B)
get_params(a::Circuit) = [get_params(eval(element)) for element in a.elements]

"""
    set_params(a::Union{CircuitElement,Circuit})

Sets the parameters for elements in a circuit. 

Currently a bit sloppy.Used in `circuit_fit`

# Attributes
- `a::Union{CircuitElement,Circuit}`: The circuit or circuit element
- `p`: The parameter list. Needs to carry tuples for elements with two parameters

# Examples
```julia
julia> circuit = r-r/q;
julia> p = [0.5,2,(0.5,0.9)]
3-element Vector{Any}:
 0.5
 2
  (0.5, 0.9)
julia> updated_circuit = EISAnalysis.set_params(circuit,p); print_circuit(updated_circuit)
0.5r
2.0r
0.5 * q ^ 0.9
```
"""
set_params(a::Resistor,p) = p*r
set_params(a::Capacitor,p) = p*c
set_params(a::Inductor,p) = p*l
set_params(a::CPE,p) = p[1]*q^p[2]
set_params(a::Warburg,p) = (a.type=="short") ? p[1]*ws^p[2] : p[1]*wo^p[2]
function set_params(a::Circuit,p) 
    b = Circuit(a.ω,a.Z,Vector(undef,length(a.elements)),a.operators,a.order,a.subcircuits)
    for i in eachindex(a.elements)
        element = set_params(eval(a.elements[i]),p[i])
        b.elements[i] = get_symbol(element)
    end
    b = rebuild(b)
    return b
end


"""
    get_symbol(a::CircuitElement)

Creates the symbol (technically expressions) for a circuit element, using globally defined variables.
    
Used for generating the elements field of a circuit.

# Attributes
- `a::CircuitElement`: The circuit element
"""
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

"""
    get_subcircuit(subelements,suboperators,suborder)

Creates a circuit from a subcircuit.

Used in `rebuild`. Currently a bit sloppy.

# Attributes
- `subelements`: elements of subcircuit
- `suboperators`: operators of subcircuit
- `suborder`: operation order of subcircuit
"""
function get_subcircuit(subelements,suboperators,suborder)
    suborder = suborder[2:end]
    suboperators = suboperators[2:end]
    #sorted_index
    si = sortperm(suborder)
    #initialize subcircuit
    if suboperators[si[1]]==-
        subcircuit = :( $(subelements[si[1]]) - $(subelements[si[1]+1]) ) 
    elseif suboperators[si[1]]==/  
        subcircuit = :( $(subelements[si[1]]) / $(subelements[si[1]+1]) ) 
    end
    #rest of the subcircuit
    for i in si[2:end]
        if suboperators[i]==-
            ##A QUICK FIX BUT NOT ELEGANT
            if i==1
                subcircuit = :( $subcircuit - $(subelements[i]) ) 
            else
                subcircuit = :( $subcircuit - $(subelements[i+1]) ) 
            end
        elseif suboperators[i]==/ 
            subcircuit = :( $subcircuit / $(subelements[i+1]) )
        end
    end
    return eval(subcircuit)
end

"""
   rebuild(circuit::Circuit)

Main function for recalculating a circuit's impedance.

Used after mutating a circuit through either ~ or set_params

# Attributes
- `circuit`: Mutated circuit to be rebuilt
"""
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