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