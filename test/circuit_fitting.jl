using EISAnalysis, Test
using EISAnalysis: unflatten_parameters,get_params

eval(initialize())

@testset "circuit_fit" begin
    randles_circuit = 0.23r-(r-0.025wo^80)/0.2c
    randles_circuit ~ [1.0,2.0,3.0]
    # p = get_params(randles_circuit)
    # @test p == [0.23,1.0,(0.025,80),0.2]
    
    # tuples = findall(x->typeof(x)!=Float64,get_params(circuit))
    # pflat = collect(Iterators.flatten(p))
    # p_shaped = unflatten_parameters(pflat,tuples)
    # @test p_shaped = [0.23,1.0,0.025,80,0.2]
end