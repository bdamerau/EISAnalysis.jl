using EISAnalysis, Test
using EISAnalysis: Resistor,Capacitor,CPE,Inductor,Warburg

eval(initialize()) #calling r, c, l, q, wo, and ws
ω = r.ω
@testset "-" begin
    #checking to make sure - holds basic circuits and circuit elements in series
    r_series = r-r
    @test r_series.Z ≈ 2r.Z
    r_doubleseries = r_series-r_series
    @test r_doubleseries.Z ≈ 2*(r_series.Z)

    c_series = c-c
    @test c_series.Z ≈ (0.5c).Z
    c_doubleseries = c_series-c_series
    @test c_doubleseries.Z ≈ 2*(c_series.Z)

    l_series = l-l
    @test l_series.Z ≈ 2l.Z
    l_doubleseries = l_series-l_series
    @test l_doubleseries.Z ≈ 2*(l_series.Z)

    println("Overloaded operator \"-\" behaving properly" )
end
@testset "/" begin
    #checking to make sure / holds basic circuits and circuit elements in parallel
    r_parallel = r/r
    @test r_parallel.Z ≈ (0.5r).Z
    r_doubleparallel = r_parallel/r_parallel
    @test r_doubleparallel.Z ≈ 0.5*(r_parallel.Z)

    c_parallel = c/c
    @test c_parallel.Z ≈ (2c).Z
    c_doubleparallel = c_parallel/c_parallel
    @test c_doubleparallel.Z ≈ 0.5*(c_parallel.Z)

    l_parallel = l/l
    @test l_parallel.Z ≈ (0.5l).Z
    l_doubleparallel = l_parallel/l_parallel
    @test l_doubleparallel.Z ≈ 0.5*(l_parallel.Z)

    println("Overloaded operator \"/\" behaving properly" )
end
@testset "*" begin
    @test (2r).Z ≈ Resistor(ω,2*r.R).Z
    @test (2c).Z ≈ Capacitor(ω,2*c.C).Z
    @test (2l).Z ≈ Inductor(ω,2*l.L).Z
    @test (2q).Z ≈ CPE(ω,2*q.Q,q.n).Z
    @test (2wo).Z ≈ Warburg(ω,"open",2*wo.A,wo.B).Z
    println("Overloaded operator \"*\" behaving properly" )
end
@testset "^" begin
    @test (q^0.5).Z ≈ CPE(ω,q.Q,0.5).Z
    @test (wo^0.5).Z ≈ Warburg(ω,"open",wo.A,0.5).Z
    println("Overloaded operator \"^\" behaving properly" )
end
@testset "~" begin
    ω_test = [1.0,2.0,3.0]
    @test (r ~ ω_test).ω ≈ ω_test
    @test (c ~ ω_test).ω ≈ ω_test
    @test (l ~ ω_test).ω ≈ ω_test
    @test (q ~ ω_test).ω ≈ ω_test
    @test (wo ~ ω_test).ω ≈ ω_test
    test_circuit = r-(r-wo)/q
    # @test (test_circuit ~ ω_test).ω ≈ ω_test
    println("Overloaded operator \"^\" behaving properly" )
end