#
# Correctness Tests
#
using EISAnalysis, Test

my_tests = ["circuit_fitting.jl", "operators.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end