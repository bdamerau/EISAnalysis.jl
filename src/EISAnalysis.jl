module EISAnalysis

export CircuitElement,Resistor,Capacitor,CPE,Inductor,Warburg,Circuit
export r,c,q,l,wo,ws

export plot_Nyquist,plot_Nyquist!,fitted_circuit

export compute_drt

using Statistics,LinearAlgebra,LsqFit
using Plots

include("circuits.jl")
include("drt.jl")
include("functions.jl")
include("library.jl")
include("symbolic_functions.jl")

end # module EISAnalysis
