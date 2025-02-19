EISAnalysis.jl
==============

This is a small package for working with Electrochemical Impedance Spectroscopy (EIS) data 
from electrochemical systems. It is currently designed to do the following:
1. Building and fitting custom Equivalent Circuit Models (ECMs) to EIS data
2. Computing Distribution of Relaxation Times (DRTs) from EIS data

This package defines a `CircuitElement` type and exports several circuit elements as variables. 
Circuits can be built quickly using overloaded base operators. In the following example, a 
Randles Circuit is constructed and its DRT function is computed.
```julia
randles_circuit = 0.23r-(r-0.025ws^80)/0.2q
randles_fit = compute_drt(randles_circuit.ω,randles_circuit.Z)
```
![test](https://github.com/user-attachments/assets/ddc77786-f392-4297-b9bd-f040b3ca7caf)
