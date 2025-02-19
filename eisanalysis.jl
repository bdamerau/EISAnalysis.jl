include("library.jl")
include("functions.jl")
include("drt.jl")

#THIS MIGHT BE BAD CODE JUST DOING TO SAVE TIME LOLOLOL
ws,wo = Warburg("short"), Warburg("open")
r,c,q,l = Resistor(),Capacitor(),CPE(),Inductor()
print(" ")
