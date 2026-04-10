module Epidemics

using ContactMatrices
using DataFrames
using OrdinaryDiffEq

export SEIR, Intervention, Vaccination, simulate

include("types.jl")
include("ode.jl")
include("output.jl")
include("simulate.jl")

end
