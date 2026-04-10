module ContactMatrices

export Groupings, ContactMatrix
export make_symmetric, reduce_groups
export setting, groupings, ndimgroups

include("groupings.jl")
include("types.jl")
include("constructors.jl")
include("show.jl")
include("symmetry.jl")
include("operations.jl")

end # module ContactMatrices
