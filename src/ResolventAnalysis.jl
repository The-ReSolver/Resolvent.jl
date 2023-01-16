module ResolventAnalysis

using LinearAlgebra

export FullSVD, TruncatedSVD
export Resolvent

include("utils.jl")
include("resolvent.jl")

end
