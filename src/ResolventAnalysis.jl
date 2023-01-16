module ResolventAnalysis

using LinearAlgebra

export FullSVD, TruncateSVD
export Resolvent

include("utils.jl")
include("resolvent.jl")

end
