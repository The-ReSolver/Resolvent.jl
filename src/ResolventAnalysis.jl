module ResolventAnalysis

using LinearAlgebra
using TSVD

export FullSVD, TruncateSVD
export svd

export DivideAndConquer, QRIteration, Lanczos

export ChannelMode, test_modeinterface, project!

include("utils.jl")
include("modes.jl")

end
