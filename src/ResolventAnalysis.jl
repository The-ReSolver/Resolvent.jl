module ResolventAnalysis

using LinearAlgebra, TSVD

export FullSVD, TruncateSVD
export svd

export DivideAndConquer, QRIteration, Lanczos, Adaptive

export ChannelMode, test_modeinterface, project!

include("utils.jl")
include("modes.jl")

end
