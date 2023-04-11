module ResolventAnalysis

# TODO: should replace with recipes to avoid heavy dependency (really want to move mode stuff and plotting to interface package)
using LinearAlgebra, TSVD, Plots

export FullSVD, TruncateSVD
export svd

export DivideAndConquer, QRIteration, Lanczos, Adaptive

export ChannelMode, test_modeinterface, project!

include("utils.jl")
include("modes.jl")

end
