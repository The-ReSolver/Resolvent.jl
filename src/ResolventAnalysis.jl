module ResolventAnalysis

using LinearAlgebra

export FullSVD, TruncateSVD
export svd

export ChannelMode, test_modeinterface

include("utils.jl")
include("modes.jl")

end
