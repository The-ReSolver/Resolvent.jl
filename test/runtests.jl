using ResolventAnalysis
using Test
using Random
using LinearAlgebra

# set up the PyCall interface to be able to use Sean's Matlab code directly
using PyCall
o = pyimport("oct2py").octave
o.addpath("./primitive_resolvent")

using ChebUtils

include("test_utils.jl")
# include("test_modes.jl")
