using ResolventAnalysis
using Test
using Random
using LinearAlgebra

# set up the PyCall interface to be able to use Sean's Matlab code directly
using PyCall
o = pyimport("oct2py").octave
o.addpath("/home/tb6g16/Documents/PhD/Misc/Resolvent Codes/Primitive Variables Resolvent")

include("test_utils.jl")
include("test_resolvent.jl")
