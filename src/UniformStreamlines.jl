module UniformStreamlines

using FastInterpolations
using LinearAlgebra: norm
using Random: shuffle!

include("Containers.jl")
include("Tracer.jl")  
include("Stream.jl")   

export stream                                         # low-level core
export evenstream                                     # high-level entry point
export colorize, streamarrows                         # post-processing
export StreamlineData, ArrowData                      # result types


end
