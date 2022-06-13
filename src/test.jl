push!(LOAD_PATH, "/home/shane/.julia/dev/")
using CUDA 
using SimDNA
using FastaLoader, DataStructures, DoubleFloats, Statsbase

include("src/SEARCH_setup.jl")
include("src/wrapper.jl")
include("src/search/helpers.jl")
include("src/search/PWM_Touzer.jl")
include("src/search/scan.jl")
include("src/search/extend.jl")
include("src/search/trim.jl")
include("src/search/e_value_filter.jl")
include("src/search/allr.jl")



println("hi")