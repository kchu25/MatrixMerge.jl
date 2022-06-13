module MatrixMerge

using CUDA, DataStructures, DoubleFloats, 
      StatsBase, HypothesisTests, SimDNA,
      FastaLoader

export good_stuff, SEARCH_setup, find_motif

include("SEARCH_setup.jl")
include("wrapper.jl")
include("constants.jl")
include("search/helpers.jl")
include("search/PWM_Touzer.jl")
include("search/scan.jl")
include("search/extend.jl")
include("search/trim.jl")
include("search/e_value_filter.jl")
include("search/allr.jl")
include("find.jl")





end
