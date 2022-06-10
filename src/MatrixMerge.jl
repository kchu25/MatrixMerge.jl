module MatrixMerge

using CUDA 

include("SEARCH_setup.jl")
include("wrapper.jl")
include("search/helpers.jl")
include("search/PWM_Touzer.jl")
include("search/scan.jl")
include("search/extend.jl")
include("search/trim.jl")
include("search/e_value_filter.jl")
include("search/allr.jl")


end
