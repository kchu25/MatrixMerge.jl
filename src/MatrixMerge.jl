module MatrixMerge

using CUDA, DataStructures, DoubleFloats, 
      StatsBase, HypothesisTests, SimDNA,
      FastaLoader, JLD2, DataFrames, Mustache,
      CSV, Gadfly, Cairo

export good_stuff, 
       SEARCH_setup, 
       find_motif!, 
       find_and_save_motif_sim,
       save_found_results_sim,
       try_to_find_motif 

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
include("save_sim/cover_calculations.jl")
include("save_sim/save_learned_results.jl")
include("save_expr/get_config.jl")
include("save_expr/template.jl")
include("save_expr/save.jl")

include("find.jl")

# function find_and_save_motif_sim(g::good_stuff, target_folder::String)
#     find_motif!(g);
#     save_found_results_sim(target_folder, g);
# println("hi")
# end

end
