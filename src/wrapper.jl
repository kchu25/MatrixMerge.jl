Base.@kwdef mutable struct SEARCH_setup{T <: Integer, S <: Real}
    # hyperparameters for greedy alignment
    pval_thresh::S=2.7e-4                          # search p-value used for defining score threshold
    e_value_thresh_1::S=5e-1                        # e-value used to eliminate extra motifs
    e_value_thresh_2::S=1e-3                       # e-value used to eliminate extra motifs
    max_score_iter::T=20                           # maximum number of iterations to maximize likelihood ratio
    ic_extension_thresh::S=0.075                   # if the expanded column's ic is lower than this, expand
    mask_len::Vector{T}=[8,4,2,1]                  # mask for trimming
    mask_thresh::Vector{S}=[.625,.525,.345,.2]    # mask threshold for trimming
    perc_no_less_than::S=0.95                       # percentage tolerance in case of false captures
    entropy_score_tol::S=1e-1                      # entropy score tolerance
    pseudocounts_r::S=.0295                         # 
    pseudocounts::S=5                              # number pseudocounts added after each msa expansion
    allr_fraction_of_columns::S=0.85               # highest scoring floor(fraction) of columns to take into account for allr 
    allr_score_each_column::S=0.95               # score for each column for allr
    bg::Vector{S}=[.25,.25,.25,.25]                # background model for PWMs
    default_max_iter::T=6                          # default max iteration for maximizing the entropy scores    
    pfm_count_gpu::T=32                            # use gpu to scan if we have number of pfms higher than this number
end

mutable struct good_stuff{T <: Integer, S <: Real}
    data::Union{Sim_DNA{T,S}, FASTA_DNA{S}}
    search::SEARCH_setup{T,S}    
    ms::Union{Nothing, motifs}                      # motifs
    ms_bg::Union{Nothing, motifs}                   # copy of the motifs to scan the background for fisher exact tests

    function good_stuff{T,S}(data::Union{Sim_DNA{T,S}, FASTA_DNA{S}},
                               search_setup::SEARCH_setup{T, S}
                               ) where {T <: Integer, S <: Real}
        new(data, search_setup, nothing, nothing)
    end
end
