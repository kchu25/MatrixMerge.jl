mutable struct motifs{T <: Integer, S <: Real}
    pfms::Vector{Matrix{S}}
    pwms::Union{Nothing,Vector{Matrix{S}}}
    thresh::Union{Nothing, Vector{S}}
    lens::Vector{T}
    num_motifs::T
    positions::Union{Nothing, Vector{Dict{T, T}}}
    scores::Union{Nothing, Vector{Dict{T, S}}}
    use_comp::Union{Nothing, Vector{Dict{T, Bool}}}
end

function init_motifs(pfms::Vector{Matrix{S}}) where {S <: Real}
    # pfms is a vector of position frequency matrices
    num_pfms = length(pfms);
    pwms = get_pwms(pfms, g.search.bg);
    thresh = get_thresh(pwms, g.search.pval_thresh, S);
    lens = type_int.([size(pfm,2) for pfm in pfms]);
    g.ms = motifs(
        pfms,
        pwms,
        thresh,
        lens,
        type_int(num_pfms),
        nothing, 
        nothing,
        nothing 
    );
end