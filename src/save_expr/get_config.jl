function map_key_update(selected_keys, sort_perm_map)
    selected_keys_update = typeof(selected_keys)();
    for i = 1:length(selected_keys)
        push!(selected_keys_update, [(sort_perm_map[j[1]], j[2]) for j in selected_keys[i]]);
    end
    selected_keys_update
end

function get_configurations(g::Union{good_stuff, Nothing})
    configurations = Dict{Vector{Tuple{Int, Bool}}, Int}();
    gap_configurations = Dict{Vector{Tuple{Int, Bool}}, Vector{Vector{Int}}}();
    for n = 1:g.data.N
        occured_n = Tuple{Int, Bool}[];
        position_n = Int[];
        for k = 1:g.ms.num_motifs
            if haskey(g.ms.positions[k],n)
                push!(occured_n, (k, g.ms.use_comp[k][n]));
                # push!(occured_n, (sort_perm_indices[k], g.ms.use_comp[k][n])); # since we've permuted by using e-values
                push!(position_n, g.ms.positions[k][n])            
            end
        end

        if !isempty(occured_n)
            sort_perm = sortperm(position_n);
            occured_n = occured_n[sort_perm];
            if haskey(configurations, occured_n) 
                configurations[occured_n] += 1;
            else configurations[occured_n] = 1;
            end
            if length(occured_n) > 1
                position_n = position_n[sort_perm];
                gap_vec = Int[];
                num_gaps = length(occured_n)-1;
                for i = 1:num_gaps
                    # which_one = inv_sort_perm_map[occured_n[i][1]];
                    gap_len = position_n[i+1]-(position_n[i]+g.ms.lens[occured_n[i][1]]-1);
                    push!(gap_vec, gap_len)
                end
                @assert all(gap_vec .≥ 0) "gap between PWM positions should not be negative"
                if haskey(gap_configurations, occured_n)
                    push!(gap_configurations[occured_n], gap_vec);
                else gap_configurations[occured_n] = [gap_vec];
                end
            end
        end
    end
    return configurations, gap_configurations
end

function make_configurations_plots_2(configurations, 
                                   gap_configurations,
                                   pics_folder,
                                   sort_perm_map; 
                                   threshold_weight=0.001,
                                   kde_bandwidth=1.25)
    config_counts = [configurations[k] for k in keys(configurations)];
    config_weights = config_counts ./ sum(config_counts);
    # take only configurations that have "significant" weights
    # might do a more rigorious way to evaluate what is significant later
    # for now, simply discard all the configuration that has weight less than threshold_weight
    selected_keys = [k for k in keys(configurations)][config_weights .> threshold_weight];
    selected_configurations = Dict(k=>configurations[k] for k in selected_keys);
    selected_gap_configurations = Dict(k=>gap_configurations[k] 
                                for k in selected_keys if length(k) > 1); 
    # re-estimate the weights                                                                      
    selected_config_counts = [selected_configurations[k] for k in selected_keys];
    selected_config_weights = selected_config_counts ./ sum(selected_config_counts);
    # sort the patterns according to their weights
    sp_ind = sortperm(selected_config_weights, rev=true);
    selected_config_weights_sorted = selected_config_weights[sp_ind];
    selected_config_counts_sorted = selected_config_counts[sp_ind];
    selected_keys_sorted = selected_keys[sp_ind];
    selected_keys_sorted_update = map_key_update(selected_keys_sorted, sort_perm_map)
    
    # plot the gap kde plots
    # from the largest weight to the smallest weight
    for ind = 1:length(sp_ind)
        key = selected_keys_sorted[ind];
        for gap_bt_ind = 1:(length(key)-1)     
            # println("ind:", ind)
            # println("length key: $(length(key))")
            # remember gap_config is a vector of vector of integers...
            # [[gap_len_1, gap_len_2]₁, [gap_len_1, gap_len_2]₂, ...]
            q=[i[gap_bt_ind] for i in selected_gap_configurations[key]];
            p=Gadfly.plot(DataFrame(q=q), x=:q, 
                          Geom.density(bandwidth=kde_bandwidth),
                          Guide.xlabel("Number of nucleotides in between"), 
                          Guide.xticks(ticks=0:10:maximum(q)),
                          kde_theme);
            # save plot
            draw(PNG(pics_folder*"/gap_$(gap_bt_ind)_mode_$ind.png", 12inch, 7inch), p)
        end        
    end
    
    # plot the weights
    # patterns = [join([j[2] ? "$(j[1])" : "rc($(j[1]))" for j in i],", ") for i in keys(selected_configurations)]

    patterns = [join([!j[2] ? "D$(j[1])" : "rc(D$(j[1]))" for j in i],", ") 
                    for i in selected_keys_sorted_update]
    df=DataFrame(pattern=patterns, weights=selected_config_weights_sorted); df = sort(df, :weights);
    p = Gadfly.plot(df, y=:pattern, x=:weights,
        Geom.bar(orientation=:horizontal),
        Guide.xlabel("Weights", orientation=:horizontal),
        Guide.ylabel("    "),
        # Guide.ylabel("Binding Patterns"),
        Guide.title("Enriched Patterns and their weights"),
        tufte_bar);
    draw(PNG(pics_folder*"/mb.png", 12inch, 12inch), p)
    return selected_keys_sorted_update, 
           selected_config_weights_sorted, 
           selected_config_counts_sorted,
           maximum(length.(selected_keys_sorted_update))-1
end    