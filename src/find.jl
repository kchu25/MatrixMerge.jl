function get_ssc(g::good_stuff) 
    ssc = sum([sum([g.ms.scores[k][n] 
                                for n in keys(g.ms.scores[k])])
                                for k = 1:g.ms.num_motifs]) 
    # println("SSC: $(round(ssc,digits=1)), # pfms: $(length(t.ms.pfms))");
    return ssc
end

function find_motif!(g::good_stuff)
    for _ = 1:g.search.default_max_iter
        scan_w_gpu!(g, g.data.data_matrix_gpu; re_evaluate_pfm=false,
                                               re_evaluate_pwm=false,
                                               re_evaluate_thresh=false);
        extend!(g; re_evaluate_pfm=true,
                   re_evaluate_pwm=false,
                   re_evaluate_thresh=false,
                   smoothing=true);
        masked_trim!(g; re_evaluate_pfm=false,
                        re_evaluate_pwm=true,
                        re_evaluate_thresh=true);
    end

    begin 
        filter_using_evalue!(g);
        allr_merge!(g; re_evaluate_pwm=false,
                       re_evaluate_thresh=false);
        masked_trim!(g; re_evaluate_pfm=false,
                        re_evaluate_pwm=true,
                        re_evaluate_thresh=true,
                        update_positions=false);
    end

    # set up a higher ic threshold as most extensions is done beforehand
    g.search.ic_extension_thresh=0.6;

    # println(" motifs --- $(length(g.ms.pfms))")

    @time begin 
        last_ssc = Inf; cur_ssc = nothing;
        for j = 1:g.search.max_score_iter
            # println("$j-th iter")
            # println("yyoyoo  111?")
            overlapping_scan!(g)
            # println("yyoyoo 222?")
            non_overlap_scan!(g; re_evaluate_pfm=false, 
                                re_evaluate_pwm=false, 
                                re_evaluate_thresh=false)
            extend!(g; re_evaluate_pfm=true,
                    re_evaluate_pwm=false,
                    re_evaluate_thresh=false,
                    smoothing=false)
            masked_trim!(g; re_evaluate_pwm=true,
                            re_evaluate_thresh=true)
            cur_ssc = get_ssc(g);
            abs(last_ssc - cur_ssc) < g.search.entropy_score_tol && break;
            last_ssc = cur_ssc;        
        end
    end
    # println(" motifs --- $(length(g.ms.pfms))")

    begin
        last_ssc = Inf; cur_ssc = nothing;

        filter_using_evalue!(g; cpu=true, non_overlap=true);

        allr_merge!(g)
        rid_of_len_less_than_six!(g)
        for _ = 1:g.search.max_score_iter
            overlapping_scan!(g)
            non_overlap_scan!(g;re_evaluate_pfm=true, 
                                re_evaluate_pwm=true, 
                                re_evaluate_thresh=true,
                                smoothing=false,
                                less_pseudocount=true)
            cur_ssc = get_ssc(g);
            abs(last_ssc - cur_ssc) < g.search.entropy_score_tol && break;
            last_ssc = cur_ssc;
        end
    end
    @info "Found $(length(g.ms.pfms)) motif(s)"
    # println(" motifs --- $(length(g.ms.pfms))")
end


function try_to_find_motif(filters, fil_size, data, target_folder_expr; 
                           number_trials=10,
                           pval_thresh_inc=5e-4,
                           eval_thresh_inc_1=2.5e-1,
                           eval_thresh_inc_2=1e-2,
                           simulated_data=false)                              
    g=nothing; motif_found = false; 
    pval_thresh = dat_t(2.7e-4);     
    eval_thresh1 = dat_t(5e-1);                                       
    eval_thresh2 = dat_t(1e-3);
    num_trial_left = number_trials;
    while(num_trial_left > 0)
        try
            g = good_stuff{int_t,dat_t}(data, 
                                        SEARCH_setup{int_t, dat_t}(),
                                        filters, fil_size
                                        );
            g.search.pval_thresh = pval_thresh;
            g.search.e_value_thresh_1 = eval_thresh1;
            g.search.e_value_thresh_2 = eval_thresh2;                                        
            find_motif!(g)
            motif_found = true;
            break
        catch e
            if isa(e, ArgumentError)
                pval_thresh += pval_thresh_inc;
                eval_thresh1 += eval_thresh_inc_1;
                eval_thresh2 += eval_thresh_inc_2;
                @info "Did not find any motif..."
                @info "Relax search p-value to $(round(pval_thresh,digits=4)) for the score threshold"
                @info "Relax search e-value cutoff-1 to $(round(eval_thresh1,digits=4))"
                @info "Relax search e-value cutoff-2 to $(round(eval_thresh2,digits=4))"
                num_trial_left -= 1;
            else
                break
            end        
        end            
    end
    # println("hi 1")
    if simulated_data
        # for simulated data, if couldn't find anything, as least save the ground truth
        @save target_folder_expr*"/gt_motif.jld2" g.data.motif
        motif_found && save_found_results_sim(target_folder_expr, g);
    else
        if motif_found
            save_result_fasta(g, target_folder_expr);
        end
    end
end


# this is just for testing jarpar motifs
function try_to_find_motif_jaspar(filters, fil_size, data, 
                           source_p_v_out, expr_logos_p_v_transfac, ref_logo_where, ref_save_where, matrix_name)
    g=nothing; motif_found = false; 
    pval_thresh = dat_t(2.7e-4);
    eval_thresh1 = dat_t(5e-1);
    eval_thresh2 = dat_t(1e-3);

    try
        g = good_stuff{int_t,dat_t}(data, 
                                    SEARCH_setup{int_t, dat_t}(),
                                    filters, fil_size
                                    );
        g.search.pval_thresh = pval_thresh;
        g.search.e_value_thresh_1 = eval_thresh1;
        g.search.e_value_thresh_2 = eval_thresh2;                                        
        find_motif!(g);
        length(g.ms.pfms) > 0  && (motif_found = true);
    catch e
        if isa(e, ArgumentError)
            println("Catched argument error; no motif found in this run.")
        end
    end

    return motif_found, save_result_fasta_jaspar(g, source_p_v_out, expr_logos_p_v_transfac, ref_logo_where, ref_save_where, matrix_name);
end
