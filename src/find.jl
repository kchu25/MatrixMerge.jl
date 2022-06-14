function get_ssc(g::good_stuff) 
    ssc = sum([sum([g.ms.scores[k][n] 
                                for n in keys(g.ms.scores[k])])
                                for k = 1:g.ms.num_motifs]) 
    # println("SSC: $(round(ssc,digits=1)), # pfms: $(length(t.ms.pfms))");
    return ssc
end

function find_motif(g::good_stuff)
    @time for _ = 1:g.search.default_max_iter
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

    @time begin 
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

    @time begin 
        last_ssc = Inf; cur_ssc = nothing;
        for j = 1:g.search.max_score_iter
            # println("$j-th iter")
            overlapping_scan!(g)
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

    @time begin
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
end

