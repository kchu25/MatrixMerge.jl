function save_found_motifs(target_folder::String, ms)
    for (ind,pfm) in enumerate(ms.pfms) 
        CSV.write(target_folder*"/d_pfm$ind.csv",  Tables.table(pfm), writeheader=false) 
    end
end


const pkg_dir = pkgdir(MatrixMerge);
const script_loc_template = pkg_dir*"/scripts_python/jinja_templates/";


function save_found_results_sim(output_folder::String, g::Union{good_stuff, Nothing})        
    if !isnothing(g)
        # save the discovered motifs for PWM plots
        save_found_motifs(output_folder, g.ms);
        # make the csv on ground truth motifs to show how they bind
        make_motif_type_file(g.data.motif, output_folder);
        # get binding info, score contribution, and perf-coeff-d; return the max lr score
        max_lr_score = gt_cover_by(g, output_folder);
        # get the performance coefficient
        perf_coeff = get_perform_coeff(g.ms, g.data.motif, g.data);        
        @info "perf coeff: $perf_coeff"
        ###### get the e-values ##########################
        g.ms.e_values = get_evalues_each(g);
        # println(g.ms.e_values...)
        ###### save the pickle ###########################
        gms = g.data.motif; dms = g.ms;
        @save output_folder*"/gt_motif.jld2" gms
        @save output_folder*"/d_motif.jld2" dms
        # @save output_folder*"/gt_data_matrix.jld2" data_matrix             
        script_loc_logo = pkg_dir*"/scripts_python/get_logos.py";
        script_loc_render = pkg_dir*"/scripts_python/render.py";  
        run(`python3 $script_loc_logo $output_folder`);
        run(`python3 $script_loc_render $output_folder $perf_coeff $max_lr_score $(g.data.prob_per_seq) $script_loc_template`);
    else
        # script_loc_render = pkg_dir*"/scripts_python/render_not_found.py";  
        # run(`python3 $script_loc_render $output_folder $(data.prob_per_seq) $script_loc_template`);
    end
end