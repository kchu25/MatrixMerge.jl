function get_pfm_from_transfac(transfac_path::String; float_type=Float32)
    f = open(transfac_path, "r"); r = readlines(f); close(f);
    rows = Vector{String}();
    start = false
    for i in r
        if i[1:2] == "01" 
            start = true
        elseif i[1:2] == "XX"
            start = false
        end
        start && push!(rows, i);    
    end
    parse_counts = [parse.(float_type, i[2:end]) for i in split.(rows, "\t")];
    count_mat = reduce(hcat, parse_counts); count_mat = count_mat .+ 0.01;
    return count_mat ./ sum(count_mat, dims=1)
end

function save_pfms_as_transfac(logo_folder::String, g::good_stuff, sort_perm::Vector{Int})
    for (i,ind) in enumerate(sort_perm)
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = g.ms.pfms[ind] .* 1000; # make it a count matrix
        for j = 1:size(g.ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -o $(logo_folder)/d$(i).png`);

        # do it for the reverse complement as well
        io = open(logo_folder*"/d$(i)_c.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = pfm_complement(g.ms.pfms[ind]) .* 1000;
        for j = 1:size(g.ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 17 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);
    end
end

function save_result_fasta(g::Union{good_stuff, Nothing}, target_folder::String)
    logo_folder_name = "logos"
    pics_folder_name = "other_pics"
    logo_folder = target_folder*"/"*logo_folder_name
    pics_folder = target_folder*"/"*pics_folder_name

    if !isnothing(g)
        # make folders
        !isdir(target_folder) && (mkpath(target_folder);)
        !isdir(logo_folder) && (mkpath(logo_folder);)
        !isdir(pics_folder) && (mkpath(pics_folder);)
        ###### save msa file for each pwm ################    
        evalues = filter_using_evalue!(g; cpu=true, non_overlap=true, get_evalues_only=true);            
        sort_perm = sortperm(evalues); # sort it according to e-values (small to big)
        sort_perm_map = Dict(item=>index for (index,item) in enumerate(sort_perm));
        evalues = round.(evalues[sort_perm], sigdigits=3);
        scores = [round(sum(values(g.ms.scores[i])),digits=1)
                    for i = 1:g.ms.num_motifs][sort_perm]; # scores of each discovered motifs
        msa_counts = [length(g.ms.positions[j]) for j = 1:g.ms.num_motifs][sort_perm];

        configurations, gap_configurations = get_configurations(g);
        selected_keys_sorted_update, weights_sorted, count_sorted, max_gap_num = make_configurations_plots_2(configurations, gap_configurations, pics_folder, sort_perm_map);
        weights_sorted = round.(weights_sorted, digits=3);

        labels = ["D$j" for j = 1:g.ms.num_motifs];
        logos = ["d$(j)" for j = 1:g.ms.num_motifs];
            
        save_pfms_as_transfac(logo_folder, g, sort_perm);

        # render the main result page
        df = DataFrame(label=labels, count=msa_counts, slrs=scores, eval=evalues, logo=logos);
        out = Mustache.render(html_template, target_folder=target_folder, logo_folder=logo_folder_name, num_seq=g.data.N, DF=df);
        ###### html stuff #########################            
        io = open(target_folder*"/summary.html", "w")
        print(io, out);
        close(io)

        # render the binding pattern page
        tpl=create_template_bindind_patterns(max_gap_num)
        out = Mustache.render(tpl, DF=get_df_for_binding_pattern(selected_keys_sorted_update, 
                                                                 logo_folder_name, 
                                                                 pics_folder_name,
                                                                 weights_sorted,
                                                                 count_sorted))
        io = open(target_folder*"/bp.html", "w")
        print(io, out);
        close(io)
    end
end

function get_rounded_eval(pval::Real)
    str = "$pval";
    if !occursin("e-", str)
        return string(round(pval, sigdigits=3));
    else
        q = split(str, "e-");
        return join([q[1][1:4], q[2]], "e-")
    end
end

function save_result_fasta_jaspar(g::Union{good_stuff, Nothing}, 
                                  target_folder::String, 
                                  source_folder_logo_transfac::String,
                                  ref_logo_where::String,
                                  ref_save_where::String,
                                  matrix_name::String
                                  )
    logo_folder_name = "logos"
    pics_folder_name = "other_pics"
    logo_folder = target_folder*"/"*logo_folder_name
    pics_folder = target_folder*"/"*pics_folder_name
    
    # for summary
    top3_evalues = Vector{Real}();
    top3_logo_link = Vector{String}();
    details_link = "";    
    #############

    if !isnothing(g) && length(g.ms.pfms) != 0
        # make folders
        !isdir(target_folder) && (mkpath(target_folder);)
        !isdir(logo_folder) && (mkpath(logo_folder);)
        !isdir(pics_folder) && (mkpath(pics_folder);)
        ###### save msa file for each pwm ################    
        evalues = filter_using_evalue!(g; cpu=true, non_overlap=true, get_evalues_only=true);            
        sort_perm = sortperm(evalues); # sort it according to e-values (small to big)
        sort_perm_map = Dict(item=>index for (index,item) in enumerate(sort_perm));
        evalues = get_rounded_eval.(evalues[sort_perm]);
        # round.(round.(evalues[sort_perm] .* 10, sigdigits=3) ./ 10, sigdigits=3);
        # println(evalues[1])
        scores = [round(sum(values(g.ms.scores[i])),digits=1)
                    for i = 1:g.ms.num_motifs][sort_perm]; # scores of each discovered motifs
        msa_counts = [length(g.ms.positions[j]) for j = 1:g.ms.num_motifs][sort_perm];

        configurations, gap_configurations = get_configurations(g);
        selected_keys_sorted_update, weights_sorted, count_sorted, max_gap_num = make_configurations_plots_2(configurations, gap_configurations, pics_folder, sort_perm_map);
        weights_sorted = round.(weights_sorted, digits=3);

        labels = ["D$j" for j = 1:g.ms.num_motifs];
        logos = ["d$(j)" for j = 1:g.ms.num_motifs];
            
        save_pfms_as_transfac(logo_folder, g, sort_perm);

        # render the main result page
        df = DataFrame(label=labels, count=msa_counts, slrs=scores, eval=evalues, logo=logos);
        out = Mustache.render(html_template, target_folder=target_folder, logo_folder=logo_folder_name, num_seq=g.data.N, DF=df);
        ###### html stuff #########################            
        io = open(target_folder*"/summary.html", "w")
        print(io, out);
        close(io)

        # render the binding pattern page
        tpl=create_template_bindind_patterns(max_gap_num)
        out = Mustache.render(tpl, DF=get_df_for_binding_pattern(selected_keys_sorted_update, 
                                                                 logo_folder_name, 
                                                                 pics_folder_name,
                                                                 weights_sorted,
                                                                 count_sorted))
        io = open(target_folder*"/bp.html", "w")
        print(io, out);
        close(io)

        # for summary #############
        top3 = min(length(evalues),3);
        # compare with jaspar logo
        pfm_jaspar = get_pfm_from_transfac(source_folder_logo_transfac; float_type=eltype(g.ms.pfms[1]));
        complement_or_not = Bool[];
        for pfm in g.ms.pfms[sort_perm][1:top3]
            allr_score,_,_,_,_ = convolve_allr(pfm_jaspar, pfm, 1000, 1000, size(pfm_jaspar,2), size(pfm,2), 0f0, 3, g.search.bg);
            allr_score_c,_,_,_,_ = convolve_allr(pfm_jaspar, pfm_complement(pfm), 1000, 1000, size(pfm_jaspar,2), size(pfm,2), 0f0, 3, g.search.bg);
            if allr_score > allr_score_c
                push!(complement_or_not, false)
            else
                push!(complement_or_not, true)
            end            
        end        

        top3_evalues = evalues[1:top3];
        top3_logo_link = [complement_or_not[i] ? ref_save_where*"/logos/d$(i)_c.png" : ref_save_where*"/logos/d$i.png" for i = 1:top3];
        details_link = ref_save_where*"/summary.html";
    end
    
    return (name=matrix_name,
            jaspar_link="https://jaspar.genereg.net/matrix/"*matrix_name, 
            jaspar_logo=ref_logo_where, 
            top3_evalues=top3_evalues, 
            top3_logo_link=top3_logo_link, 
            details_link=details_link)
end
