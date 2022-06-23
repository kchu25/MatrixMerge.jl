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

