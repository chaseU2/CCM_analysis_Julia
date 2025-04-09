function run_ccm_analysis(
    input_data;
    L=110,
    E=2,
    tau=1,
    THRESHOLD=0.8,
    save_plots=false,
    save_protocol=false,
    output_dir=nothing,
    show_plots=true
)
    """
    Run Convergent Cross Mapping (CCM) analysis on time series data.
    
    Parameters:
    -----------
    input_data : Union{String, DataFrame}
        Input data as file path or DataFrame
    L : Int
        Time series length (default: 110)
    E : Int
        Embedding dimension (default: 2)
    tau : Int
        Time delay (default: 1)
    THRESHOLD : Float64
        Threshold for significant interactions (default: 0.8)
    save_plots : Bool
        Whether to save convergence plots (default: false)
    save_protocol : Bool
        Whether to save protocol file (default: false)
    output_dir : Union{Nothing, String}
        Directory to save results (default: nothing)
    show_plots : Bool
        Whether to display plots interactively (default: true)
    
    Returns:
    --------
    DataFrame
        DataFrame with significant interactions
    """
    
    # Create output directory if needed
    if (save_plots || save_protocol) && output_dir === nothing
        output_dir = "ccm_results_" * Dates.format(now(), "yyyymmdd_HHMMSS")
    end
    
    if (save_plots || save_protocol) && !isnothing(output_dir)
        mkpath(output_dir)
    end
    
    # Helper function to create shadow manifold
    function shadow_manifold(time_series_Y, L, E, tau)
        shadow_M = Dict{Int, Vector{Float64}}()
        max_t = min(L-1, length(time_series_Y)-1)
        for t in (E - 1) * tau:max_t
            lag = Float64[]
            for t2 in 0:E-1
                idx = t - t2 * tau + 1
                if idx >= 1 && idx <= length(time_series_Y)
                    push!(lag, time_series_Y[idx])
                else
                    push!(lag, NaN)
                end
            end
            shadow_M[t] = lag
        end
        return shadow_M
    end

    # Helper function to calculate distance matrix
    function vec_dist_matrx(shadow_M)
        steps = Int[]
        vecs = Vector{Float64}[]
        for (k, v) in shadow_M
            push!(steps, k)
            push!(vecs, v)
        end
        distance_metrics = pairwise(Euclidean(), hcat(vecs...), dims=2)
        return distance_metrics, steps
    end

    # Helper function to find nearest neighbors (more robust version)
    function nearest_dist_and_step(timepoint_oi, steps, dist_matr, E)
        index_timepoint = findfirst(==(timepoint_oi), steps)
        if index_timepoint === nothing
            return Int[], Float64[]
        end
        
        dist_timepoint = @view dist_matr[index_timepoint, :]
        
        # Find all valid indices (excluding self and NaN distances)
        valid_indices = [i for i in 1:length(dist_timepoint) if i != index_timepoint && !isnan(dist_timepoint[i])]
        
        if length(valid_indices) < E
            return Int[], Float64[]
        end
        
        # Get indices of E nearest neighbors
        partial_sort_order = partialsortperm(dist_timepoint[valid_indices], 1:E, rev=false)
        nearest_indis = valid_indices[partial_sort_order]
        
        nearest_timesteps = steps[nearest_indis]
        nearest_distances = dist_timepoint[nearest_indis]
        
        return nearest_timesteps, nearest_distances
    end

    # Helper function for prediction
    function prediction(timepoint_oi, time_series_X, shadow_m, E)
        dist_matrix, steps = vec_dist_matrx(shadow_m)
        non_zero = 0.000001
        nearest_timesteps, nearest_distances = nearest_dist_and_step(timepoint_oi, steps, dist_matrix, E)
        
        if isempty(nearest_timesteps)
            return NaN, NaN
        end
        
        u = exp.(-nearest_distances ./ max(non_zero, minimum(nearest_distances)))
        w = u ./ sum(u)
        X_true = time_series_X[timepoint_oi + 1]
        X_cor = time_series_X[nearest_timesteps .+ 1]
        X_hat = sum(w .* X_cor)
        return X_true, X_hat
    end

    # Main causality calculation function
    function find_causality(time_series_X, time_series_Y, L, E, tau)
        My = shadow_manifold(time_series_Y, L, E, tau)
        X_true_list = Float64[]
        X_hat_list = Float64[]
        
        if all(==(time_series_X[1]), time_series_X) || all(==(time_series_Y[1]), time_series_Y)
            return 0.0, 1.0
        end
            
        for t in keys(My)
            X_true, X_hat = prediction(t, time_series_X, My, E)
            if !isnan(X_true) && !isnan(X_hat)
                push!(X_true_list, X_true)
                push!(X_hat_list, X_hat)
            end
        end
        
        if isempty(X_true_list) || isempty(X_hat_list)
            return 0.0, 1.0
        end
        
        r = cor(X_true_list, X_hat_list)
        p = isnan(r) ? 1.0 : 0.0
        return isnan(r) ? 0.0 : r, p
    end

    # Load and prepare data
    if input_data isa String
        ccm_data_matrix = readdlm(input_data, '\t', header=true)
        header = string.(ccm_data_matrix[2][1, :])
        data = ccm_data_matrix[1][2:end, :]
        ccm_data = DataFrame(data, header)
        rename!(ccm_data, Symbol.(header))
    elseif input_data isa DataFrame
        ccm_data = deepcopy(input_data)
    else
        error("input_data must be either a file path (String) or a DataFrame")
    end
    
    actual_L = min(L, size(ccm_data, 1))
    
    # Initialize final results container
    final_matrix = DataFrame(Species = String[])
    
    # Create protocol file if requested
    if save_protocol
        protocol_path = joinpath(output_dir, "protocol.txt")
        open(protocol_path, "w") do protocol
            write(protocol, "CCM Analysis Protocol - $(now())\n")
            write(protocol, "Parameters used:\n")
            write(protocol, "- L (time series length): $actual_L\n")
            write(protocol, "- E (embedding dimension): $E\n")
            write(protocol, "- tau (time delay): $tau\n")
            write(protocol, "- THRESHOLD: $THRESHOLD\n\n")
            
            # Initialize results dataframe
            species_names = names(ccm_data)
            results_df = DataFrame(zeros(length(species_names), length(species_names)), :auto)
            rename!(results_df, species_names)
            results_df[!, :Species] = species_names
            
            col_order = [:Species; Symbol.(species_names)]
            results_df = results_df[:, col_order]
            
            # Calculate initial CCM matrix
            write(protocol, "Calculating initial CCM matrix...\n")
            for (i, species1) in enumerate(species_names), (j, species2) in enumerate(species_names)
                r, p = find_causality(ccm_data[!, species1], ccm_data[!, species2], actual_L, E, tau)
                results_df[i, j+1] = max(0, round(r, digits=4))
            end
            
            # Convergence analysis
            L_range = 5:5:actual_L-1
            convergence_results = DataFrame(Species_X=String[], Species_Y=String[], X_to_Y_ρ=Float64[], Y_to_X_ρ=Float64[])

            analyzed_pairs = Set{Tuple{String, String}}()
            
            write(protocol, "Starting convergence analysis...\n")
            for (i, species1) in enumerate(species_names), (j, species2) in enumerate(species_names)
                if i == j || (species1, species2) in analyzed_pairs
                    continue
                end
                        
                score = results_df[i, j+1]
                if score <= THRESHOLD
                    continue
                end
                        
                push!(analyzed_pairs, (species1, species2))
                push!(analyzed_pairs, (species2, species1))
                    
                x_to_y = Float64[]
                y_to_x = Float64[]
                for L_val in L_range
                    r_xy, _ = find_causality(ccm_data[!, species1], ccm_data[!, species2], L_val, E, tau)
                    r_yx, _ = find_causality(ccm_data[!, species2], ccm_data[!, species1], L_val, E, tau)
                    push!(x_to_y, isnan(r_xy) ? 0.0 : r_xy)
                    push!(y_to_x, isnan(r_yx) ? 0.0 : r_yx)
                end
                    
                plt = plot(L_range, x_to_y, label="$species1 → $species2", color=:blue, marker=:circle)
                plot!(plt, L_range, y_to_x, label="$species2 → $species1", color=:red, marker=:square)
                xlabel!(plt, "Library Size (L)")
                ylabel!(plt, "Cross-Map Skill (ρ)")
                title!(plt, "CCM Convergence: $species1 ↔ $species2\nScore = $(round(score, digits=2))")
                plt[:grid] = true
                
                if show_plots
                    display(plt)
                end
                
                println("\nAnalyzing: $species1 ↔ $species2")
                println("Final ρ: $species1→$species2: $(round(x_to_y[end], digits=2)), $species2→$species1: $(round(y_to_x[end], digits=2))")
                    
                # Get user input
                decision = "0"
                while true
                    print("Convergence in? (1=both, 2=$species1→$species2 only, 3=$species2→$species1 only, 0=none): ")
                    decision = readline()
                    if decision in ["0", "1", "2", "3"]
                        break
                    end
                    println("Invalid input! Please enter 0-3")
                end
                
                if save_plots
                    plot_path = joinpath(output_dir, "convergence_$(species1)_vs_$(species2).png")
                    savefig(plt, plot_path)
                end
                    
                if decision == "1"
                    push!(convergence_results, (species1, species2, x_to_y[end], y_to_x[end]))
                elseif decision == "2"
                    push!(convergence_results, (species1, species2, x_to_y[end], 0.0))
                elseif decision == "3"
                    push!(convergence_results, (species1, species2, 0.0, y_to_x[end]))
                end
                
                write(protocol, "\nPair: $species1 <-> $species2\n")
                write(protocol, "ρ values: $(round(x_to_y[end], digits=2)), $(round(y_to_x[end], digits=2))\n")
                write(protocol, "Decision: $decision\n")
            end
            
            # Create final results without heatmap
            if !isempty(convergence_results)
                significant_species = unique(vcat(convergence_results.Species_X, convergence_results.Species_Y))
                sort!(significant_species)
                
                final_matrix = DataFrame(zeros(length(significant_species), length(significant_species)), :auto)
                rename!(final_matrix, significant_species)
                final_matrix[!, :Species] = significant_species
                
                col_order = [:Species; Symbol.(significant_species)]
                final_matrix = final_matrix[:, col_order]
                
                for row in eachrow(convergence_results)
                    if row.X_to_Y_ρ > THRESHOLD
                        col_idx = findfirst(==(row.Species_Y), names(final_matrix))
                        row_idx = findfirst(==(row.Species_X), final_matrix.Species)
                        final_matrix[row_idx, col_idx] = row.X_to_Y_ρ
                    end
                    if row.Y_to_X_ρ > THRESHOLD
                        col_idx = findfirst(==(row.Species_X), names(final_matrix))
                        row_idx = findfirst(==(row.Species_Y), final_matrix.Species)
                        final_matrix[row_idx, col_idx] = row.Y_to_X_ρ
                    end
                end
            else
                write(protocol, "\nNo significant results above threshold found.\n")
                final_matrix = DataFrame()
            end
        end
    else
        # Run analysis without protocol saving
        species_names = names(ccm_data)
        results_df = DataFrame(zeros(length(species_names), length(species_names)), :auto)
        rename!(results_df, species_names)
        results_df[!, :Species] = species_names
        
        col_order = [:Species; Symbol.(species_names)]
        results_df = results_df[:, col_order]
        
        # Calculate initial CCM matrix
        for (i, species1) in enumerate(species_names), (j, species2) in enumerate(species_names)
            r, p = find_causality(ccm_data[!, species1], ccm_data[!, species2], actual_L, E, tau)
            results_df[i, j+1] = max(0, round(r, digits=4))
        end
        
        # Convergence analysis
        L_range = 5:5:actual_L-1
        convergence_results = DataFrame(Species_X=String[], Species_Y=String[], X_to_Y_ρ=Float64[], Y_to_X_ρ=Float64[])

        analyzed_pairs = Set{Tuple{String, String}}()
        
        for (i, species1) in enumerate(species_names), (j, species2) in enumerate(species_names)
            if i == j || (species1, species2) in analyzed_pairs
                continue
            end
                    
            score = results_df[i, j+1]
            if score <= THRESHOLD
                continue
            end
                    
            push!(analyzed_pairs, (species1, species2))
            push!(analyzed_pairs, (species2, species1))
                
            x_to_y = Float64[]
            y_to_x = Float64[]
            for L_val in L_range
                r_xy, _ = find_causality(ccm_data[!, species1], ccm_data[!, species2], L_val, E, tau)
                r_yx, _ = find_causality(ccm_data[!, species2], ccm_data[!, species1], L_val, E, tau)
                push!(x_to_y, isnan(r_xy) ? 0.0 : r_xy)
                push!(y_to_x, isnan(r_yx) ? 0.0 : r_yx)
            end
                
            plt = plot(L_range, x_to_y, label="$species1 → $species2", color=:blue, marker=:circle)
            plot!(plt, L_range, y_to_x, label="$species2 → $species1", color=:red, marker=:square)
            xlabel!(plt, "Library Size (L)")
            ylabel!(plt, "Cross-Map Skill (ρ)")
            title!(plt, "CCM Convergence: $species1 ↔ $species2\nScore = $(round(score, digits=2))")
            plt[:grid] = true
            
            if show_plots
                display(plt)
            end
                
            println("\nAnalyzing: $species1 ↔ $species2")
            println("Final ρ: $species1→$species2: $(round(x_to_y[end], digits=2)), $species2→$species1: $(round(y_to_x[end], digits=2))")
                
            # Get user input
            decision = "0"
            while true
                print("Convergence in? (1=both, 2=$species1→$species2 only, 3=$species2→$species1 only, 0=none): ")
                decision = readline()
                if decision in ["0", "1", "2", "3"]
                    break
                end
                println("Invalid input! Please enter 0-3")
            end
            
            if save_plots
                plot_path = joinpath(output_dir, "convergence_$(species1)_vs_$(species2).png")
                savefig(plt, plot_path)
            end
                
            if decision == "1"
                push!(convergence_results, (species1, species2, x_to_y[end], y_to_x[end]))
            elseif decision == "2"
                push!(convergence_results, (species1, species2, x_to_y[end], 0.0))
            elseif decision == "3"
                push!(convergence_results, (species1, species2, 0.0, y_to_x[end]))
            end
        end
        
        # Create final results without heatmap
        if !isempty(convergence_results)
            significant_species = unique(vcat(convergence_results.Species_X, convergence_results.Species_Y))
            sort!(significant_species)
            
            final_matrix = DataFrame(zeros(length(significant_species), length(significant_species)), :auto)
            rename!(final_matrix, significant_species)
            final_matrix[!, :Species] = significant_species
            
            col_order = [:Species; Symbol.(significant_species)]
            final_matrix = final_matrix[:, col_order]
            
            for row in eachrow(convergence_results)
                if row.X_to_Y_ρ > THRESHOLD
                    col_idx = findfirst(==(row.Species_Y), names(final_matrix))
                    row_idx = findfirst(==(row.Species_X), final_matrix.Species)
                    final_matrix[row_idx, col_idx] = row.X_to_Y_ρ
                end
                if row.Y_to_X_ρ > THRESHOLD
                    col_idx = findfirst(==(row.Species_X), names(final_matrix))
                    row_idx = findfirst(==(row.Species_Y), final_matrix.Species)
                    final_matrix[row_idx, col_idx] = row.Y_to_X_ρ
                end
            end
        else
            final_matrix = DataFrame()
        end
    end
    
    if (save_plots || save_protocol) && !isnothing(output_dir)
        println("\nAnalysis complete. Results saved in: $output_dir")
    else
        println("\nAnalysis complete.")
    end
    
    return final_matrix
end