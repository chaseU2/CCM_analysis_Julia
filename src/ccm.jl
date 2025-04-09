using DataFrames
using LinearAlgebra
using Statistics
using Dates
using DelimitedFiles
using Plots
using Distances

function robust_prompt(prompt_msg::String, options::Vector{String}=String[]; default::String="")
    """
    Robust input handling that works in both REPL and script execution
    """
    # Clear input buffer (compatible with all Julia versions)
    if isinteractive()
        ccall(:jl_process_events, Nothing, ())
        while bytesavailable(stdin) > 0
            read(stdin, Char)
        end
    end
    
    while true
        try
            print(prompt_msg)
            response = strip(readline())
            
            # Use default if empty response provided
            if isempty(response) && !isempty(default)
                return default
            end
            
            # Validate against options if provided
            if isempty(options) || response ∈ options
                return response
            else
                println("Invalid input! Valid options: ", join(options, ", "))
            end
        catch e
            println("Input error - please try again")
        end
    end
end

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
    # Create output directory if needed
    if (save_plots || save_protocol) && output_dir === nothing
        output_dir = "ccm_results_" * Dates.format(now(), "yyyymmdd_HHMMSS")
    end
    
    if (save_plots || save_protocol) && !isnothing(output_dir)
        mkpath(output_dir)
    end

    # Shadow manifold construction
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

    # Distance matrix calculation
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

    # Nearest neighbors finding
    function nearest_dist_and_step(timepoint_oi, steps, dist_matr, E)
        index_timepoint = findfirst(==(timepoint_oi), steps)
        if index_timepoint === nothing
            return Int[], Float64[]
        end
        
        dist_timepoint = @view dist_matr[index_timepoint, :]
        valid_indices = [i for i in 1:length(dist_timepoint) if i != index_timepoint && !isnan(dist_timepoint[i])]
        
        if length(valid_indices) < E
            return Int[], Float64[]
        end
        
        partial_sort_order = partialsortperm(dist_timepoint[valid_indices], 1:E, rev=false)
        nearest_indis = valid_indices[partial_sort_order]
        
        nearest_timesteps = steps[nearest_indis]
        nearest_distances = dist_timepoint[nearest_indis]
        
        return nearest_timesteps, nearest_distances
    end

    # Prediction function
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

    # Causality calculation
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

    # Data loading
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
    species_names = names(ccm_data)
    
    # Initialize results
    results_df = DataFrame(zeros(length(species_names), length(species_names)), :auto)
    rename!(results_df, species_names)
    results_df[!, :Species] = species_names
    results_df = results_df[:, [:Species; Symbol.(species_names)]]

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
        i == j && continue
        (species1, species2) in analyzed_pairs && continue
        
        score = results_df[i, j+1]
        score <= THRESHOLD && continue

        push!(analyzed_pairs, (species1, species2), (species2, species1))
        
        # Calculate convergence
        x_to_y = [find_causality(ccm_data[!, species1], ccm_data[!, species2], L_val, E, tau)[1] for L_val in L_range]
        y_to_x = [find_causality(ccm_data[!, species2], ccm_data[!, species1], L_val, E, tau)[1] for L_val in L_range]
        
        # Plotting
        plt = plot(L_range, x_to_y, label="$species1 → $species2", color=:blue, marker=:circle)
        plot!(plt, L_range, y_to_x, label="$species2 → $species1", color=:red, marker=:square)
        xlabel!(plt, "Library Size (L)")
        ylabel!(plt, "Cross-Map Skill (ρ)")
        title!(plt, "CCM Convergence: $species1 ↔ $species2\nScore = $(round(score, digits=2))")
        plt[:grid] = true
        
        show_plots && display(plt)
        save_plots && !isnothing(output_dir) && savefig(plt, joinpath(output_dir, "convergence_$(species1)_vs_$(species2).png"))

        # User decision with robust input
        decision = robust_prompt("""
        Convergence for $species1 ↔ $species2?
        1 = both directions
        2 = $species1→$species2 only
        3 = $species2→$species1 only
        0 = none
        Enter choice (0-3):""", ["0", "1", "2", "3"])

        # Store results based on decision
        if decision == "1"
            push!(convergence_results, (species1, species2, x_to_y[end], y_to_x[end]))
        elseif decision == "2"
            push!(convergence_results, (species1, species2, x_to_y[end], 0.0))
        elseif decision == "3"
            push!(convergence_results, (species1, species2, 0.0, y_to_x[end]))
        end
    end

    # Generate final results
    final_matrix = if !isempty(convergence_results)
        sig_species = unique(vcat(convergence_results.Species_X, convergence_results.Species_Y))
        sort!(sig_species)
        
        df = DataFrame(zeros(length(sig_species), length(sig_species)), :auto)
        rename!(df, sig_species)
        df[!, :Species] = sig_species
        df = df[:, [:Species; Symbol.(sig_species)]]
        
        for row in eachrow(convergence_results)
            if row.X_to_Y_ρ > THRESHOLD
                df[findfirst(==(row.Species_X), df.Species), findfirst(==(row.Species_Y), names(df))] = row.X_to_Y_ρ
            end
            if row.Y_to_X_ρ > THRESHOLD
                df[findfirst(==(row.Species_Y), df.Species), findfirst(==(row.Species_X), names(df))] = row.Y_to_X_ρ
            end
        end
        df
    else
        DataFrame()
    end

    # Save protocol if requested
    if save_protocol && !isnothing(output_dir)
        protocol_path = joinpath(output_dir, "protocol.txt")
        open(protocol_path, "w") do io
            write(io, "CCM Analysis Protocol - $(now())\n\n")
            write(io, "Parameters:\n")
            write(io, "- L: $actual_L\n- E: $E\n- tau: $tau\n- Threshold: $THRESHOLD\n\n")
            write(io, "Significant Interactions:\n")
            
            if !isempty(convergence_results)
                for row in eachrow(convergence_results)
                    write(io, "$(row.Species_X) ↔ $(row.Species_Y): ")
                    write(io, "X→Y ρ=$(round(row.X_to_Y_ρ, digits=3)), ")
                    write(io, "Y→X ρ=$(round(row.Y_to_X_ρ, digits=3))\n")
                end
            else
                write(io, "No significant interactions found above threshold\n")
            end
        end
    end

    return final_matrix
end