# CCMAnalysis.jl

Convergent Cross Mapping (CCM) analysis package for causal inference in time series data.

## Installation

```julia
# Install directly from GitHub
using Pkg
Pkg.add(url="https://github.com/chaseU2/CCM_analysis_Julia.git")

# Activate the package
using CCMAnalysis
```

## How to call the funktion:

```julia
# Run CCM analysis with comprehensive parameters
results = run_ccm_analysis(
    "path/to/your/data.tsv",  # Tab/CSV file (columns = variables, rows = timepoints)
    L = 11o,                  # Library size (time series length)
    E = 2,                    # Embedding dimension
    tau = 1,                  # Time delay steps
    THRESHOLD = 0.8,          # Minimum ρ for significant causation
    save_plots = true,        # Saves convergence plots as PNG
    save_protocol = true,     # Generates analysis log
    output_dir = "ccm_results", # Custom output folder
    show_plots = true        # Set true for interactive plots
)
