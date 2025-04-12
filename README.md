# CCMAnalysis.jl

Package implementing Convergent Cross Mapping for causality inference in dynamical systems as defined by ["Sugihara et al (2012)"](https://www.science.org/doi/10.1126/science.1227079).

## Installation

```julia
# Install directly from GitHub
using Pkg
Pkg.add(url="https://github.com/chaseU2/CCM_analysis_Julia.git")

# Activate the package
using CCM_analysis_Julia
```

## How to call the funktion:

```julia
# Run CCM analysis with comprehensive parameters
results = run_ccm_analysis(
    "path/to/your/data.tsv",  # Tab/CSV file (columns = variables, rows = timepoints)
    L = 110,                  # Library size (time series length)
    E = 2,                    # Embedding dimension
    tau = 1,                  # Time delay steps
    THRESHOLD = 0.8,          # Minimum ρ for significant causation
    save_plots = true,        # Saves convergence plots as PNG
    save_protocol = true,     # Generates analysis log
    output_dir = "ccm_results", # Custom output folder
    show_plots = true        # Set true for interactive plots
)
```


### User-Provided Data

If you want to use your own data for testing or analysis, please ensure that your data file is in the **same format** as [`ccm_test_data.txt`](https://github.com/chaseU2/ccm-analysis-tool/blob/master/ccm_test_data.txt). Importantly, **all missing values (NAs)** in the dataset **must be replaced with zeros** to ensure proper functionality of the algorithm.

You can upload your data file in the same format, with columns and rows matching the original data structure.




## Dependencies and Acknowledgements

Parts of this project are based on the Convergent Cross Mapping (CCM) implementation from the following repository from Prince Javier :

- [Convergent Cross Mapping GitHub Repository Prince Javier ](https://github.com/PrinceJavier/causal_ccm.git)

I have utilized parts of the CCM algorithm from this repository to help analyze causality in time series data in my own project.
