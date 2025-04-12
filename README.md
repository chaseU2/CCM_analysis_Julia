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



## Interpreting Convergence Plots

The package generates diagnostic convergence plots that reveal causal relationships between variables. Below are the three characteristic patterns to analyze:


When analyzing convergence plots, you'll be prompted:

```julia
Analyzing: A ↔ B
Final ρ: A→B: 0.82, B→A: 0.41

Convergence in? (1=both, 2=A→B only, 3=B→A only, 0=none): 
```
p is in this case the calculated crossmap score 

### ⚠️ Input Instructions
If the console doesn't respond to your first number:
1. Press **Enter** to submit
2. Type the number **again**
3. Press **Enter** again

This ensures proper input handling.

### Options Table
| Key | Action                  |
|-----|-------------------------|
| 1   | Keep both directions    |
| 2   | Keep only A→B           |
| 3   | Keep only B→A           |
| 0   | Discard both            |


### 1. No Significant Causality

![No Causal Relationship](https://raw.githubusercontent.com/chaseU2/CCM_analysis_Julia/main/src/Screenshot%205.png)

- Neither directional curve (X→Y in blue, Y→X in red) show a clear convergence to a final crossmap score
- Example use case: Independent systems

### 2. Unidirectional Causality
![Unidirectional Causality](https://raw.githubusercontent.com/chaseU2/CCM_analysis_Julia/main/src/Screenshot%204.png)

- One direction (X→Y) converges to a final crossmap score
- Reverse direction (Y→X) does ot show a clear convegence
- Interpretation: X drives Y but not vice versa
- Pay attention to whether convergence is present in the X→Y or Y→X plot, and enter 2 or 3 accordingly

### 3. Bidirectional Causality
![Bidirectional Causality](https://raw.githubusercontent.com/chaseU2/CCM_analysis_Julia/main/src/Screenshot%202.png)

- Both directions show positive convergence
- Typical of feedback systems
- Convergence rates may differ (e.g., X→Y stronger than Y→X)

---

## Generating and Customizing Plots

### Basic Plot Generation
```julia
results = run_ccm_analysis(
    your_dataframe,
    save_plots = true,       # Required to save images
    output_dir = "analysis",  # Custom folder name
    show_plots = true,       # Interactive display
    plot_theme = :default    # :dark, :ggplot2, etc.
)
```
---

## Dependencies and Acknowledgements

Parts of this project are based on the Convergent Cross Mapping (CCM) implementation from the following repository from Prince Javier :

- [Convergent Cross Mapping GitHub Repository Prince Javier ](https://github.com/PrinceJavier/causal_ccm.git)

I have utilized parts of the CCM algorithm from this repository to help analyze causality in time series data in my own project.
