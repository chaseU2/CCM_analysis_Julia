using Documenter, CCMAnalysis

makedocs(
    sitename = "CCMAnalysis",
    modules = [CCMAnalysis],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "https://github.com/chaseU2/CCM_analysis_Julia.git"
)