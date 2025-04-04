module CCM_analysis_Julia

using DelimitedFiles
using DataFrames
using Statistics
using LinearAlgebra
using Distances
using Plots
using StatsBase
using Dates
using Printf
using Measures

include("ccm.jl")  # Hier kommt Ihr Hauptcode

export run_ccm_analysis  # Exportieren Sie die Hauptfunktion

end
