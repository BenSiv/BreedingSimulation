"""
Building Genetic map based on recombination observations in breeding population
"""


cd("/home/bensiv/Documents/GitHub/BreedingSimulation")

using Pkg
Pkg.activate(".")

using CSV, DataFrames, DataFramesMeta, Plots

population_df = CSV.read("PopulationGenotypes.csv", DataFrame)

F2mat = Matrix(population_df[!,5:end])
heatmap(F2mat)

Parentmat = Matrix(population_df[!,2:4])
Parentmat = Parentmat[:,[1,3,2]]

F2Identity = hcat([argmax.(eachrow(F2mat[:,col] .== Parentmat)) for col in 1:size(F2mat)[2]]...)
heatmap(F2Identity)

savefig("F2Identity.png")
    
SumCols = sum.(eachcol(F2Identity))
F2Identity_sorted = F2Identity[:,sortperm(SumCols)]
heatmap(F2Identity_sorted)

savefig("F2Identity_sorted.png")

using DataStructures
CentiMorgans_diff = []
for row in 2:size(F2Identity_sorted)[1]
    RecVec = F2Identity_sorted[row,:] .== F2Identity_sorted[row-1,:]
    RecCount = DataStructures.counter(RecVec)
    cM = (RecCount[0] / sum(values(RecCount))) * (size(F2Identity_sorted)[2] / size(F2Identity_sorted)[1]) * 100
    push!(CentiMorgans_diff, cM)
end

CentiMorgans = [0.0]
for cM in CentiMorgans_diff
    push!(CentiMorgans, cM + last(CentiMorgans))
end

MarkerPositions = CSV.read("MarkerPositions.csv", DataFrame)

scatter(MarkerPositions.Position, CentiMorgans, legend = false, grid = false, xlabel = "Phisical position (bp)", ylabel = "Genetic position (cM)")

savefig("GeneticMap.png")

GeneticMap = MarkerPositions
GeneticMap.CentiMorgans = CentiMorgans

CSV.write("GeneticMap.csv", GeneticMap)