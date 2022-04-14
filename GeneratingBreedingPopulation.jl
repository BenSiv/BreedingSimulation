"""
This script generates a breeding population genotypes.
"""

cd("/home/bensiv/Documents/GitHub/BreedingSimulation")

using Pkg
Pkg.activate(".")

using YAML
PopParams = YAML.load_file("PopParams.yml")

using CSV, DataFrames, DataFramesMeta, Pipe

# Generating breeding population
nMarkers = PopParams["number_of_markers"]
nF2s = PopParams["number_of_F2s"]
# marker names
markers = string.(fill("Ma",nMarkers), 1:nMarkers)
genotypes = [-1, 0, 0, 1]

# generate pair of parents with homeozygote markers (polimorphic from each other)
parent1 = rand([-1,1], nMarkers)
parent2 = parent1 .*-1


population_df = DataFrame(Markers = markers, Parent1 = parent1, Parent2 = parent2)
population_df.F1 = fill(0, nrow(population_df))

# generate a recombination frequency map for the breeding population
using Distributions, StatsBase
using Plots, StatsPlots

# n hotspots at the edges of the chromosome
ChromosomeLength = Int(PopParams["chromosome_length"]) # bp
Edge1 = (ChromosomeLength*(PopParams["edge_portion"]/100)):(ChromosomeLength*PopParams["edge_portion"])
Edge2 = reverse(ChromosomeLength .- Edge1)

nHotSpots1 = Int(round(PopParams["number_of_hotspots"] / 2))
nHotSpots2 = Int(PopParams["number_of_hotspots"] - nHotSpots1)

Edge1HotSpots = rand(Edge1, nHotSpots1)
Edge2HotSpots = rand(Edge2, nHotSpots2)
HotSpots = [Edge1HotSpots ; Edge2HotSpots]

HotSpotsModel = Distributions.MixtureModel(Normal.(HotSpots, fill(ChromosomeLength*PopParams["hotspots_diffusion"], length(HotSpots))),)
bp = 1:ChromosomeLength

#Recombination likelihood model
RecLikelihood = pdf.(HotSpotsModel, bp) .* ChromosomeLength

RecLikly_plot = Plots.plot(bp, RecLikelihood, legend=nothing)
Plots.title!("Recombination Likelihood %")

savefig(RecLikly_plot, "Recombination_Likelihood.png")

SampleFrom = vcat(fill.(bp, Int.(round.(RecLikelihood)))...)

# A table with the recombination locations for each F2 plant
F2RecLocs = DataFrame(F2ID = Int[], RecLocs = Array{Int}[])
for plant in 1:nF2s
    # number of recombination occurrences
    nRec = rand(0:PopParams["max_recombination"])
    push!(F2RecLocs, [plant, sort(rand(SampleFrom, nRec))])
end

# generate random marker locations with uniform distribution across the chromosome
MarkerPos = sort(rand(bp, nMarkers))

MarkerPositions = DataFrame(MarkerID = markers, Position = MarkerPos)
CSV.write("MarkerPositions.csv", MarkerPositions)

# generate F2 plants
for plant in 1:nF2s
    ParentIdentities = rand([2,3,4,4], length(F2RecLocs.RecLocs[plant])+1)
    IdentityIndexs = [0 ; F2RecLocs.RecLocs[plant] ; ChromosomeLength]

    GenotypeCollection = []
    for (ind,col) in zip(1:(length(IdentityIndexs)-1), ParentIdentities)
        MarkersInRange = findall(in(IdentityIndexs[ind]+1:IdentityIndexs[ind+1]).(MarkerPositions.Position))
        push!(GenotypeCollection, population_df[MarkersInRange,col])
    end

    population_df[!, string(plant+4)] = vcat(GenotypeCollection...)
end

F2mat = Matrix(population_df[!,5:end])
heatmap(F2mat)
savefig("GenotypesHeatmap.png")


using CSV
CSV.write("PopulationGenotypes.csv", population_df)