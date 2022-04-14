"""
This is an attempt to create breeding simulation
"""

cd("/home/bensiv/Documents/GitHub/BreedingSimulation")

using Pkg
Pkg.activate(".")

using DataFrames, DataFramesMeta, Pipe

# Generating breeding population
number_of_markers = 1000
number_of_F2s = 200
# marker names
markers = string.(fill("Ma",number_of_markers), 1:number_of_markers)
genotypes = [-1, 0, 0, 1]

# generate pair of parents with homeozygote markers (polimorphic from each other)
parent1 = rand([-1,1], number_of_markers)
parent2 = parent1 .*-1


population_df = DataFrame(Markers = markers, Parent1 = parent1, Parent2 = parent2)
population_df.F1 = fill(0, nrow(population_df))

# generate a recombination frequency map for the breeding population
using Distributions, StatsBase
using Plots, StatsPlots

# 10 hotspots at the edges of the chromosome
ChromosomeLength = Int(1e6) # bp
Edge1 = (ChromosomeLength/1600):ChromosomeLength/16
Edge2 = reverse(ChromosomeLength .- Edge1)

Edge1HotSpots = rand(Edge1, 5)
Edge2HotSpots = rand(Edge2, 5)
HotSpots = [Edge1HotSpots ; Edge2HotSpots]

HotSpotsModel = Distributions.MixtureModel(Normal.(HotSpots, fill(200, length(HotSpots))),)
bp = 1:ChromosomeLength

#Recombination likelihood model
RecLikelihood = pdf.(HotSpotsModel, bp) .* 1e5

Plots.plot(bp, RecLikelihood, legend=nothing)
Plots.ylabel!("\$f_X(x)\$")
Plots.xlabel!("\$x\$")
Plots.title!("Gaussian mixture PDF")

SampleFrom = vcat(fill.(bp, Int.(round.(RecLikelihood)))...)

# A table with the recombination locations for each F2 plant
F2RecLocs = DataFrame(F2ID = Int[], RecLocs = Array{Int}[])
nF2s = 100
for plant in 1:nF2s
    # number of recombination occurrences
    nRec = rand(0:20)
    push!(F2RecLocs, [plant, sort(rand(SampleFrom, nRec))])
end

# generate random marker locations with uniform distribution across the chromosome
MarkerPos = sort(rand(bp, number_of_markers))
density(MarkerPos, xlimits = (1,ChromosomeLength))

MarkerPositions = DataFrame(MarkerID = markers, Position = MarkerPos)

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


heatmap(Matrix(population_df[!,5:end]))

