using DmspElfinConjunction, DMSP
using Dates, TimeseriesUtilities
using DimensionalData, DataFrames, DataFramesMeta
using Accessors
using AlgebraOfGraphics
using CairoMakie
using StatsBase, CategoricalArrays
import Beforerr
using Beforerr: easy_save
using LaTeXStrings

Beforerr.DEFAULT_FORMATS = [:png, :pdf]

include(joinpath(@__DIR__, "../src/plot.jl"))
include(joinpath(@__DIR__, "../src/workload.jl"))
include(joinpath(@__DIR__, "../src/statplot.jl"))
include(joinpath(@__DIR__, "../src/analysis.jl"))
