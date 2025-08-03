module ReproduceEntropyStableDSEM

using JLD2
using OrdinaryDiffEq
using TimerOutputs
using LaTeXStrings
using Suppressor
using StableSpectralElements
using LinearAlgebra

export run_driver

export EulerDriver
include("euler_refinement.jl")

export EulerPRefinementDriver
include("euler_p_refinement.jl")
end
