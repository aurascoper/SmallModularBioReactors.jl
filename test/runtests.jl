using Test
using LinearAlgebra
if !isdefined(Main, :SmallModularBioReactors)
    try
        @eval using SmallModularBioReactors
    catch
        include(joinpath(@__DIR__, "..", "src", "SmallModularBioReactors.jl"))
        @eval using .SmallModularBioReactors
    end
end

@testset "SmallModularBioReactors.jl Complete Verification Suite" begin
    include(joinpath(@__DIR__, "contracts", "test_contracts.jl"))
    include(joinpath(@__DIR__, "transport", "test_transport.jl"))
    include(joinpath(@__DIR__, "biological", "test_biological.jl"))
    include(joinpath(@__DIR__, "inference", "test_inference.jl"))
    include(joinpath(@__DIR__, "integration", "test_integration.jl"))
end

