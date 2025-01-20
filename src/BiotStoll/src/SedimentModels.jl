# SedimentModels2.jl (Referencing AbstractSediments2 correctly)
module SedimentModels

export Sediment, create_sediment_with_fluid, predefined_sediment

# Define a single Sediment struct
struct Sediment
    name::String
    d::Float64
    ρr::Float64
    Kr::Float64
    k::Float64
    γ::Float64
    δl::Float64
    δs::Float64
    ρF::Float64
    KF::Float64
    η::Float64
    ϕ::Float64
    β::Float64
    a::Float64
    Kb::Float64
    μ::Float64
end

function create_sediment_with_fluid(name::String, d, ρr, Kr, k, γ, δl, δs; ρF=998.0, KF=2.195e9, η=1.0e-3)
    ϕ = -log2(d)
    β = 0.3105 + 0.0552 * ϕ
    a = d * 1e-3 / 3 * β / (1 - β)
    Kb = 10^(2.70932 - 4.25391 * β) * 1e8
    ρ = β * ρF + (1 - β) * ρr
    cs = 4.40*β^(-2.69)
    μ = ρ * cs^2
    return Sediment(name, d, ρr, Kr, k, γ, δl, δs, ρF, KF, η, ϕ, β, a, Kb, μ)
end

function predefined_sediment(name::String; ρF=998.0, KF=2.195e9, η=1.0e-3)
    if name == "MediumSilt"
        return create_sediment_with_fluid("MediumSilt", 0.0221, 2650, 3.6e10, 4.22e-12, 1.25, 0.15, 0.15; ρF=ρF, KF=KF, η=η)
    elseif name == "VeryFineSand"
        return create_sediment_with_fluid("VeryFineSand", 0.0884, 2650, 3.6e10, 2.26e-11, 1.25, 0.15, 0.15; ρF=ρF, KF=KF, η=η)
    elseif name == "MediumSand"
        return create_sediment_with_fluid("MediumSand", 0.354, 2650, 3.6e10, 1.15e-10, 1.25, 0.15, 0.15; ρF=ρF, KF=KF, η=η)
    else
        error("Sediment '$name' not recognized.")
    end
end

"""
    list_sediments()

List all available predefined sediments.
"""
function list_sediments()
    return ["MediumSilt", "VeryFineSand", "MediumSand"]
end

"""
    select_sediment(name::String)

Select a predefined sediment with default fluid properties.
"""
function select_sediment(name::String)::Union{Sediment, Nothing}
    try
        return predefined_sediment(name)
    catch e
        println("Error: ", e)
        return nothing
    end
end


end  # End of module
