using PlotlyJS
using FileIO
using Images
using ImageView
using ArchGDAL
using StatsBase
#using CairoMakie
#using Plots
using DelimitedFiles
include("gridProb.jl")
using .gridProb
#using GLMakie
arc_probs = gridProb.arc_probs

function check(x,i) # Debugging function to check nugget is always going down hill 
    if(x[i] < x[i+1])
        throw("You've still got a bug")
    else
        if(i+1 == length(x)) return 0 else return check(x,i+1) end
    end
end

#Section 1: Open TIF File

file = ArchGDAL.read(joinpath(@__DIR__,"christmasislandz48.tif"))

data = ArchGDAL.getband(file,1) 

elv_data = ArchGDAL.read(file)

elv_data[elv_data .<= 0] .= 0

elv_data = reshape(elv_data,4200,3600)

struct Sample
    x::Int64
    y::Int64
    type::String
end

function flatten(td_array,length,width)::Array{Float64,1} # Builtin flatten gives unsuitable results
    z = nothing
    for x in 1:length
        for y in 1:width
            if (z === nothing) z = [td_array[x,y]] else z = [z;td_array[x,y]] end
        end
    end
    return convert(Array{Float64,1}, z)
end

function choose(elv::Matrix{Float32},arc_probs::Matrix{Float32},mode::Bool)::Array{Int64, 1} # returns delta y and delta x

    current_elevation = elv[2,2]
    relative_elevation = Nothing
    #Mode 0 -> Loam, Mode 1 -> Rev-Loam 
    if(mode) relative_elevation = elv .- current_elevation else relative_elevation = -current_elevation .+ elv end

    if minimum(relative_elevation) == 0 # Everything is uphill
        return [0,0]
    else
        relative_elevation[relative_elevation .> 0] .= 0
        relative_elevation = abs.(relative_elevation)

        weights = (relative_elevation ./ sum(relative_elevation)).^1 # Linear weighting - possible improvement here
        weights = weights .* arc_probs # Multiple weights by arc length
        weights = weights ./ sum(weights)

        choice = sample([[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]], Weights(flatten(weights,3,3)))
    
        if(choice == [0,0]) throw("Something impossible was chosen") end 
        if(elv[2+choice[2],2+choice[1]] > current_elevation) throw("Something higher was chosen in rev-loam mode or something lower was chosen in loam mode") end

        return choice

    end
end


function trace(start_x::Int64,start_y::Int64,elv::Array{Float32, 2},arc_probs::Matrix{Float32},mode::Bool)::Matrix{Float32}

    x = start_x
    y = start_y
    choice = [1,1]
    steps = 0

    map = zeros(size(elv))
    
    route = Nothing

    while choice != [0,0]
        adjacent_heights = elv[y-1:y+1,x-1:x+1]
        choice = choose(adjacent_heights,arc_probs,mode)
        x += choice[1]
        y += choice[2]
        if(route === Nothing) route = [x y elv[y,x]] else route = [route;[x y elv[y,x]]] end
        steps += 1
        map[y,x] = 1
    end

    return map
end

function monte_carlo(x::Int64,y::Int64,elv::Array{Float32,2},arc_probs::Matrix{Float32},n_simulations::Int64,mode::Bool)::Matrix{Float64}
    heat_map = zeros(size(elv))
    for _ in 1:n_simulations
        Trace = trace(x,y,elv,arc_probs,mode)
        heat_map .+= Trace
    end
    heat_map
end

function loam(Samples::Vector{Sample},elv::Array{Float32,2},arc_probs::Matrix{Float32},simulations_per_sample::Int64)::Matrix{Float64}
    heat_map = zeros(size(elv))
    for sample in Samples
        heat_map += monte_carlo(sample.x,sample.y,elv,arc_probs,simulations_per_sample,false)
    end
    heat_map
end

function loam(sample::Sample,elv::Array{Float32,2},arc_probs::Matrix{Float32},n_simulations::Int64)::Matrix{Float64}

    nthreads = Threads.nthreads() 
    sims_per_sample::Int64 = convert(Int64, n_simulations/nthreads::Float64)

    ThreadSums = Vector{Matrix{Float64}}(undef, nthreads)

    Threads.@threads for i in 1:nthreads
        ThreadSums[i] = monte_carlo(sample.x,sample.y,elv,arc_probs,sims_per_sample,false)
    end

    total_heat_map = zeros(size(elv))

    for sum in ThreadSums total_heat_map += sum end

    return(total_heat_map)
end

#1587
#1108
X = 1587
Y = 1108

elv_data = elv_data[Y-200:Y+200,X-200:X+200]

heat_map = monte_carlo(201,201,elv_data,arc_probs,10000,true)

z = elv_data
heat_map = heat_map
r, c = size(z)
x = 1:c
y = 1:r

fig=Plot([surface(x=x, y=y, z=z, surfacecolor=heat_map, 
        colorscale=colors.gist_earth, colorbar_thickness=25, colorbar_len=0.75
        ), surface(x=x, y=y, z=-1*ones(size(z)), surfacecolor=heat_map,
                        colorscale=colors.gist_earth,
                        showscale=false)],
Layout(width=600, height=600, font_size=10, scene_camera_eye=attr(x=1.6, y=1.6, z=1)))

fig