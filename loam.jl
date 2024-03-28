module loam

export Sample
export Loam
export revLoam

struct Sample
    x::Int64
    y::Int64
    type::String
    weighting::Float64
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
    if(mode) relative_elevation = elv .- current_elevation else relative_elevation = -1*elv .+ current_elevation end

    if minimum(relative_elevation) == 0 # Everything is uphill
        return [0,0]
    else
        relative_elevation[relative_elevation .> 0] .= 0
        relative_elevation = abs.(relative_elevation)

        grid_weights = (relative_elevation ./ sum(relative_elevation))#.^0.5 # Linear weighting - possible improvement here
        grid_weights = grid_weights .* arc_probs # Multiple weights by arc length
        grid_weights = grid_weights ./ sum(grid_weights)

        choice = sample([[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]], Weights(flatten(grid_weights,3,3)))
    
        if(choice == [0,0]) throw("Something impossible was chosen") end

        if(mode)
            if(elv[2+choice[2],2+choice[1]] > current_elevation) throw("Something higher was chosen in rev-loam mode") end
        else
            if(elv[2+choice[2],2+choice[1]] < current_elevation) throw("Something lower was chosen in loam mode") end
        end 

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

function Loam(Samples::Vector{Sample},elv::Array{Float32,2},arc_probs::Matrix{Float32},simulations_per_sample::Int64)::Matrix{Float64}
    heat_map = zeros(size(elv))
    for sample in Samples
        heat_map += monte_carlo(sample.x,sample.y,elv,arc_probs,simulations_per_sample,false)
    end
    heat_map
end

function revLoam(Samples::Vector{Sample},elv::Array{Float32,2},arc_probs::Matrix{Float32},simulations_per_sample::Int64)::Matrix{Float64}
    heat_map = zeros(size(elv))
    for sample in Samples
        heat_map += monte_carlo(sample.x,sample.y,elv,arc_probs,simulations_per_sample,true)
    end
    heat_map
end

function Loam(sample::Sample,elv::Array{Float32,2},arc_probs::Matrix{Float32},n_simulations::Int64)::Matrix{Float64}

    nthreads = Threads.nthreads() 
    sims_per_sample::Int64 = convert(Int64, n_simulations/nthreads)

    ThreadSums = Vector{Matrix{Float64}}(undef, nthreads)

    Threads.@threads for i in 1:nthreads
        ThreadSums[i] = monte_carlo(sample.x,sample.y,elv,arc_probs,sims_per_sample,false)
    end

    total_heat_map = zeros(size(elv))

    for sum in ThreadSums total_heat_map += sum end

    return(total_heat_map)
end

function revLoam(sample::Sample,elv::Array{Float32,2},arc_probs::Matrix{Float32},n_simulations::Int64)::Matrix{Float64}

    nthreads = Threads.nthreads() 
    sims_per_sample::Int64 = convert(Int64, n_simulations/nthreads)

    ThreadSums = Vector{Matrix{Float64}}(undef, nthreads)

    Threads.@threads for i in 1:nthreads
        ThreadSums[i] = monte_carlo(sample.x,sample.y,elv,arc_probs,sims_per_sample,true)
    end

    total_heat_map = zeros(size(elv))

    for sum in ThreadSums total_heat_map += sum end

    return(total_heat_map)
end

end # End of module 