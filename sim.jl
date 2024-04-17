using PlotlyJS
using FileIO
using Images
using ImageView
using ArchGDAL
using StatsBase
using DelimitedFiles
using TOML
include("gridProb.jl")
using .gridProb
include("loam.jl")
#using .loam

function loadData(file_name)

    file = ArchGDAL.read(file_name)

    elv_data = ArchGDAL.read(file)

    elv_data[elv_data .<= 0] .= 0

    elv_data = reshape(elv_data,size(elv_data,1),size(elv_data,2))

    return(elv_data)
end

config_file = open("config.toml", "r")

config = TOML.parse(read(config_file, String))

close(config_file)


#Section 1: Loam parameters

elv_data = loadData(joinpath(@__DIR__,config["sim"]["file"]))#"island.tif"))

X = config["sim"]["centre_of_square_x"]
Y = config["sim"]["centre_of_square_y"]

#Try loading res data
try
    res = floor(Int,config["sim"]["square_size"]/2)
    elv_data = elv_data[Y-res:Y+res,X-res:X+res]
catch e
    @warn "res unspecified, defaulting to entire TIF image.\
    If image is large this may result in very slow simulations\
    or the inability to display the simulation."
end

#Load deposit information (if it exists)
try
    deposit_x = config["deposit"]["center_x"]
    deposit_y = config["deposit"]["center_y"]
    length = config["deposit"]["length"]
    trend = config["deposit"]["trend"]
    dip = config["deposit"]["dip"]

    type = "-"
    try
        type = config["deposit"]["type"]
    catch
        @warn "Type of deposit not decleared, defaulting to normal"
    end
    
    Reef = reef(size(elv_data,2)+deposit_x,size(elv_data,1)+deposit_y,50,250,70,elv_data,"-")

catch e
    print("Deposit configuration not found (or configuration incomplete, defaulting to singe sample)")


end

Reef = reef(res-20,res-30,50,250,70,elv_data,"-")

heat_map = revLoam(Reef,elv_data,arc_probs,100000,false)

max = maximum(heat_map)

for sample in Reef
        heat_map[sample.y,sample.x] = max+1 # Set the reef to the maximum value as to make it evident
end

z = elv_data
x_bounds, y_bounds = size(z)
x = 1:x_bounds
y = 1:y_bounds

heat_map = log.(heat_map .+ 1) # Display log-probablity

fig=Plot([surface(x=x, y=y, z=z, surfacecolor=heat_map, 
        colorscale=colors.gist_earth, colorbar_thickness=25, colorbar_len=0.75
        ), surface(x=x, y=y, z=100*ones(size(z)), surfacecolor=heat_map,
                        colorscale=colors.gist_earth,
                        showscale=false)],
Layout(width=600, height=600, font_size=10, scene_camera_eye=attr(x=1.6, y=1.6, z=1)))

fig # Display figure