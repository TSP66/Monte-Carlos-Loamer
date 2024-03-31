using PlotlyJS
using FileIO
using Images
using ImageView
using ArchGDAL
using StatsBase
using DelimitedFiles
include("gridProb.jl")
using .gridProb
include("loam.jl")
#using .loam

#Section 1: Open TIF File

file = ArchGDAL.read(joinpath(@__DIR__,"christmasislandz48.tif"))

data = ArchGDAL.getband(file,1) 

elv_data = ArchGDAL.read(file)

elv_data[elv_data .<= 0] .= 0

elv_data = reshape(elv_data,4200,3600)

#1587
#1108
X = 1587
Y = 1108

res = 240

elv_data = elv_data[Y-res:Y+res,X-res:X+res]

#heat_map = monte_carlo(201,201,elv_data,arc_probs,30000,true)
#heat_map = revLoam(Sample(res+1,res+1,"-",0.0),elv_data,arc_probs,100000,false)
#heat_map = Loam([Sample(res+25,res+25,"-",1), Sample(res+30,res+50,"-",1)],elv_data,arc_probs,250000,true)

heat_map = revLoam(reef(res+20,res+20,100,125,"Normal"),elv_data,arc_probs,30000,false)

z = elv_data
r, c = size(z)
x = 1:c
y = 1:r

heat_map = log.(heat_map .+ 1)

fig=Plot([surface(x=x, y=y, z=z, surfacecolor=heat_map, 
        colorscale=colors.gist_earth, colorbar_thickness=25, colorbar_len=0.75
        ), surface(x=x, y=y, z=-1*ones(size(z)), surfacecolor=heat_map,
                        colorscale=colors.gist_earth,
                        showscale=false)],
Layout(width=600, height=600, font_size=10, scene_camera_eye=attr(x=1.6, y=1.6, z=1)))

fig