module deposits

export Sample
export vertical_reef
export reef

struct Sample
    x::Int64
    y::Int64
    type::String
    weight::Float64
end

function vertical_reef(centre_x::Int64,centre_y::Int64,trend::Int64,length::Int64,type::String) # centre_x,centre_y are the centres of the reef

    samples = Vector{Sample}([])
    positions = nothing

    for point in Int(round(-length/2)):Int(round(length/2))

        x = Int(round(centre_x + point * cos(trend*pi/180)))
        y = Int(round(centre_y + point * sin(trend*pi/180)))

        if (positions === nothing)
            samp = Sample(x,y,type,1)
            samples = [samp]
            positions = [[x,y]]
        else
            if([x,y] ∉ positions) # Avoid duplicates
                samp = Sample(x,y,type,1)
                push!(samples,samp)
                push!(positions,[x,y])
            
            end
        end
    end
    samples

end

function reef(centre_x::Int64,centre_y::Int64,trend::Int64,length::Int64,dip::Int64,elv::Array{Float32,2},type::String)

    if(dip == 90) return(vertical_reef(centre_x,centre_y,trend,length,type)) end # Way more efficient 

    samples = Vector{Sample}([])

    plane(x,y) = tan(dip*pi/180)*((x-centre_x)*sin(trend*pi/180)^2+(y-centre_y)*cos(trend*pi/180)^2)+elv[centre_y,centre_x]

    min_x = Int(round(centre_x + -length/2 * cos(trend*pi/180)))
    max_x = Int(round(centre_x + length/2 * cos(trend*pi/180)))
    min_y = Int(round(centre_y + -length/2 * sin(trend*pi/180)))
    max_y = Int(round(centre_y + length/2 * sin(trend*pi/180)))

    for y in min_y:max_y
        for x in min_x:max_x
            if(abs(plane(x,y)-elv[y,x])<1) 
                push!(samples,Sample(x,y,type,1.0)) 
            end
        end
    end

    samples

end

end