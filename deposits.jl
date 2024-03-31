module deposits

export Sample
export reef

struct Sample
    x::Int64
    y::Int64
    type::String
    weight::Float64
end

function reef(centre_x::Int64,centre_y::Int64,trend::Int64,length::Int64,type::String) # centre_x,centre_y are the centres of the reef

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

end