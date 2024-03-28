#Calcuate probabilities of nugget travelling into the 8 adjacent blocks

#Weight probability by arc-length of each section

module gridProb

export arc_probs

one = 0.9999999999999 # one = 1 generates errors as asin(1) is undefined

IOCF(x)::Float32 = sqrt(1/(1-x^2))*sqrt(1-x^2)*asin(x)

arc_length(x_1,x_2)::Float32 = IOCF(x_2) - IOCF(x_1)

arc_probs = Matrix{Float32}(undef, 3, 3)

arc_probs[2,2] = 0

arc_probs[2,1] = 2*arc_length(-one,-sqrt(8)/3) # Center left block
arc_probs[1,1] = arc_length(-sqrt(8)/3,-1/3) # Top left block
arc_probs[1,2] = arc_length(-1/3,1/3) # Top middle block (equal to arc_probs[2,1])

arc_probs[1,3] = arc_probs[1,1]
arc_probs[2,3] = arc_probs[2,1]
arc_probs[3,3] = arc_probs[1,1]
arc_probs[3,2] = arc_probs[2,1]
arc_probs[3,1] = arc_probs[1,1]

#Normalize
arc_probs = arc_probs./sum(arc_probs) # Same as arc_probs = arc_probs./(2π)

end