using PyPlot


function dx(r)
    # Parameters
    ωo = 6
    
    θ = r[1]
    ω = r[2]
    return [ω;-ωo*sin(θ)]
end
    

function RK(x)
    Δ = 1E-2        # s
    
    K1 = dx(x)
    K2 = dx(x .+ 0.5*Δ*K1)
    K3 = dx(x .+ 0.5*Δ*K2)
    K4 = dx(x .+ Δ*K3)

    return x + Δ*(K1 .+ 2*K2 .+ 2*K3 + K4)/6
end

function play(x)
    y = []
    z = []
    for it in 1:1000
                
        push!(y,x[1]*180/π)
        push!(z,x[2])
        
        x = RK(x)
    end
    
    return y,z
end

y,z = play([0.1;0])
plot(y,z)
# plot(y)
# xlabel("time [s]")
# ylabel("Angle [degrees]")
# 
# ax1 = gca()
# ax2 = ax1.twinx()
# ax2.plot(z,color="orange")
# ax2.set_ylabel("Position [m]")

show()
