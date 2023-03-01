using Plots


function dx(r)
    # Parameters
    x = r[1]
    ẋ = r[2]
    θ = r[3]
    ω = r[4]
    ω² = ω^2

    m = 0.1     # kg
    M = 1       # kg
    g = 9.8067  # m/s²
    l = 0.3     # m
    F = 0       # N

    βx = 0.1   # Linear friction coefficient
    βθ = 0.1   # Angular friction coefficient

    d = l*(M + m*sin(θ)^2)
    f1 = (F + m*sin(θ)*(l*ω² - g*cos(θ)) + βθ*m*cos(θ)*ω - βx*ẋ)*l/d
    f2 = (-cos(θ)*(F + m*l*ω²*sin(θ)) + g*(M + m)*sin(θ) - βθ*(m + M)*ω + βx*cos(θ)*ẋ)/d

    return [ẋ;f1;ω;f2]
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
    for it in 1:5000

        push!(y,x[3])
        push!(z,x[1])

        x = RK(x)
    end

    return y,z
end


y,z = play([0;0;0.01;0])

anim = @animate for i in 1:1000
    θ = y[i]
    x = z[i]
    r1 = [x,x+sin(θ)]
    r2 = [0,cos(θ)]
    plot(r1,r2,markershape = :square, markersize=[8, 1],xlims=(-1.2,1.2),ylims=(-1.2,1.2))
end

gif(anim,"test.gif",fps=50)
#xlabel("time [s]")
#ylabel("Angle [degrees]")

#ax1 = gca()
#ax2 = ax1.twinx()
#ax2.plot(z,color="orange")
#ax2.set_ylabel("Position [m]")

#show()
