using Plots

function dx(r,F)
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

    βx = 0.1   # Linear friction coefficient
    βθ = 0.1   # Angular friction coefficient

    d = l*(M + m*sin(θ)^2)
    f1 = (F + m*sin(θ)*(l*ω² - g*cos(θ)) + βθ*m*cos(θ)*ω - βx*ẋ)*l/d
    f2 = (-cos(θ)*(F + m*l*ω²*sin(θ)) + g*(M + m)*sin(θ) - βθ*(m + M)*ω + βx*cos(θ)*ẋ)/d

    return [ẋ;f1;ω;f2]
end

function RK(x,F)
    Δ = 1E-3        # s

    K1 = dx(x,F)
    K2 = dx(x .+ 0.5*Δ*K1,F)
    K3 = dx(x .+ 0.5*Δ*K2,F)
    K4 = dx(x .+ Δ*K3,F)

    return x + Δ*(K1 .+ 2*K2 .+ 2*K3 + K4)/6
end

function play(x,N)
    m = 0.1     # kg
    M = 1       # kg
    g = 9.8067  # m/s²
    l = 0.3     # m

    y = []
    z = []

    F = 0
    for m in 1:30
        for it in 1:50
            push!(y,x[3])
            push!(z,x[1])

            x = RK(x,F)
        end
        for it in 1:4
            push!(y,x[3])
            push!(z,x[1])

            x = RK(x,0)
        end

        h2 = 0.001
        h1 = g*(m+M)-h2^2/(4*l*M)
        F = h1*x[3]+h2*x[4]
    end

    return y,z
end

y,z = play([0;0;0.01;0],20000)

if (true)
    println("Plots")

    anim = @animate for i in 1:length(y)
        θ = y[i]
        x = z[i]
        r1 = [x,x+sin(θ)]
        r2 = [0,cos(θ)]
        plot(r1,r2,markershape = :square, markersize=[8, 1],ylims=(-1.2,1.2),xlims=(-2.5,2.5))
    end

    gif(anim,"test.gif",fps=450)
else
    println("PyPlot")
    using PyPlot

    plot(y)
    xlabel("time [s]")
    ylabel("Angle [degrees]")

    ax1 = gca()
    ax2 = ax1.twinx()
    ax2.plot(z,color="orange")
    ax2.set_ylabel("Position [m]")

    show()
end
