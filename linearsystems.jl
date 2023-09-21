using PyPlot
using LinearAlgebra


# 4th-order Runge Kutta
function RungeKutta(A,r,Δ)
    K1 = A*r
    K2 = A*(r .+ 0.5*Δ*K1)
    K3 = A*(r .+ 0.5*Δ*K2)
    K4 = A*(r .+ Δ*K3)

    return r + Δ*(K1 .+ 2K2 .+ 2K3 + K4)/6
end

function trajectory(rule,Δt,rₒ)
    r = rₒ
    x = []
    y = []

    for it in 1:1000
        r = RungeKutta(rule,r,Δt)

        push!(x,r[1])
        push!(y,r[2])
    end

    return x,y
end


equilibrium_case = [
    ("Stable Node", [-2 0; 0 -3]),
    ("Saddle Point", [1 0; 0 -1]),
    ("Focus or Spiral Point", [-2 -1; 1 -2]),
    ("Improper Node", [1 -16; 1 -7]),
    ("Center or Vortex", [0 -1; 1 0]),
    ("Degenerate Node", [-1 0; 1 0]),
    ("Source, Sink, or Proper Node", [-2 0; 0 -2])
]

for num in 1:7
    tit,A = equilibrium_case[num]

    p = tr(A)
    q = det(A)

    # Plot results
    println(tit)
    println("p....: ",p)
    println("q....: ",q)
    println("p^2..: ",p^2)
    println("4q...: ",4q)
    println("s1,s2: ",eigvals(A))
    println("-------------------")

    X = [0.1 0.1 -0.1 -0.1]
    Y = [0.1 -0.1 0.1 -0.1]
    C = ["blue","orange","green","red"]

    for i in 1:4
        x,y = trajectory(A,0.01,[X[i];Y[i]])
        plot([X[i]],[Y[i]],color=C[i],"o")
        plot(x,y,color=C[i])
    end

    axis([-0.2, 0.2, -0.2, 0.2])
    grid(true)

    X = LinRange(-0.2,0.2,10)
    Y = LinRange(-0.2,0.2,10)

    dx = zeros(10,10)
    dy = zeros(10,10)

    for x in 1:10, y in 1:10
        r = [X[x];Y[y]]
        d = RungeKutta(A,r,0.01) - r

        dx[y,x] = d[1]
        dy[y,x] = d[2]
    end

    quiver(X,Y,dx,dy,angles="uv")
    title(tit)
    show()
end
