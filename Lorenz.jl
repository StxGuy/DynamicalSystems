using Plots

# 4th-order Runge Kutta
function RungeKutta(dr,r,t,Δ)
    K1 = dr(r, t)
    K2 = dr(r .+ 0.5*Δ*K1, t + 0.5*Δ)
    K3 = dr(r .+ 0.5*Δ*K2, t + 0.5*Δ)
    K4 = dr(r .+ Δ*K3, t + Δ)

    return r + Δ*(K1 .+ 2*K2 .+ 2*K3 + K4)/6
end

# Lorenz attractor
function lorenz(u,t)
    σ,ρ,β = [10.0,28.0,8.0/3]
    x,y,z = u

    dx = σ*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z

    return [dx,dy,dz]
end

# Calculate the trajector from 0 to T
# time step Δt, and initial condition uₒ
function trajectory(rule,T,Δt,uₒ)
    N = Int(T÷Δt)

    plt = plot3d(1,xlim=(-30,30),ylim=(-30,30),zlim=(0,60),title="Lorenz",legend=false,marker=2)

    u = uₒ
    anim = @animate for nt in 1:N
        t = (nt-1)*Δt
        u = RungeKutta(rule,u,t,Δt)

        push!(plt,u[1],u[2],u[3])
    end every 10

    return anim
end

#--- MAIN ---
A = trajectory(lorenz,50,0.02,[1,1,1])
gif(A,"test.gif",fps=50)


