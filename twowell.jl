using PyPlot

fnl(x,v) = 4*x*(1-x^2) - v      # Double well with dissipation
#fnl(x) = -2*x                  # Simple harmonic potential

# 4-th order Runge Kutta for two coupled
# differential equations
function RK(xᵢ,vᵢ,Δ,func)
    k₁ₓ = Δ*(vᵢ)
    k₁ᵥ = Δ*func(xᵢ,vᵢ)

    k₂ₓ = Δ*(vᵢ + 0.5*k₁ᵥ)
    k₂ᵥ = Δ*func(xᵢ + 0.5*k₁ₓ, vᵢ + 0.5*k₁ᵥ)

    k₃ₓ = Δ*(vᵢ + 0.5*k₂ᵥ)
    k₃ᵥ = Δ*func(xᵢ + 0.5*k₂ₓ, vᵢ + 0.5*k₂ᵥ)

    k₄ₓ = Δ*(vᵢ + k₃ᵥ)
    k₄ᵥ = Δ*func(xᵢ + k₃ₓ, vᵢ + 0.5*k₃ᵥ)

    Δₓ = (k₁ₓ + 2*k₂ₓ + 2*k₃ₓ + k₄ₓ)/6
    Δᵥ = (k₁ᵥ + 2*k₂ᵥ + 2*k₃ᵥ + k₄ᵥ)/6

    return (xᵢ + Δₓ),(vᵢ + Δᵥ)
end

# Let the system evolve for 1500 steps
# given initial conditions x₀,v₀
function evolve(x₀,v₀)
    Δ = 0.01

    x = x₀
    v = v₀
    xs = []
    vs = []
    for i in 1:1500
        x,v = RK(x,v,Δ,fnl)

        push!(xs,x)
        push!(vs,v)
    end

    return xs,vs
end

# After letting the system evolve, find
# which attractor the system relaxes to
function endpoint(x₀,v₀)
    x,v = evolve(x₀,v₀)

    if (abs(x[end]-1) < abs(x[end]+1))
        return 1
    else
        return -1
    end
end

# Create a map of endpoints for different
# initial conditions
function scan(Nx,Nv)
    M = zeros(Nx,Nv)

    for nx in 1:Nx
        x₀ = -3 + 6*(nx-1)/(Nx-1)
        for nv in 1:Nv
            v₀ = -8 + 16*(nv-1)/(Nv-1)
            M[nx,nv] = endpoint(x₀,v₀)
        end
    end

    return M
end

# Plot phase space with appropriate
# parameters
function createImage(N,M)
    Δn = Int(N÷4)
    Δm = Int(M÷4)

    imshow(scan(N,M)',origin="lower",interpolation="bicubic")
    xticks([0,Δn,2Δn,3Δn,4Δn],["-3","-1.5","0","1.5","3"])
    yticks([0,Δm,2Δm,3Δm,4Δm],["-8","-4","0","4","8"])
    xlabel("x")
    ylabel("dx/dt")
    show()
end

#==========================#
# MAIN
#==========================#
createImage(200,200)

