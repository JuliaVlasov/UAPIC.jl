export Particle, Particles

struct Particle

    px :: Float64
    py :: Float64
    vx :: ComplexF64
    vy :: ComplexF64
    ex :: ComplexF64
    ey :: ComplexF64
    bx :: Float64
    t  :: Float64

end

mutable struct Particles

    nbpart :: Int64

    x :: Array{Float64,2}
    v :: Array{Float64,2}
    e :: Array{Float64,2}
    b :: Vector{Float64}
    t :: Vector{Float64}
    w :: Float64

    function Particles( nbpart :: Int64, w :: Float64 )

        x = zeros(Float64, (2,nbpart))
        v = zeros(Float64, (2,nbpart))
        e = zeros(Float64, (2,nbpart))
        b = zeros(Float64, nbpart)
        t = zeros(Float64, nbpart)

        new( nbpart, x, v, e, b, t , w )

    end

end

include("read_particles.jl")
include("plasma.jl")
include("compute_rho.jl")
include("interpolation.jl")

