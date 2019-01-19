
export Particles

mutable struct Particles

    nbpart :: Int64

    px :: Vector{Float64}
    py :: Vector{Float64}
    vx :: Vector{Float64}
    vy :: Vector{Float64}
    ex :: Vector{Float64}
    ey :: Vector{Float64}
    bx :: Vector{Float64}
    t  :: Vector{Float64}
    w  :: Float64

    function Particles( nbpart :: Int64, w :: Float64 )

        px = zeros(Float64, nbpart)
        py = zeros(Float64, nbpart)
        vx = zeros(Float64, nbpart)
        vy = zeros(Float64, nbpart)
        ex = zeros(Float64, nbpart)
        ey = zeros(Float64, nbpart)
        bx = zeros(Float64, nbpart)
        t  = zeros(Float64, nbpart)

        new( nbpart, px, py, vx, vy, ex, ey, bx, t , w )

    end

end

include("read_particles.jl")
include("plasma.jl")
include("compute_rho.jl")
include("interpolation.jl")

