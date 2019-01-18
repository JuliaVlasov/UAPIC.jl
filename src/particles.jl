
export Particles

mutable struct Particles

    nbpart :: Int64

    dx :: Vector{Float32}
    dy :: Vector{Float32}
    ix :: Vector{Int32}
    iy :: Vector{Int32}
    vx :: Vector{Float64}
    vy :: Vector{Float64}
    ex :: Vector{Float64}
    ey :: Vector{Float64}
    bz :: Vector{Float64}
    p  :: Float64

    function Particles( nbpart :: Int64, p :: Float64 )

        dx = zeros(Float32, nbpart)
        dy = zeros(Float32, nbpart)
        ix = zeros(Int32,   nbpart)
        iy = zeros(Int32,   nbpart)
        vx = zeros(Float64, nbpart)
        vy = zeros(Float64, nbpart)
        ex = zeros(Float64, nbpart)
        ey = zeros(Float64, nbpart)
        bz = zeros(Float64, nbpart)

        new( nbpart, dx, dy, ix, iy, vx, vy, ex, ey, bz, p )

    end

end

include("read_particles.jl")
include("plasma.jl")
include("compute_rho.jl")
include("interpolation.jl")

