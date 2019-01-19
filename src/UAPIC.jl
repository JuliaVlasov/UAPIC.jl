module UAPIC

using FFTW

export Mesh

struct Mesh

    xmin :: Float64
    xmax :: Float64
    nx   :: Int32
    dx   :: Float64
    ymin :: Float64
    ymax :: Float64
    ny   :: Int32
    dy   :: Float64

    function Mesh( xmin, xmax, nx, ymin, ymax, ny )

        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny

        new( xmin, xmax, nx, dx, ymin, ymax, ny, dy )

    end

end


export MeshFields

mutable struct MeshFields

    mesh :: Mesh

    ex :: Array{Float64,2}
    ey :: Array{Float64,2}
    bz :: Array{Float64,2}
    ρ  :: Array{Float64,2}

    function MeshFields( mesh :: Mesh )
	    
	     nx, ny = mesh.nx, mesh.ny

         ex = zeros(Float64, (nx+1, ny+1))
         ey = zeros(Float64, (nx+1, ny+1))
         bz = zeros(Float64, (nx+1, ny+1))
         ρ  = zeros(Float64, (nx+1, ny+1))

         new( mesh, ex, ey, bz, ρ)

    end

end

export UA

struct UA
  
    nτ :: Int64
    ε  :: Float64
    τ  :: Vector{Float64}
    lτ :: Vector{Float64}
    pτ :: FFTW.cFFTWPlan{Complex{Float64},-1,false,1}

    function UA( nτ, ε )

        dτ = 2π / nτ
        
        lτ  = zeros(Float64, nτ)
        lτ .= vcat(0:nτ÷2-1, -nτ÷2:-1) 
        
        τ   = zeros(Float64, nτ)
        τ  .= [ i*dτ for i=0:nτ-1 ]

        pτ  = plan_fft(τ)
        new( nτ, ε, τ, lτ, pτ )

    end

end
    

include("integrate.jl")
include("gnuplot.jl")
include("poisson.jl")
include("particles.jl")

end # module
