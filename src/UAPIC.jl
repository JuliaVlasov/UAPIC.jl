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
  
    ntau :: Int64
    tau  :: Vector{Float64}
    ltau :: Vector{Float64}
    ptau :: FFTW.cFFTWPlan{Complex{Float64},-1,false,1}

    function UA( ntau )

        dtau = 2π / ntau
        
        ltau  = zeros(Float64, ntau)
        ltau .= vcat(0:ntau÷2-1, -ntau÷2:-1) 
        
        tau   = zeros(Float64, ntau)
        tau  .= [ i*dtau for i=0:ntau-1 ]

        ptau  = plan_fft(tau)
        new( ntau, tau, ltau, ptau )

    end

end
    

include("integrate.jl")
include("gnuplot.jl")
include("poisson.jl")
include("particles.jl")

end # module
