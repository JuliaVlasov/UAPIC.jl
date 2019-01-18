module UAPIC

export Mesh

struct Mesh

    xmin :: Float64
    xmax :: Float64
    nx   :: Int32
    ymin :: Float64
    ymax :: Float64
    ny   :: Int32
    dx   :: Float64
    dy   :: Float64

    function Mesh( xmin, xmax, nx, ymin, ymax, ny )

        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny

        new( xmin, xmax, nx, ymin, ymax, ny, dx, dy )

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

include("gnuplot.jl")
include("poisson.jl")

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

export plasma

function plasma( mesh :: Mesh, nbpart :: Int64 )

    kx     = 0.5
    alpha  = 0.05

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx = mesh.xmax - mesh.xmin
    dimy = mesh.ymax - mesh.ymin
    
    println(dimx * dimy)
    weight = (dimx * dimy) / nbpart
    println( dimx, dimy )
    println( weight )

    particles = Particles(nbpart, weight )
    
    k = 1
    while (k<=nbpart)

        xi   = rand() * dimx
        yi   = rand() * dimy 
        zi   = (2.0 + alpha) * rand()
        temm = 1.0 + sin(yi) + alpha * cos(kx*xi)

        if (temm >= zi)
            particles.ix[k] = trunc(Int32, xi/dimx*nx)
            particles.iy[k] = trunc(Int32, yi/dimy*ny)
            particles.dx[k] = Float32(xi/dx - particles.ix[k])
            particles.dy[k] = Float32(yi/dy - particles.iy[k])
            k = k + 1
        end

    end

    k = 1
    while (k<=nbpart)

        xi   = (rand()-0.5)*10
        yi   = (rand()-0.5)*10
        zi   = rand()
        temm = ( exp(-((xi-2)^2 + yi^2)/2)
               + exp(-((xi+2)^2 + yi^2)/2))/2

        if (temm >= zi)
            particles.vx[k] = xi
            particles.vy[k] = yi
            k = k + 1
        end

    end

    particles

end 

include("compute_rho.jl")
include("interpolation.jl")

end # module
