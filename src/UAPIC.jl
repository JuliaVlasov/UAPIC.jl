module UAPIC

export Mesh

struct Mesh

    xmin :: Float64
    xmax :: Float64
    nx   :: Int32
    ymin :: Float64
    ymax :: Float64
    ny   :: Float32
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

    ex :: Array{Float64,2}
    ey :: Array{Float64,2}
    bz :: Array{Float64,2}
    ro :: Array{Float64,2}

    function MeshFields( nx :: Int64, ny :: Int64 )

         ex = zeros(Float64, (nx+1, ny+1))
         ey = zeros(Float64, (nx+1, ny+1))
         bz = zeros(Float64, (nx+1, ny+1))
         ro = zeros(Float64, (nx+1, ny+1))

         new( ex, ey, bz, ro)

    end

end

export Particle

mutable struct Particle

    num :: Int64
    dpx :: Vector{Float32}
    dpy :: Vector{Float32}
    idx :: Vector{Int32}
    idy :: Vector{Int32}
    vpx :: Vector{Float64}
    vpy :: Vector{Float64}
    epx :: Vector{Float64}
    epy :: Vector{Float64}
    bpz :: Vector{Float64}
    p   :: Float64

    function Particle( nbpart :: Int64, p :: Float64 )

        dpx = zeros(Float32, nbpart)
        dpy = zeros(Float32, nbpart)
        idx = zeros(Int32,   nbpart)
        idy = zeros(Int32,   nbpart)
        vpx = zeros(Float64, nbpart)
        vpy = zeros(Float64, nbpart)
        epx = zeros(Float64, nbpart)
        epy = zeros(Float64, nbpart)
        bpz = zeros(Float64, nbpart)

        new( nbpart, dpx, dpy, idx, idy, vpx, vpy, epx, epy, bpz, p )

    end

end

export plasma

function plasma( mesh :: Mesh )

    eps    = 1.e-12
    vth    = 1.0
    nbpart = 204800
    kx     = 0.5
    alpha  = 0.05

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx = mesh.xmax - mesh.xmin
    dimy = mesh.ymax - mesh.ymin
    
    weight = (dimx * dimy) / nbpart

    ele = Particle(nbpart, weight )
    
    k = 1
    while (k<=nbpart)

        xi   = rand() * dimx
        yi   = rand() * dimy 
        zi   = (2.0 + alpha) * rand()
        temm = 1.0 + sin(yi) + alpha * cos(kx*xi)

        if (temm>=zi)
            ele.idx[k] = floor(Int32, xi/dimx*nx)
            ele.idy[k] = floor(Int32, yi/dimy*ny)
            ele.dpx[k] = Float32(xi/dx - ele.idx[k])
            ele.dpy[k] = Float32(yi/dy - ele.idy[k])
            k = k + 1
        end
    end

    k = 1

    while (k<=nbpart)
        xi   = (rand()-0.5)*10.0
        yi   = (rand()-0.5)*10.0
        zi   = rand()
        temm = (exp(-((xi-2.0)^2+(yi-0.)^2)/2.0)+exp(-((xi+2)^2+yi^2)/2))/2

        if (temm>=zi)
            ele.vpx[k] = xi
            ele.vpy[k] = yi
            k = k + 1
        end
    end



end 

end # module
