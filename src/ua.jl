import FFTW
using LinearAlgebra

export UA

mutable struct UA
  
    ntau :: Int64
    ε    :: Float64
    tau  :: Vector{Float64}
    ltau :: Vector{Float64}
    ftau :: Array{ComplexF64,1}
    ptau :: FFTW.cFFTWPlan{ComplexF64,-1,false,1}

    function UA( ntau, ε )

        dtau = 2π / ntau
        
        ltau  = zeros(Float64, ntau)
        ltau .= vcat(0:ntau÷2-1, -ntau÷2:-1) 
        
        tau   = zeros(Float64, ntau)
        tau  .= [ i*dtau for i=0:ntau-1 ]

        ftau  = zeros(ComplexF64, ntau)

        ptau  = FFTW.plan_fft(ftau)

        new( ntau, ε, tau, ltau, ftau, ptau )

    end

end

export update_particles_e!

function update_particles_e!( particles :: Particles, 
                              et        :: Array{Float64,3}, 
                              fields    :: MeshFields, 
                              ua        :: UA, 
                              xt        :: Array{ComplexF64,3}) 
    for n=1:ua.ntau
    
        for m=1:particles.nbpart
            xt1 = real(xt[n,1,m])
            xt2 = real(xt[n,2,m])
            particles.x[1,m] = xt1
            particles.x[2,m] = xt2
        end
        
        interpol_eb_m6!( particles, fields )
        
        for m=1:particles.nbpart
            et[1,m,n] = particles.e[1,m]
            et[2,m,n] = particles.e[2,m]
        end
    
    end
    
end

export update_particles_x!

function update_particles_x!( particles :: Particles, 
                              fields    :: MeshFields, 
                              ua        :: UA, 
                              xt        :: Array{ComplexF64,3}) 

    nbpart = particles.nbpart
    nx     = fields.mesh.nx
    ny     = fields.mesh.ny
    dx     = fields.mesh.dx
    dy     = fields.mesh.dy
    xmin   = fields.mesh.xmin
    xmax   = fields.mesh.xmax 
    ymin   = fields.mesh.ymin
    ymax   = fields.mesh.ymax
    dimx   = xmax - xmin
    dimy   = ymax - ymin

    for m = 1:nbpart

        t = particles.t[m]

        mul!(ua.ftau, ua.ptau, view(xt,:,1,m))

        for n = 1:ua.ntau
            ua.ftau[n] *= exp(1im*ua.ltau[n]*t/ua.ε)/ua.ntau
        end 

        xt1 = real(sum(ua.ftau))

        mul!(ua.ftau, ua.ptau, view(xt,:,2,m))

        for n = 1:ua.ntau
            ua.ftau[n] *= exp(1im*ua.ltau[n]*t/ua.ε)/ua.ntau
        end 

        xt2 = real(sum(ua.ftau))

        particles.x[1,m] = xt1
        particles.x[2,m] = xt2

    end

    calcul_rho_m6!( fields, particles )

end
