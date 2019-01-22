module UAPIC

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
    pl   :: Array{ComplexF64,2}
    ql   :: Array{ComplexF64,2}
    r    :: Array{ComplexF64,2}
    r̃    :: Array{ComplexF64,2}
    rtau :: FFTW.cFFTWPlan{ComplexF64,-1,false,2}

    function UA( ntau, ε, nbpart )

        dtau = 2π / ntau
        
        ltau  = zeros(Float64, ntau)
        ltau .= vcat(0:ntau÷2-1, -ntau÷2:-1) 
        
        tau   = zeros(Float64, ntau)
        tau  .= [ i*dtau for i=0:ntau-1 ]

        ftau  = zeros(ComplexF64, ntau)

        ptau  = FFTW.plan_fft(ftau)

        pl = zeros(ComplexF64, (ntau, nbpart))
        ql = zeros(ComplexF64, (ntau, nbpart))

        r  = zeros(ComplexF64, (ntau,2))
        r̃  = zeros(ComplexF64, (ntau,2))

        rtau = plan_fft(r,1)

        new( ntau, ε, tau, ltau, ftau, ptau, pl, ql, r, r̃, rtau )

    end

end

include("meshfields.jl")
include("integrate.jl")
include("gnuplot.jl")
include("poisson.jl")
include("particles.jl")

export preparation!


function preparation!( ua         :: UA, 
                       dt         :: Float64, 
                       particles  :: Particles, 
                       xt         :: Array{ComplexF64,3}, 
                       yt         :: Array{ComplexF64,3}) 

    ε      = ua.ε
    ntau   = ua.ntau
    nbpart = particles.nbpart
    tau    = ua.tau
    ltau   = ua.ltau
    
    for m = 1:nbpart

        x1 = particles.x[1,m]
        x2 = particles.x[2,m]

        particles.b[m] = 1 + 0.5 * sin(x1) * sin(x2)
        particles.t[m] = dt * particles.b[m]

        ua.pl[1,m] = particles.t[m]
        ua.ql[1,m] = particles.t[m]^2 / 2

        t = particles.t[m]
        b = particles.b[m]

        for n=2:ntau
            elt = exp(-1im*ua.ltau[n]*t/ε) 
            ua.pl[n,m] = ε * 1im*(elt-1)/ua.ltau[n]
            ua.ql[n,m] = ε * (ε*(1-elt) -1im*ua.ltau[n]*t)/ua.ltau[n]^2
        end

        ex,  ey  = particles.e[1:2,m]
        vx,  vy  = particles.v[1:2,m]
        vxb, vyb = vx/b, vy/b

        for n = 1:ntau

            τ = tau[n]

            h1 = ε * (sin(τ) * vxb - cos(τ) * vyb)
            h2 = ε * (sin(τ) * vyb + cos(τ) * vxb)

            xt1 = x1 + h1 + ε * vyb
            xt2 = x2 + h2 - ε * vxb

            xt[n,1,m] = xt1
            xt[n,2,m] = xt2

            interv=(1+0.5*sin(xt1)*sin(xt2)-b)/ε

            exb = ((  cos(τ)*vy - sin(τ)*vx) * interv + ex)/b
            eyb = (( -cos(τ)*vx - sin(τ)*vy) * interv + ey)/b

            ua.r[n,1] = cos(τ)* exb - sin(τ) * eyb
            ua.r[n,2] = sin(τ)* exb + cos(τ) * eyb

        end

        mul!(ua.r̃, ua.rtau, ua.r)

        for n = 2:ntau
            ua.r̃[n,1] = -1im * ua.r̃[n,1]/ltau[n]
            ua.r̃[n,2] = -1im * ua.r̃[n,2]/ltau[n]
        end

        ldiv!(ua.r, ua.rtau, ua.r̃)

        for n = 1:ntau
            yt[n,1,m] = vx + (ua.r[n,1] - ua.r[1,1]) * ε
            yt[n,2,m] = vy + (ua.r[n,2] - ua.r[1,2]) * ε
        end

    end

end

export update_particles_e!

function update_particles_e!( particles :: Particles, 
                              et        :: Array{Float64,3}, 
                              fields    :: MeshFields, 
                              ua        :: UA, 
                              xt        :: Array{ComplexF64,3}) 

    interpol_eb_m6!( et, fields, xt, particles.nbpart, ua.ntau) 

#    for n=1:ua.ntau
#    
#        for m=1:particles.nbpart
#            xt1 = real(xt[n,1,m])
#            xt2 = real(xt[n,2,m])
#            #particles.x[1,m] = xt1
#            #particles.x[2,m] = xt2
#        end
#        
#        #interpol_eb_m6!( particles, fields )
#        
#        #for m=1:particles.nbpart
#        #    et[1,m,n] = particles.e[1,m]
#        #    et[2,m,n] = particles.e[2,m]
#        #end
#    
#    end
    
end

export update_particles_x!

function update_particles_x!( particles :: Particles, 
                              fields    :: MeshFields, 
                              ua        :: UA, 
                              xt        :: Array{ComplexF64,3}) 

    #nbpart = particles.nbpart
    #nx     = fields.mesh.nx
    #ny     = fields.mesh.ny
    #dx     = fields.mesh.dx
    #dy     = fields.mesh.dy
    #xmin   = fields.mesh.xmin
    #xmax   = fields.mesh.xmax 
    #ymin   = fields.mesh.ymin
    #ymax   = fields.mesh.ymax
    #dimx   = xmax - xmin
    #dimy   = ymax - ymin

    #for m = 1:nbpart

    #    t = particles.t[m]

    #    mul!(ua.ftau, ua.ptau, view(xt,:,1,m))

    #    for n = 1:ua.ntau
    #        ua.ftau[n] *= exp(1im*ua.ltau[n]*t/ua.ε)/ua.ntau
    #    end 

    #    xt1 = real(sum(ua.ftau))

    #    mul!(ua.ftau, ua.ptau, view(xt,:,2,m))

    #    for n = 1:ua.ntau
    #        ua.ftau[n] *= exp(1im*ua.ltau[n]*t/ua.ε)/ua.ntau
    #    end 

    #    xt2 = real(sum(ua.ftau))

    #    particles.x[1,m] = xt1
    #    particles.x[2,m] = xt2

    #end

    #calcul_rho_m6!( fields, particles )
    calcul_rho_m6!( fields, particles, xt, ua )

end

end # module
