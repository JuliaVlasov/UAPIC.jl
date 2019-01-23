using Test
using UAPIC
using FFTW
using LinearAlgebra


function test_pic2d( ntau )

    kx       = 0.50
    ky       = 1.0
    dimx     = 2π/kx
    dimy     = 2π/ky 
    nx       = 128	
    ny       = 64 
    tfinal   = 1.0 

    t = 0.

    xmin, xmax = 0.0, dimx
    ymin, ymax = 0.0, dimy

    mesh = Mesh( xmin, xmax, nx, ymin, ymax, ny )

    dx, dy = mesh.dx, mesh.dy

    dt = π / 2 / (2^3) #Used to determine N in MRC, not the macro-step
    tfinal = π / 2

    nstep  = trunc(Int64, tfinal/dt)

    fields = MeshFields( mesh )
    
    particles = read_particles( "particles.dat", mesh )

    nbpart = particles.nbpart

    poisson! = Poisson(mesh)

    ε = 0.1
    
    ua = UA( ntau, ε, nbpart )

    tau  = ua.tau
    ltau = ua.ltau

    et = zeros(Float64, (2, nbpart, ntau))

    xt  = zeros(ComplexF64, (ntau, 2, nbpart))
    x̃t  = zeros(ComplexF64, (ntau, 2, nbpart))
    yt  = zeros(ComplexF64, (ntau, 2, nbpart))
    ỹt  = zeros(ComplexF64, (ntau, 2, nbpart))

    ftau = plan_fft(xt,  1)

    fx = zeros(ComplexF64, (ntau, 2, nbpart))
    fy = zeros(ComplexF64, (ntau, 2, nbpart))

    gx = zeros(ComplexF64, (ntau, 2, nbpart))
    gy = zeros(ComplexF64, (ntau, 2, nbpart))

    calcul_rho_m6!( fields, particles )

    nrj =  poisson!( fields )

    interpol_eb_m6!( particles, fields )

    for istep = 1:2

        preparation!( ua, dt, particles, xt, yt) 

        update_particles_e!( particles, et, fields, ua, xt)

        #  prediction

        compute_f!( fx, fy, ua, particles, xt, yt, et )

        mul!(x̃t, ftau, xt) 

        ua_step!( xt, x̃t, ua, particles, fx )

        ifft!(xt,1)

        mul!(ỹt, ftau, yt) 

        ua_step!( yt, ỹt, ua, particles, fy )

        ifft!(yt,1) 

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields ) 

        update_particles_e!( particles, et, fields, ua, xt)

        # correction

        compute_f!( gx, gy, ua, particles, xt, yt, et )

        ua_step!( xt, x̃t, ua, particles, fx ) 

        ua_step!( yt, ỹt, ua, particles, fy )

        for m=1:nbpart

            t = particles.t[m]

            for n=1:ntau

                xt[n,1,m] += ua.ql[n,m] * (gx[n,1,m] - fx[n,1,m]) / t
                xt[n,2,m] += ua.ql[n,m] * (gx[n,2,m] - fx[n,2,m]) / t
                yt[n,1,m] += ua.ql[n,m] * (gy[n,1,m] - fy[n,1,m]) / t
                yt[n,2,m] += ua.ql[n,m] * (gy[n,2,m] - fy[n,2,m]) / t

            end

        end

        ifft!(xt,1)

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields )

        update_particles_e!( particles, et, fields, ua, xt)

        for m=1:nbpart

            t = particles.t[m]

            px, py = 0.0, 0.0
            for n = 1:ntau
                elt = exp(1im*ltau[n]*t/ε) 
                px += yt[n,1,m]/ntau * elt
                py += yt[n,2,m]/ntau * elt
            end

            particles.v[1,m] = real(cos(t/ε)*px+sin(t/ε)*py)
            particles.v[2,m] = real(cos(t/ε)*py-sin(t/ε)*px)

        end


    end

    @show @views sum(particles.v[1,:]), sum(particles.v[2,:])

    true

end 

const ntau = 16

@time test_pic2d( ntau )
