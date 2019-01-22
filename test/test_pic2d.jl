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

    gx = zeros(ComplexF64, (ntau, 2))
    gy = zeros(ComplexF64, (ntau, 2))
    g̃x = zeros(ComplexF64, (ntau, 2))
    g̃y = zeros(ComplexF64, (ntau, 2))

    calcul_rho_m6!( fields, particles )

    nrj =  poisson!( fields )

    interpol_eb_m6!( particles, fields )

    for istep = 1:2

        preparation!( ua, dt, particles, xt, yt) 

        update_particles_e!( particles, et, fields, ua, xt)

        #  prediction--
        for m=1:nbpart

            b = particles.b[m]
            for n=1:ntau

                xt1 = real(xt[n,1,m])
                xt2 = real(xt[n,2,m])

                yt1 = yt[n,1,m] 
                yt2 = yt[n,2,m] 

                τ = tau[n]

                fx[n,1,m] = ( cos(τ) * yt1 + sin(τ) * yt2)/b
                fx[n,2,m] = (-sin(τ) * yt1 + cos(τ) * yt2)/b

                interv = (1 + 0.5*sin(xt1)*sin(xt2)-b)/ε

                tmp1 = et[1,m,n]+(  cos(τ)*yt2 - sin(τ)*yt1)*interv
                tmp2 = et[2,m,n]+(- cos(τ)*yt1 - sin(τ)*yt2)*interv

                fy[n,1,m] = (cos(τ)*tmp1-sin(τ)*tmp2)/b
                fy[n,2,m] = (sin(τ)*tmp1+cos(τ)*tmp2)/b

            end

        end

        fft!(fx,1)
        fft!(fy,1)

        mul!(x̃t, ftau, xt) 

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                elt = exp(-1im*ltau[n]*t/ε) 
                xt[n,1,m] = elt * x̃t[n,1,m] + ua.pl[n,m] * fx[n,1,m]
                xt[n,2,m] = elt * x̃t[n,2,m] + ua.pl[n,m] * fx[n,2,m]
            end
        end

        ifft!(xt,1)

        mul!(ỹt, ftau, yt) 

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                elt = exp(-1im*ltau[n]*t/ε) 
                yt[n,1,m] = elt * ỹt[n,1,m] + ua.pl[n,m]*fy[n,1,m]
                yt[n,2,m] = elt * ỹt[n,2,m] + ua.pl[n,m]*fy[n,2,m]
            end
        end

        ifft!(yt,1) 

        update_particles_x!( particles, fields, ua, xt)
        nrj = poisson!( fields ) 
        update_particles_e!( particles, et, fields, ua, xt)

        # correction
        for m=1:nbpart

            b = particles.b[m]
            t = particles.t[m]

            for n = 1:ntau

                τ = tau[n]

                xt1 = xt[n,1,m] 
                xt2 = xt[n,2,m] 
                yt1 = yt[n,1,m] 
                yt2 = yt[n,2,m] 
                
                gx[n,1]=( cos(τ)*yt1+sin(τ)*yt2)/b
                gx[n,2]=(-sin(τ)*yt1+cos(τ)*yt2)/b
    
                interv=(1+0.5*sin(real(xt1))*sin(real(xt2))-b)/ε
    
                tmp1 = et[1,m,n]+( cos(τ)*yt2-sin(τ)*yt1)*interv
                tmp2 = et[2,m,n]+(-cos(τ)*yt1-sin(τ)*yt2)*interv
    
                gy[n,1] = (cos(τ)*tmp1-sin(τ)*tmp2)/b
                gy[n,2] = (sin(τ)*tmp1+cos(τ)*tmp2)/b

            end

            mul!(g̃x, ua.rtau, gx)
            mul!(g̃y, ua.rtau, gy)

            for n=1:ntau

                elt = exp(-1im*ltau[n]*t/ε) 

                fx1 = fx[n,1,m]
                fx2 = fx[n,2,m]
                fy1 = fy[n,1,m]
                fy2 = fy[n,2,m]

                x̃t[n,1,m] = elt*x̃t[n,1,m] + ua.pl[n,m] * fx1 + ua.ql[n,m] * (g̃x[n,1] - fx1) / t
                x̃t[n,2,m] = elt*x̃t[n,2,m] + ua.pl[n,m] * fx2 + ua.ql[n,m] * (g̃x[n,2] - fx2) / t
                ỹt[n,1,m] = elt*ỹt[n,1,m] + ua.pl[n,m] * fy1 + ua.ql[n,m] * (g̃y[n,1] - fy1) / t
                ỹt[n,2,m] = elt*ỹt[n,2,m] + ua.pl[n,m] * fy2 + ua.ql[n,m] * (g̃y[n,2] - fy2) / t

            end
        end

        ldiv!(xt,ftau,x̃t)

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields )

        update_particles_e!( particles, et, fields, ua, xt)

        for m=1:nbpart

            t = particles.t[m]

            px = 0.0
            py = 0.0
            for n = 1:ntau
                elt = exp(1im*ltau[n]*t/ε) 
                px += ỹt[n,1,m]/ntau * elt
                py += ỹt[n,2,m]/ntau * elt
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
