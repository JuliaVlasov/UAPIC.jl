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
    ε = 0.1
    
    ua = UA( ntau, ε )

    tau  = view(ua.tau,:)
    ltau = view(ua.ltau,:)

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

    calcul_rho_m6!( fields, particles )

    nrj =  poisson!( fields )

    gnuplot("fields.dat", fields) 

    interpol_eb_m6!( particles, fields )

    auxpx = zeros(Float64, (2,nbpart))

    pl = zeros(ComplexF64, (ntau, nbpart))
    ql = zeros(ComplexF64, (ntau, nbpart))

    et = zeros(Float64, (2, nbpart, ntau))

    r  = zeros(ComplexF64, (ntau,2))
    r̃  = zeros(ComplexF64, (ntau,2))

    rtau = plan_fft(r,1)

    xt  = zeros(ComplexF64, (ntau, 2, nbpart))
    x̃t  = zeros(ComplexF64, (ntau, 2, nbpart))
    yt  = zeros(ComplexF64, (ntau, 2, nbpart))
    ỹt  = zeros(ComplexF64, (ntau, 2, nbpart))

    ftau = plan_fft(xt,  1)

    fx1 = zeros(ComplexF64, (ntau, nbpart))
    fy1 = zeros(ComplexF64, (ntau, nbpart))
    fx2 = zeros(ComplexF64, (ntau, nbpart))
    fy2 = zeros(ComplexF64, (ntau, nbpart))

    gx = zeros(ComplexF64, (ntau, 2))
    gy = zeros(ComplexF64, (ntau, 2))
    g̃x = zeros(ComplexF64, (ntau, 2))
    g̃y = zeros(ComplexF64, (ntau, 2))

    for m = 1:nbpart
        auxpx[1,m] = particles.x[1,m]
        auxpx[2,m] = particles.x[2,m]
    end

    for istep = 1:2

        # preparation
        for m = 1:nbpart

            particles.b[m] = 1 + 0.5 * sin(auxpx[1,m]) * sin(auxpx[2,m])
            particles.t[m]  = dt * particles.b[m]

            pl[1,m] = particles.t[m]
            ql[1,m] = particles.t[m]^2 / 2

            t = particles.t[m]
            b = particles.b[m]

            for n=2:ntau
                elt = exp(-1im*ltau[n]*t/ε) 
                pl[n,m] = ε * 1im*(elt-1)/ltau[n]
                ql[n,m] = ε * (ε*(1-elt) -1im*ltau[n]*t)/ltau[n]^2
            end

            # preparation initial data

            ex  = particles.e[1,m]
            ey  = particles.e[2,m]
            vx  = particles.v[1,m]
            vy  = particles.v[2,m]
            vxb = vx/b
            vyb = vy/b

            for n = 1:ntau

                τ = tau[n]

                h1 = ε * (sin(τ) * vxb - cos(τ) * vyb)
                h2 = ε * (sin(τ) * vyb + cos(τ) * vxb)

                xt1 = auxpx[1,m] + h1 + ε * vyb
                xt2 = auxpx[2,m] + h2 - ε * vxb

                xt[n,1,m] = xt1
                xt[n,2,m] = xt2

                interv=(1+0.5*sin(xt1)*sin(xt2)-b)/ε

                exb = ((  cos(τ)*vy - sin(τ)*vx) * interv + ex)/b
                eyb = (( -cos(τ)*vx - sin(τ)*vy) * interv + ey)/b

                r[n,1] = cos(τ)* exb - sin(τ) * eyb
                r[n,2] = sin(τ)* exb + cos(τ) * eyb

            end

            mul!(r̃,rtau,r)

            for n = 2:ntau
                r̃[n,1] = -1im * r̃[n,1]/ltau[n]
                r̃[n,2] = -1im * r̃[n,2]/ltau[n]
            end

            ldiv!(r, rtau, r̃)

            for n = 1:ntau
                yt[n,1,m] = vx + (r[n,1] - r[1,1]) * ε
                yt[n,2,m] = vy + (r[n,2] - r[1,2]) * ε
            end

        end

        update_particles_e!( particles, et, fields, ua, xt)

        #  prediction--
        for m=1:nbpart

            b = particles.b[m]
            for n=1:ntau

                yt1 = yt[n,1,m] 
                yt2 = yt[n,2,m] 

                τ = tau[n]

                fx1[n,m] = ( cos(τ) * yt1 + sin(τ) * yt2)/b
                fx2[n,m] = (-sin(τ) * yt1 + cos(τ) * yt2)/b

                interv = (1 + 0.5*sin(real(xt[n,1,m]))*sin(real(xt[n,2,m]))-b)/ε

                tmp1 = et[1,m,n]+(  cos(τ)*yt2 - sin(τ)*yt1)*interv
                tmp2 = et[2,m,n]+(- cos(τ)*yt1 - sin(τ)*yt2)*interv

                fy1[n,m] = (cos(τ)*tmp1-sin(τ)*tmp2)/b
                fy2[n,m] = (sin(τ)*tmp1+cos(τ)*tmp2)/b

            end

        end

        fft!(fx1,1)
        fft!(fx2,1)
        fft!(fy1,1)
        fft!(fy2,1)

        mul!(x̃t, ftau, xt) 

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                elt = exp(-1im*ltau[n]*t/ε) 
                xt[n,1,m] = elt * x̃t[n,1,m] + pl[n,m] * fx1[n,m]
                xt[n,2,m] = elt * x̃t[n,2,m] + pl[n,m] * fx2[n,m]
            end
        end

        ifft!(xt,1)

        mul!(ỹt, ftau, yt) 

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                elt = exp(-1im*ltau[n]*t/ε) 
                yt[n,1,m] = elt * ỹt[n,1,m] + pl[n,m]*fy1[n,m]
                yt[n,2,m] = elt * ỹt[n,2,m] + pl[n,m]*fy2[n,m]
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

            mul!(g̃x, rtau, gx)
            mul!(g̃y, rtau, gy)

            for n=1:ntau

                elt = exp(-1im*ltau[n]*t/ε) 

                x̃t[n,1,m] = ( elt*x̃t[n,1,m] + pl[n,m] * fx1[n,m] + ql[n,m] * (g̃x[n,1] - fx1[n,m]) / t )
                x̃t[n,2,m] = ( elt*x̃t[n,2,m] + pl[n,m] * fx2[n,m] + ql[n,m] * (g̃x[n,2] - fx2[n,m]) / t )
                ỹt[n,1,m] = ( elt*ỹt[n,1,m] + pl[n,m] * fy1[n,m] + ql[n,m] * (g̃y[n,1] - fy1[n,m]) / t )
                ỹt[n,2,m] = ( elt*ỹt[n,2,m] + pl[n,m] * fy2[n,m] + ql[n,m] * (g̃y[n,2] - fy2[n,m]) / t )

            end
        end

        ldiv!(xt,ftau,x̃t)

        update_particles_x!( particles, fields, ua, xt)

        for m = 1:nbpart
            auxpx[1,m] = particles.x[1,m]
            auxpx[2,m] = particles.x[2,m]
        end

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
