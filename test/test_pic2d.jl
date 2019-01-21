using Test
using UAPIC
using FFTW
using LinearAlgebra

function update_particles!( particles :: Particles, 
                            fields    :: MeshFields, 
                            ua        :: UA, 
                            xt1       :: Array{ComplexF64,2}, 
                            xt2       :: Array{ComplexF64,2})

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

        mul!(ua.ftau, ua.ptau, view(xt1,:,m))

        for n = 1:ntau
            ua.ftau[n] *= exp(1im*ua.ltau[n]*t/ua.ε)/ua.ntau
        end 

        xxt1 = real(sum(ua.ftau))

        mul!(ua.ftau, ua.ptau, view(xt2,:,m))

        for n = 1:ntau
            ua.ftau[n] *= exp(1im*ua.ltau[n]*t/ua.ε)/ua.ntau
        end 

        xxt2 = real(sum(ua.ftau))

        particles.x[1,m] = xxt1
        particles.x[2,m] = xxt2

    end

    calcul_rho_m6!( fields, particles )

end

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

    tau  = ua.tau
    ltau = ua.ltau

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

    x̃t1 = zeros(ComplexF64, (ntau, nbpart))
    ỹt1 = zeros(ComplexF64, (ntau, nbpart))
    x̃t2 = zeros(ComplexF64, (ntau, nbpart))
    ỹt2 = zeros(ComplexF64, (ntau, nbpart))

    etx = zeros(Float64, (nbpart, ntau))
    ety = zeros(Float64, (nbpart, ntau))

    r  = zeros(ComplexF64, (ntau,2))
    r̃  = zeros(ComplexF64, (ntau,2))

    rtau = plan_fft(r,1)

    xt1 = zeros(ComplexF64, (ntau, nbpart))
    xt2 = zeros(ComplexF64, (ntau, nbpart))
    yt1 = zeros(ComplexF64, (ntau, nbpart))
    yt2 = zeros(ComplexF64, (ntau, nbpart))

    ftau = plan_fft(xt1, 1)

    gx1 = zeros(ComplexF64, (ntau, nbpart))
    fx1 = zeros(ComplexF64, (ntau, nbpart))
    gy1 = zeros(ComplexF64, (ntau, nbpart))
    fy1 = zeros(ComplexF64, (ntau, nbpart))

    gx2 = zeros(ComplexF64, (ntau, nbpart))
    fx2 = zeros(ComplexF64, (ntau, nbpart))
    gy2 = zeros(ComplexF64, (ntau, nbpart))
    fy2 = zeros(ComplexF64, (ntau, nbpart))

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

            for i=2:ntau
                pl[i,m] = ε * 1im*(exp(-1im*ltau[i]*t/ε)-1)/ltau[i]
                ql[i,m] = ε * (ε*(1-exp(-1im*ltau[i]*t/ε)) -1im*ltau[i]*t)/ltau[i]^2
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

                xxt1 = auxpx[1,m] + h1 + ε * vyb
                xxt2 = auxpx[2,m] + h2 - ε * vxb

                xt1[n,m] = xxt1
                xt2[n,m] = xxt2

                interv=(1+0.5*sin(xxt1)*sin(xxt2)-b)/ε

                exb = ((  cos(τ)*vy - sin(τ)*vx) * interv + ex)/b
                eyb = (( -cos(τ)*vx - sin(τ)*vy) * interv + ey)/b

                r[n,1] = cos(τ)* exb - sin(τ) * eyb
                r[n,2] = sin(τ)* exb + cos(τ) * eyb

            end

            mul!(r̃,rtau,r)

            for n = 2:ntau
                r̃[n,1] = -1im * r̃[n,1]/ltau[n]
            end
            for n = 2:ntau
                r̃[n,2] = -1im * r̃[n,2]/ltau[n]
            end

            ldiv!(r, rtau, r̃)

            for n = 1:ntau
                yt1[n,m] = vx + (r[n,1] - r[1,1]) * ε
            end
            for n = 1:ntau
                yt2[n,m] = vy + (r[n,2] - r[1,2]) * ε
            end

        end


        for n=1:ntau

            for m=1:nbpart
                xxt1 = real(xt1[n,m])
                xxt2 = real(xt2[n,m])
                particles.x[1,m] = xxt1
                particles.x[2,m] = xxt2
            end

            interpol_eb_m6!( particles, fields )

            for m=1:nbpart
                etx[m,n] = particles.e[1,m]
                ety[m,n] = particles.e[2,m]
            end

        end

        #  prediction--
        for m=1:nbpart

            b = particles.b[m]
            for n=1:ntau

                τ = tau[n]

                fx1[n,m] = ( cos(τ) * yt1[n,m] + sin(τ) * yt2[n,m])/b
                fx2[n,m] = (-sin(τ) * yt1[n,m] + cos(τ) * yt2[n,m])/b

                interv = (1 + 0.5*sin(real(xt1[n,m]))*sin(real(xt2[n,m]))-b)/ε

                tmp1 = etx[m,n]+(  cos(τ)*yt2[n,m] - sin(τ)*yt1[n,m])*interv
                tmp2 = ety[m,n]+(- cos(τ)*yt1[n,m] - sin(τ)*yt2[n,m])*interv

                fy1[n,m] = (cos(τ)*tmp1-sin(τ)*tmp2)/b
                fy2[n,m] = (sin(τ)*tmp1+cos(τ)*tmp2)/b

            end

        end

        fft!(fx1,1)
        fft!(fx2,1)
        fft!(fy1,1)
        fft!(fy2,1)

        mul!(x̃t1, ftau, xt1) 
        mul!(x̃t2, ftau, xt2) 

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                xt1[n,m] = (exp.(-1im*ltau[n]*t/ε) * x̃t1[n,m] + pl[n,m] * fx1[n,m])
                xt2[n,m] = (exp.(-1im*ltau[n]*t/ε) * x̃t2[n,m] + pl[n,m] * fx2[n,m])
            end
        end

        ifft!(xt1,1)
        ifft!(xt2,1)

        mul!(ỹt1, ftau, yt1) 
        mul!(ỹt2, ftau, yt2) 

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                yt1[n,m] = (exp(-1im*ltau[n]*t/ε) * ỹt1[n,m] + pl[n,m]*fy1[n,m])
                yt2[n,m] = (exp(-1im*ltau[n]*t/ε) * ỹt2[n,m] + pl[n,m]*fy2[n,m])
            end
        end

        ifft!(yt1,1) 
        ifft!(yt2,1) 

        update_particles!( particles, fields, ua, xt1, xt2)
        nrj = poisson!( fields ) 

        for n=1:ntau

            for m=1:nbpart
                xxt1 = real(xt1[n,m])
                xxt2 = real(xt2[n,m])
                particles.x[1,m] = xxt1
                particles.x[2,m] = xxt2
            end

            interpol_eb_m6!( particles, fields )

            for m=1:nbpart
                etx[m,n] = particles.e[1,m]
                ety[m,n] = particles.e[2,m]
            end

        end

        # correction
        for m=1:nbpart

            b = particles.b[m]

            for n = 1:ntau

                τ = tau[n]
                
                gx1[n,m]=( cos(τ)*yt1[n,m]+sin(τ)*yt2[n,m])/b
                gx2[n,m]=(-sin(τ)*yt1[n,m]+cos(τ)*yt2[n,m])/b
    
                interv=(1+0.5*sin(real(xt1[n,m]))*sin(real(xt2[n,m]))-b)/ε
    
                temp1 = etx[m,n]+( cos(τ)*yt2[n,m]-sin(τ)*yt1[n,m])*interv
                temp2 = ety[m,n]+(-cos(τ)*yt1[n,m]-sin(τ)*yt2[n,m])*interv
    
                gy1[n,m] = (cos(τ)*temp1-sin(τ)*temp2)/b
                gy2[n,m] = (sin(τ)*temp1+cos(τ)*temp2)/b

            end

        end

        fft!(gx1,1)
        fft!(gy1,1)
        fft!(gx2,1)
        fft!(gy2,1)

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                xt1[n,m] = ( exp(-1im*ltau[n]*t/ε)*x̃t1[n,m]
                          + pl[n,m] * fx1[n,m] + ql[n,m] * 
                          (gx1[n,m] - fx1[n,m]) / t )
                xt2[n,m] = ( exp(-1im*ltau[n]*t/ε)*x̃t2[n,m]
                          + pl[n,m] * fx2[n,m] + ql[n,m] * 
                          (gx2[n,m] - fx2[n,m]) / t )
            end
        end

        ifft!(xt1,1)
        ifft!(xt2,1)

        for m=1:nbpart
            t = particles.t[m]
            for n=1:ntau
                yt1[n,m] = ( exp(-1im*ltau[n]*t/ε)*ỹt1[n,m]
                          + pl[n,m]*fy1[n,m]
                          + ql[n,m]*(gy1[n,m]-fy1[n,m])/t)
                yt2[n,m] = ( exp(-1im*ltau[n]*t/ε)*ỹt2[n,m]
                          + pl[n,m]*fy2[n,m]
                          + ql[n,m]*(gy2[n,m]-fy2[n,m])/t)
            end
        end

        ifft!(yt1,1) 
        ifft!(yt2,1) 

        update_particles!( particles, fields, ua, xt1, xt2)

        @simd for m = 1:nbpart
            @inbounds auxpx[1,m] = particles.x[1,m]
            @inbounds auxpx[2,m] = particles.x[2,m]
        end

        nrj = poisson!( fields )

        for n=1:ntau

            for m=1:nbpart
                xxt1 = real(xt1[n,m])
                xxt2 = real(xt2[n,m])
                particles.x[1,m] = xxt1
                particles.x[2,m] = xxt2
            end

            interpol_eb_m6!( particles, fields )

            for m=1:nbpart
                etx[m,n] = particles.e[1,m]
                ety[m,n] = particles.e[2,m]
            end

        end

        fft!(yt1,1) 
        fft!(yt2,1) 

        for m=1:nbpart
            t = particles.t[m]

            px = 0.0
            py = 0.0
            for n = 1:ntau
                elt = exp(1im*ltau[n]*t/ε) 
                px += yt1[n,m]/ntau * elt
                py += yt2[n,m]/ntau * elt
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
