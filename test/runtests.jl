using Test
using UAPIC
using FFTW
using LinearAlgebra

#include("test_poisson.jl")
#include("test_particles.jl")

"""
UA scheme for 4d VP in Fluid-scaling with b(x)
Update b(x(tn)) every step
"""

macro apply_bc()

    return esc( quote

        while ( xxt1 > xmax )
            xxt1 = xxt1 - dimx
        end

        while ( xxt1 < xmin )
            xxt1 = xxt1 + dimx
        end

        while ( xxt2 > ymax )
            xxt2 = xxt2  - dimy
        end

        while ( xxt2  < ymin )
            xxt2 = xxt2  + dimy
        end

    end)

end 



function test_pic2d( ntau )

    @show nstepmax = 20000	
    @show kx       = 0.50
    @show ky       = 1.0
    @show dimx     = 2π/kx
    @show dimy     = 2π/ky 
    @show nx       = 128	
    @show ny       = 64 
    @show tfinal   = 1.0 

    time = 0.
    ep   = 0.1
    dtau = 2π / ntau
    
    ltau  = zeros(Float64, ntau)
    ltau .= vcat(0:ntau÷2-1, -ntau÷2:-1) 
    
    tau   = zeros(Float64, ntau)
    tau  .= [ i*dtau for i=0:ntau-1 ]
    
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

    @show nrj =  poisson!( fields )

    gnuplot("fields.dat", fields) 

    interpol_eb_m6!( particles, fields )

    @show ep   
    @show dt   
    @show dx   
    @show dy   
    @show dimx 
    @show dimy 
    @show ntau 
    
    auxpx = zeros(Float64, (2,nbpart))

    bx = zeros(Float64, nbpart)
    ds = zeros(Float64, nbpart)
    pl = zeros(ComplexF64, (ntau, nbpart))
    ql = zeros(ComplexF64, (ntau, nbpart))

    x̃t1 = zeros(ComplexF64, (ntau, nbpart))
    ỹt1 = zeros(ComplexF64, (ntau, nbpart))
    x̃t2 = zeros(ComplexF64, (ntau, nbpart))
    ỹt2 = zeros(ComplexF64, (ntau, nbpart))

    Et1 = zeros(Float64, (nbpart, ntau))
    Et2 = zeros(Float64, (nbpart, ntau))

    tilde = zeros(ComplexF64, ntau)

    fft_plan = plan_fft(tilde)

    r1    = zeros(ComplexF64, ntau)
    r2    = zeros(ComplexF64, ntau)
    xt1   = zeros(ComplexF64, (ntau, nbpart))
    xt2   = zeros(ComplexF64, (ntau, nbpart))
    yt1   = zeros(ComplexF64, (ntau, nbpart))
    yt2   = zeros(ComplexF64, (ntau, nbpart))

    fx1 = zeros(ComplexF64, (ntau, 2, nbpart))
    fx0 = zeros(ComplexF64, (ntau, 2, nbpart))
    fy1 = zeros(ComplexF64, (ntau, 2, nbpart))
    fy0 = zeros(ComplexF64, (ntau, 2, nbpart))

    for m = 1:nbpart
        auxpx[1,m] = (particles.dx[m]+particles.ix[m]) * dx
        auxpx[2,m] = (particles.dy[m]+particles.iy[m]) * dy
    end

    for istep = 1:1

        # preparation
        for m = 1:nbpart

            bx[m]   = 1 + 0.5 * sin(auxpx[1,m]) * sin(auxpx[2,m])
            ds[m]   = dt * bx[m]
            pl[1,m] = ds[m]
            ql[1,m] = ds[m]^2 / 2

            for i=2:ntau
                pl[i,m] = ep * 1im*(exp(-1im*ltau[i]*ds[m]/ep)-1)/ltau[i]
                ql[i,m] = ep * (ep*(1-exp(-1im*ltau[i]*ds[m]/ep))
                                -1im*ltau[i]*ds[m])/ltau[i]^2
            end

            # preparation initial data

            ex  = particles.ex[m]
            ey  = particles.ey[m]
            vx  = particles.vx[m]
            vy  = particles.vy[m]
            b   = bx[m]
            vxb = vx/b
            vyb = vy/b

            for n = 1:ntau

                h1 = ep * (sin(tau[n]) * vxb - cos(tau[n]) * vyb)
                h2 = ep * (sin(tau[n]) * vyb + cos(tau[n]) * vxb)

                xxt1 = auxpx[1,m] + h1 + ep * vyb
                xxt2 = auxpx[2,m] + h2 - ep * vxb

                xt1[n,m] = xxt1
                xt2[n,m] = xxt2

                interv=(1+0.5*sin(xxt1)*sin(xxt2)-b)/ep

                exb =((  cos(tau[n])*vy - sin(tau[n])*vx)
                       * interv + ex)/b

                eyb =(( - cos(tau[n])*vx - sin(tau[n])*vy)
                        * interv + ey)/b

                r1[n] = cos(tau[n])* exb - sin(tau[n]) * eyb
                r2[n] = sin(tau[n])* exb + cos(tau[n]) * eyb

            end

            mul!(tilde,fft_plan,r1)
            tilde[1] = 0.0
            for n = 2:ntau
                tilde[n] = -1im * tilde[n]/ltau[n]
            end
            ldiv!(r1,fft_plan, tilde)

            mul!(tilde,fft_plan,r2)
            tilde[1] = 0.0
            for n = 2:ntau
                tilde[n] = -1im * tilde[n]/ltau[n]
            end
            ldiv!(r2,fft_plan,tilde)

            for n = 1:ntau
                yt1[n,m] = vx + (r1[n] - r1[1]) * ep
                yt2[n,m] = vy + (r2[n] - r2[1]) * ep
            end

        end
        @show sum(yt1) + sum(yt2)


        for n=1:ntau

            for m=1:nbpart
                xxt1 = real(xt1[n,m])
                xxt2 = real(xt2[n,m])
                @apply_bc()
                particles.ix[m] = trunc(Int32, xxt1/dimx*nx)
                particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
                particles.iy[m] = trunc(Int32, xxt2/dimy*ny)
                particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
            end

            interpol_eb_m6!( particles, fields )

            for m=1:nbpart
                Et1[m,n] = particles.ex[m]
                Et2[m,n] = particles.ey[m]
            end

        end
        @show sum(particles.ix .+ particles.dx)
        @show sum(particles.iy .+ particles.dy)
        @show sum(Et1) + sum(Et2)

        #  !--time iteration
        #  prediction--
        for m=1:nbpart

            b = bx[m]
            for n=1:ntau

                fx0[n,1,m] = (   cos(tau[n]) * yt1[n,m] 
                               + sin(tau[n]) * yt2[n,m])/b
                fx0[n,2,m] = ( - sin(tau[n]) * yt1[n,m] 
                               + cos(tau[n]) * yt2[n,m])/b

                interv = (1 + 0.5*sin(real(xt1[n,m]))*sin(real(xt2[n,m]))-b)/ep

                temp1 = Et1[m,n]+(   cos(tau[n])*yt2[n,m]
                                   - sin(tau[n])*yt1[n,m])*interv

                temp2 = Et2[m,n]+(- cos(tau[n])*yt1[n,m]
                                  - sin(tau[n])*yt2[n,m])*interv

                fy0[n,1,m] = (cos(tau[n])*temp1-sin(tau[n])*temp2)/b
                fy0[n,2,m] = (sin(tau[n])*temp1+cos(tau[n])*temp2)/b

            end

        end

        fft!(fx0,1)
        fft!(fy0,1)

        x̃t1 .= fft(xt1,1) 
        x̃t2 .= fft(xt2,1) 
        ỹt1 .= fft(yt1,1) 
        ỹt2 .= fft(yt2,1) 

        for m=1:nbpart, n=1:ntau
            xt1[n,m] = (exp.(-1im*ltau[n]*ds[m]/ep) * x̃t1[n,m]
                         + pl[n,m] * fx0[n,1,m])
            xt2[n,m] = (exp.(-1im*ltau[n]*ds[m]/ep) * x̃t2[n,m]
                         + pl[n,m] * fx0[n,2,m])
        end

        ifft!(xt1,1)
        ifft!(xt2,1)

        for m=1:nbpart, n=1:ntau
            yt1[n,m] = (exp(-1im*ltau[n]*ds[m]/ep) * ỹt1[n,m]
                         + pl[n,m]*fy0[n,1,m])
            yt2[n,m] = (exp(-1im*ltau[n]*ds[m]/ep) * ỹt2[n,m]
                         + pl[n,m]*fy0[n,2,m])
        end

        ifft!(yt1,1) 
        ifft!(yt2,1) 

        @show sum(yt1) + sum(yt2)

        for m = 1:nbpart
            time = ds[m]
            tilde .= fft(view(xt1,:,m),1)
            tilde ./= ntau 
            tilde .*= exp.(1im*ltau*time/ep)
            xxt1 = real(sum(tilde))
            tilde .= fft(view(xt2,:,m),1)
            tilde ./= ntau 
            tilde .*= exp.(1im*ltau*time/ep)
            xxt2 = real(sum(tilde))
            @apply_bc( )
            particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
            particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
            particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
            particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
        end

        calcul_rho_m6!( fields, particles )
        @show nrj = poisson!( fields ) 

        for n=1:ntau
            for m=1:nbpart
                xxt1 = real(xt1[n,m])
                xxt2 = real(xt2[n,m])
                @apply_bc( )
                particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
                particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
                particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
                particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
            end

            interpol_eb_m6!( particles, fields )

            for m=1:nbpart
                Et1[m,n] = particles.ex[m]
                Et2[m,n] = particles.ey[m]
            end

        end

        @show sum(Et1) + sum(Et2)

        # correction
        for m=1:nbpart, n=1:ntau

            fx1[n,1,m]=( cos(tau[n])*yt1[n,m]+sin(tau[n])*yt2[n,m])/bx[m]
            fx1[n,2,m]=(-sin(tau[n])*yt1[n,m]+cos(tau[n])*yt2[n,m])/bx[m]

            interv=(1+0.5*sin(real(xt1[n,m]))*sin(real(xt2[n,m]))-bx[m])/ep

            temp1 = Et1[m,n]+( cos(tau[n])*yt2[n,m]
                              -sin(tau[n])*yt1[n,m])*interv
            temp2 = Et2[m,n]+(-cos(tau[n])*yt1[n,m]
                              -sin(tau[n])*yt2[n,m])*interv

            fy1[n,1,m] = (cos(tau[n])*temp1-sin(tau[n])*temp2)/bx[m]
            fy1[n,2,m] = (sin(tau[n])*temp1+cos(tau[n])*temp2)/bx[m]

        end

        fft!(fx1,1)
        fft!(fy1,1)

        for m=1:nbpart, n=1:ntau
            xt1[n,m] = ( exp(-1im*ltau[n]*ds[m]/ep)*x̃t1[n,m]
                      + pl[n,m] * fx0[n,1,m] + ql[n,m] * 
                      (fx1[n,1,m] - fx0[n,1,m]) / ds[m] )
            xt2[n,m] = ( exp(-1im*ltau[n]*ds[m]/ep)*x̃t2[n,m]
                      + pl[n,m] * fx0[n,2,m] + ql[n,m] * 
                      (fx1[n,2,m] - fx0[n,2,m]) / ds[m] )
        end

        ifft!(xt1,1)
        ifft!(xt2,1)

        for m=1:nbpart, n=1:ntau
            yt1[n,m] = ( exp(-1im*ltau[n]*ds[m]/ep)*ỹt1[n,m]
                      + pl[n,m]*fy0[n,1,m]
                      + ql[n,m]*(fy1[n,1,m]-fy0[n,1,m])/ds[m])
            yt2[n,m] = ( exp(-1im*ltau[n]*ds[m]/ep)*ỹt2[n,m]
                      + pl[n,m]*fy0[n,2,m]
                      + ql[n,m]*(fy1[n,2,m]-fy0[n,2,m])/ds[m])
        end

        ifft!(yt1,1) 
        ifft!(yt2,1) 

        @show sum(yt1) + sum(yt2)

        for m=1:nbpart

            time = ds[m]

            mul!(tilde, fft_plan, view(xt1,:,m))
            tilde ./= ntau 
            tilde .*= exp.(1im*ltau*time/ep)
            xxt1    = real(sum(tilde))

            mul!(tilde, fft_plan, view(xt2,:,m))
            
            tilde ./= ntau 
            tilde .*= exp.(1im*ltau*time/ep)
            xxt2    = real(sum(tilde))

            @apply_bc( )

            particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
            particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
            particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
            particles.dy[m] = Float32(xxt2/dy - particles.iy[m])

            auxpx[1,m] = xxt1
            auxpx[2,m] = xxt2

        end


        calcul_rho_m6!( fields, particles )

        @show nrj = poisson!( fields )

        for n=1:ntau

            for m=1:nbpart
                xxt1 = real(xt1[n,m])
                xxt2 = real(xt2[n,m])
                @apply_bc( )
                particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
                particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
                particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
                particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
            end

            interpol_eb_m6!( particles, fields )

            for m=1:nbpart
                Et1[m,n] = particles.ex[m]
                Et2[m,n] = particles.ey[m]
            end

        end

        @show sum(Et1) + sum(Et2)

        fft!(yt1,1) 
        fft!(yt2,1) 

        @simd for m=1:nbpart

            px = sum(view(yt1,:,m)/ntau .* exp.(1im*ltau*ds[m]/ep))
            py = sum(view(yt2,:,m)/ntau .* exp.(1im*ltau*ds[m]/ep))

            particles.vx[m] = real(cos(ds[m]/ep)*px+sin(ds[m]/ep)*py)
            particles.vy[m] = real(cos(ds[m]/ep)*py-sin(ds[m]/ep)*px)

        end

        @show sum(particles.vx), sum(particles.vy)

    end


    true

end 

const ntau = 16

@time test_pic2d( ntau )
