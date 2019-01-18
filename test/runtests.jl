using Test
using UAPIC
using FFTW

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

    tildex = zeros(ComplexF64, (ntau, 2, nbpart))
    tildey = zeros(ComplexF64, (ntau, 2, nbpart))

    Et    = zeros(Float64, (nbpart, ntau, 2))

    tilde = zeros(ComplexF64, (ntau, 2))
    r     = zeros(ComplexF64, (ntau, 2))
    xt    = zeros(ComplexF64, (ntau, 2, nbpart))
    yt    = zeros(ComplexF64, (ntau, 2, nbpart))

    fx1 = zeros(ComplexF64, (ntau, 2, nbpart))
    fx0 = zeros(ComplexF64, (ntau, 2, nbpart))
    fy1 = zeros(ComplexF64, (ntau, 2, nbpart))
    fy0 = zeros(ComplexF64, (ntau, 2, nbpart))

    auxpx[1,:] .= (particles.dx+particles.ix) * dx
    auxpx[2,:] .= (particles.dy+particles.iy) * dy

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

                xt1 = auxpx[1,m] + h1 + ep * vyb
                xt2 = auxpx[2,m] + h2 - ep * vxb

                xt[n,1,m] = xt1
                xt[n,2,m] = xt2

                interv=(1+0.5*sin(xt1)*sin(xt2)-b)/ep

                exb =((  cos(tau[n])*vy - sin(tau[n])*vx)
                       * interv + ex)/b

                eyb =(( - cos(tau[n])*vx - sin(tau[n])*vy)
                        * interv + ey)/b

                r[n,1] = cos(tau[n])* exb - sin(tau[n]) * eyb
                r[n,2] = sin(tau[n])* exb + cos(tau[n]) * eyb

            end

            fft!(r,1)
            for n = 2:ntau
                r[n,:] .= -1im * r[n,:]/ltau[n]
            end
            r[1,:] .= 0.0
            ifft!(r,1)

            yt[:,1,m] .= vx .+ (r[:,1] .- r[1,1]) * ep
            yt[:,2,m] .= vy .+ (r[:,2] .- r[1,2]) * ep

        end
        @show sum(yt)


        for n=1:ntau
            for m=1:nbpart
                xxt1, xxt2 = real(xt[n,1:2,m])
                @apply_bc()
                particles.ix[m] = trunc(Int32, xxt1/dimx*nx)
                particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
                particles.iy[m] = trunc(Int32, xxt2/dimy*ny)
                particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
            end
            interpol_eb_m6!( particles, fields )
            Et[:,n,1] .= particles.ex
            Et[:,n,2] .= particles.ey
        end
        @show sum(particles.ix .+ particles.dx)
        @show sum(particles.iy .+ particles.dy)
        @show sum(Et)

        #  !--time iteration
        #  prediction--
        for m=1:nbpart

            b = bx[m]
            for n=1:ntau

                fx0[n,1,m] = (   cos(tau[n]) * yt[n,1,m] 
                               + sin(tau[n]) * yt[n,2,m])/b
                fx0[n,2,m] = ( - sin(tau[n]) * yt[n,1,m] 
                               + cos(tau[n]) * yt[n,2,m])/b

                interv = (1 + 0.5*sin(real(xt[n,1,m]))*sin(real(xt[n,2,m]))-b)/ep

                temp1 = Et[m,n,1]+(  cos(tau[n])*yt[n,2,m]
                                       - sin(tau[n])*yt[n,1,m])*interv
                temp2 = Et[m,n,2]+(- cos(tau[n])*yt[n,1,m]
                                       - sin(tau[n])*yt[n,2,m])*interv

                fy0[n,1,m] = (cos(tau[n])*temp1-sin(tau[n])*temp2)/b
                fy0[n,2,m] = (sin(tau[n])*temp1+cos(tau[n])*temp2)/b

            end

        end

        fft!(fx0,1)
        fft!(fy0,1)

        tildex .= fft(xt,1) 
        tildey .= fft(yt,1) 

        for m=1:nbpart, i=1:2, n=1:ntau
            xt[n,i,m] = (exp.(-1im*ltau[n]*ds[m]/ep) * tildex[n,i,m]
                         + pl[n,m] * fx0[n,i,m])
        end

        ifft!(xt,1)

        for m=1:nbpart, i=1:2, n=1:ntau
            yt[n,i,m] = (exp(-1im*ltau[n]*ds[m]/ep) * tildey[n,i,m]
                         + pl[n,m]*fy0[n,i,m])
        end

        ifft!(yt,1) 

        @show sum(yt)

        for m = 1:nbpart
            time = ds[m]
            tilde .= fft(view(xt,:,:,m),1)
            tilde ./= ntau 
            tilde .*= exp.(1im*ltau*time/ep)
            xxt1, xxt2 = real(sum(tilde,dims=1))
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
                xxt1, xxt2 = real(xt[n,:,m])
                @apply_bc( )
                particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
                particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
                particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
                particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
            end

            interpol_eb_m6!( particles, fields )

            Et[:,n,1]=particles.ex
            Et[:,n,2]=particles.ey

        end

        @show sum(Et)

        # correction
        for m=1:nbpart, n=1:ntau

            fx1[n,1,m]=( cos(tau[n])*yt[n,1,m]+sin(tau[n])*yt[n,2,m])/bx[m]
            fx1[n,2,m]=(-sin(tau[n])*yt[n,1,m]+cos(tau[n])*yt[n,2,m])/bx[m]

            interv=(1 + 0.5*sin(real(xt[n,1,m]))*sin(real(xt[n,2,m]))
                       -bx[m])/ep

            temp1 = Et[m,n,1]+( cos(tau[n])*yt[n,2,m]
                               -sin(tau[n])*yt[n,1,m])*interv
            temp2 = Et[m,n,2]+(-cos(tau[n])*yt[n,1,m]
                               -sin(tau[n])*yt[n,2,m])*interv

            fy1[n,1,m] = (cos(tau[n])*temp1-sin(tau[n])*temp2)/bx[m]
            fy1[n,2,m] = (sin(tau[n])*temp1+cos(tau[n])*temp2)/bx[m]

        end

        fft!(fx1,1)
        fft!(fy1,1)

        for m=1:nbpart, i=1:2, n=1:ntau
        @inbounds xt[n,i,m] = ( exp(-1im*ltau[n]*ds[m]/ep)*tildex[n,i,m]
                      + pl[n,m] * fx0[n,i,m] + ql[n,m] * 
                      (fx1[n,i,m] - fx0[n,i,m]) / ds[m] )
        end

        ifft!(xt,1)

        for m=1:nbpart, i=1:2, n=1:ntau
        @inbounds yt[n,i,m] = ( exp(-1im*ltau[n]*ds[m]/ep)*tildey[n,i,m]
                      + pl[n,m]*fy0[n,i,m]
                      + ql[n,m]*(fy1[n,i,m]-fy0[n,i,m])/ds[m])
        end

        ifft!(yt,1) 

        @show sum(yt)

        for m=1:nbpart

            time = ds[m]
            tilde  .= fft(xt[:,:,m],1)
            tilde ./= ntau 
            tilde .*= exp.(1im*ltau*time/ep)

            xxt1, xxt2 = real(sum(tilde, dims=1))

            @apply_bc( )

            particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
            particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
            particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
            particles.dy[m] = Float32(xxt2/dy - particles.iy[m])

        end

        auxpx[1,:] .= (particles.dx + particles.ix) * dx
        auxpx[2,:] .= (particles.dy + particles.iy) * dy

        calcul_rho_m6!( fields, particles )

        @show nrj = poisson!( fields )

        for n=1:ntau

            for m=1:nbpart
                xxt1, xxt2 = real(xt[n,1:2,m])
                @apply_bc( )
                particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
                particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
                particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
                particles.dy[m] = Float32(xxt2/dy - particles.iy[m])
            end

            interpol_eb_m6!( particles, fields )

            Et[:,n,1] .= particles.ex
            Et[:,n,2] .= particles.ey

        end

        @show sum(Et)

        fft!(yt,1) 

        @simd for m=1:nbpart

            px = sum(view(yt,:,1,m)/ntau .* exp.(1im*ltau*ds[m]/ep))
            py = sum(view(yt,:,2,m)/ntau .* exp.(1im*ltau*ds[m]/ep))

            particles.vx[m] = real(cos(ds[m]/ep)*px+sin(ds[m]/ep)*py)
            particles.vy[m] = real(cos(ds[m]/ep)*py-sin(ds[m]/ep)*px)

        end

        @show sum(particles.vx), sum(particles.vy)

    end


    true

end 

const ntau = 16

@time test_pic2d( ntau )
