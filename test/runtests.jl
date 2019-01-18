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

macro energyuse()

    return esc( quote

        tilde[:,1] .= fft(xt[:,1,m])
        tilde[:,2] .= fft(xt[:,2,m])

        temp[1,:] .= 0.0
        for n=1:ntau
            temp[1,:] .+= tilde[n,:]/ntau * exp(1im*ltau[n]*time/ep)
        end

        xxt1, xxt2 = real(temp[1,:])

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

    Et = zeros(Float64, (2, ntau, nbpart))

    tilde = zeros(ComplexF64, (ntau, 2))
    temp  = zeros(ComplexF64, (ntau, 2))
    h     = zeros(ComplexF64, (ntau, 2))
    r     = zeros(ComplexF64, (ntau, 2))
    fx    = zeros(ComplexF64, (ntau, 2))
    fy    = zeros(ComplexF64, (ntau, 2))
    xt    = zeros(ComplexF64, (ntau, 2, nbpart))
    yt    = zeros(ComplexF64, (ntau, 2, nbpart))

    fxtemp1 = zeros(ComplexF64, (ntau, 2, nbpart))
    fxtemp0 = zeros(ComplexF64, (ntau, 2, nbpart))
    fytemp1 = zeros(ComplexF64, (ntau, 2, nbpart))
    fytemp0 = zeros(ComplexF64, (ntau, 2, nbpart))

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

            vxb = particles.vx[m]/bx[m]
            vyb = particles.vy[m]/bx[m]

            for n = 1:ntau
                h[n,1] = ep * (sin(tau[n]) * vxb
                             - cos(tau[n]) * vyb)
                h[n,2] = ep * (sin(tau[n]) * vyb 
                             + cos(tau[n]) * vxb)
                xt[n,1,m] = auxpx[1,m] .+ h[n,1] .- h[1,1]
                xt[n,2,m] = auxpx[2,m] .+ h[n,2] .- h[1,2]
            end


            for n = 1:ntau

                interv=(1+0.5*sin(real(xt[n,1,m]))
                             *sin(real(xt[n,2,m]))-bx[m])/ep

                exb =((  cos(tau[n])*particles.vy[m]
                       - sin(tau[n])*particles.vx[m])
                       * interv + particles.ex[m])/bx[m]

                eyb =(( - cos(tau[n])*particles.vx[m]
                        - sin(tau[n])*particles.vy[m])
                        * interv + particles.ey[m])/bx[m]

                r[n,1] = cos(tau[n])* exb - sin(tau[n]) * eyb
                r[n,2] = sin(tau[n])* exb + cos(tau[n]) * eyb

            end

            fft!(r,1)
            for n = 2:ntau
                r[n,:] .= -1im * r[n,:]/ltau[n]
            end
            r[1,:] .= 0.0
            ifft!(r,1)

            yt[:,1,m] .= particles.vx[m] .+ (r[:,1] .- r[1,1]) * ep
            yt[:,2,m] .= particles.vy[m] .+ (r[:,2] .- r[1,2]) * ep

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
            Et[1,n,:] .= particles.ex
            Et[2,n,:] .= particles.ey
        end
        @show sum(particles.ix .+ particles.dx)
        @show sum(particles.iy .+ particles.dy)
        @show sum(Et)
        @show sum(Et)

        #  !--time iteration
        #  prediction--
        for m=1:nbpart

            for n=1:ntau

                fx[n,1]   = (   cos(tau[n]) * yt[n,1,m] 
                              + sin(tau[n]) * yt[n,2,m])/bx[m]
                fx[n,2]   = ( - sin(tau[n]) * yt[n,1,m] 
                              + cos(tau[n]) * yt[n,2,m])/bx[m]

                interv    = (1 + 0.5*sin(real(xt[n,1,m]))*sin(real(xt[n,2,m]))-bx[m])/ep

                temp[n,1] = Et[1,n,m]+(  cos(tau[n])*yt[n,2,m]
                                       - sin(tau[n])*yt[n,1,m])*interv
                temp[n,2] = Et[2,n,m]+(- cos(tau[n])*yt[n,1,m]
                                       - sin(tau[n])*yt[n,2,m])*interv

                fy[n,1]   = (cos(tau[n])*temp[n,1]-sin(tau[n])*temp[n,2])/bx[m]
                fy[n,2]   = (sin(tau[n])*temp[n,1]+cos(tau[n])*temp[n,2])/bx[m]

            end

            fft!(fx,1)
            fxtemp0[:,:,m] .= fx/ntau
            fft!(fy,1)
            fytemp0[:,:,m] .= fy/ntau

            tildex[:,:,m] .= fft(xt[:,:,m],1) 

            for n=1:ntau
                temp[n,:] .= (exp.(-1im*ltau[n]*ds[m]/ep) .* tildex[n,:,m]/ntau 
                             .+ pl[n,m] 
                             .* fxtemp0[n,:,m])
            end

            xt[:,1,m] .= bfft(temp[:,1])
            xt[:,2,m] .= bfft(temp[:,2])

            tildey[:,:,m] .= fft(yt[:,:,m],1) 

            for n=1:ntau
                temp[n,:] .= (exp(-1im*ltau[n]*ds[m]/ep) .* tildey[n,:,m]/ntau 
                             .+ pl[n,m]*fytemp0[n,:,m])
            end

            yt[:,1,m] .= bfft(temp[:,1]) 
            yt[:,2,m] .= bfft(temp[:,2]) 

        end
        @show sum(yt)

        for m = 1:nbpart
            time = ds[m]
            @energyuse( )
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

            Et[1,n,:]=particles.ex
            Et[2,n,:]=particles.ey

        end

        @show sum(Et)

        # correction
        for m=1:nbpart

            for n=1:ntau

                fx[n,1]=( cos(tau[n])*yt[n,1,m]+sin(tau[n])*yt[n,2,m])/bx[m]
                fx[n,2]=(-sin(tau[n])*yt[n,1,m]+cos(tau[n])*yt[n,2,m])/bx[m]

                interv=(1 + 0.5*sin(real(xt[n,1,m]))*sin(real(xt[n,2,m]))
                          -bx[m])/ep

                temp[n,1] = Et[1,n,m]+( cos(tau[n])*yt[n,2,m]
                                       -sin(tau[n])*yt[n,1,m])*interv
                temp[n,2] = Et[2,n,m]+(-cos(tau[n])*yt[n,1,m]
                                       -sin(tau[n])*yt[n,2,m])*interv

                fy[n,1] = (cos(tau[n])*temp[n,1]-sin(tau[n])*temp[n,2])/bx[m]
                fy[n,2] = (sin(tau[n])*temp[n,1]+cos(tau[n])*temp[n,2])/bx[m]

            end

            tilde[:,1] .= fft(fx[:,1]) 
            tilde[:,2] .= fft(fx[:,2]) 

            fxtemp1[:,:,m] = tilde/ntau

            tilde[:,1] .= fft(fy[:,1]) 
            tilde[:,2] .= fft(fy[:,2]) 

            fytemp1[:,:,m] = tilde/ntau

            for n=1:ntau
                temp[n,:] = ( exp(-1im*ltau[n]*ds[m]/ep)*tildex[n,:,m]/ntau
                              .+ pl[n,m] * fxtemp0[n,:,m] .+ ql[n,m] * 
                                (fxtemp1[n,:,m] .- fxtemp0[n,:,m]) / ds[m] )
            end

            xt[:,1,m] .= bfft(temp[:,1])
            xt[:,2,m] .= bfft(temp[:,2])

            for n=1:ntau
                temp[n,:] = ( exp(-1im*ltau[n]*ds[m]/ep)*tildey[n,:,m]/ntau 
                             .+ pl[n,m]*fytemp0[n,:,m]
                             .+ ql[n,m]*(fytemp1[n,:,m]-fytemp0[n,:,m])/ds[m])
            end

            yt[:,1,m] .= bfft(temp[:,1]) 
            yt[:,2,m] .= bfft(temp[:,2]) 

        end

        @show sum(yt)

        for m=1:nbpart

            time = ds[m]
            @energyuse()
            @apply_bc( )
            particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
            particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
            particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
            particles.dy[m] = Float32(xxt2/dy - particles.iy[m])

        end

        auxpx[1,:] .= (particles.dx + particles.ix) * dx
        auxpx[2,:].= (particles.dy + particles.iy) * dy

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

            Et[1,n,:] .= particles.ex
            Et[2,n,:] .= particles.ey

        end

        @show sum(Et)

        for m=1:nbpart

            tilde[:,1] .= fft(yt[:,1,m]) 
            tilde[:,2] .= fft(yt[:,2,m]) 

            temp[1,:] .= 0.

            for n=1:ntau
                temp[1,:] .+= tilde[n,:]/ntau * exp(1im*ltau[n]*ds[m]/ep)
            end

            particles.vx[m] = real(cos(ds[m]/ep)*temp[1,1]
                                  +sin(ds[m]/ep)*temp[1,2])
            particles.vy[m] = real(cos(ds[m]/ep)*temp[1,2]
                                  -sin(ds[m]/ep)*temp[1,1])

        end

        @show sum(particles.vx), sum(particles.vy)

    end


    true

end 

const ntau = 16

@time test_pic2d( ntau )
