using Test
using UAPIC
using FFTW

include("test_poisson_2d.jl")

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

        tilde[1,:] .= fft(xt[1,:,m])
        tilde[2,:] .= fft(xt[2,:,m])

        temp[:,1] .= 0.0
        for n=1:ntau
            temp[:,1] .+= tilde[:,n]/ntau * exp(1im*ltau[n]*time/ep)
        end

        xxt1, xxt2 = real(temp[:,1])

    end)

end 

function test_pic2d( ntau )

    nstepmax = 20000	
    kx       = 0.50
    ky       = 1.0
    dimx     = 2π/kx
    dimy     = 2π/ky 
    nx       = 128	
    ny       = 128 
    tfinal   = 1.0 

    time = 0.
    ep   = 0.1
    dtau = 2π / ntau
    
    m    = ntau÷2
    ltau = vcat(0:m-1, -m:-1)
    
    tau  = [ i*dtau for i=0:ntau-1 ]
    
    xmin, xmax = 0.0, dimx
    ymin, ymax = 0.0, dimy

    mesh = Mesh( xmin, xmax, nx, ymin, ymax, ny )

    dx, dy = mesh.dx, mesh.dy

    dt = π / 2 / (2^3) #Used to determine N in MRC, not the macro-step
    tfinal = π / 2

    nstep  = trunc(Int64, tfinal/dt)

    fields = MeshFields( mesh )
    
    nbpart = 2048

    particles = plasma( mesh, nbpart )


    poisson! = Poisson(mesh)

    calcul_rho_m6!( fields, particles )

    println("∑ |rho| = ", sum(abs.(fields.ρ))*dx*dy)

    x = range(0, stop=dimx, length=nx+1) |> collect
    y = range(0, stop=dimy, length=ny+1) |> collect

    fields.ρ .= (1 .+ 0.05 * cos.(0.5 * x )) .+ sin.(y')

    println( poisson!( fields ) )

    gnuplot("fields.dat", fields) 

    interpol_eb_m6!( particles, fields )

    println(" mesh fields     : ", sum(view(fields.ex,1:nx,1:ny))*mesh.dx*mesh.dy, 
                             "\t", sum(view(fields.ey,1:nx,1:ny))*mesh.dx*mesh.dy)
    println(" particle fields : ", sum(particles.ex) , 
                             "\t", sum(particles.ey) )

    println(" ep   = $ep    ")
    println(" dt   = $dt    ")
    println(" dx   = $dx    ")
    println(" dy   = $dy    ")
    println(" dimx = $dimx  ")
    println(" dimy = $dimy  ")
    println(" npp  = $nbpart")
    println(" ntau = $ntau  ")
    
    auxpx = zeros(Float64, nbpart)
    auxpy = zeros(Float64, nbpart)

    auxpx = (particles.dx+particles.ix) * dx
    auxpy = (particles.dy+particles.iy) * dy

    bx = zeros(Float64, nbpart)
    ds = zeros(Float64, nbpart)
    pl = zeros(ComplexF64, (ntau, nbpart))
    ql = zeros(ComplexF64, (ntau, nbpart))

    tildex = zeros(ComplexF64, (2, ntau, nbpart))
    tildey = zeros(ComplexF64, (2, ntau, nbpart))

    Et = zeros(ComplexF64, (2, ntau, nbpart))

    tilde = zeros(ComplexF64, (2, ntau))
    temp  = zeros(ComplexF64, (2, ntau))
    h     = zeros(ComplexF64, (2, ntau))
    r     = zeros(ComplexF64, (2, ntau))
    fx    = zeros(ComplexF64, (2, ntau))
    fy    = zeros(ComplexF64, (2, ntau))
    xt    = zeros(ComplexF64, (2, ntau, nbpart))
    yt    = zeros(ComplexF64, (2, ntau, nbpart))

    fxtemp1 = zeros(ComplexF64, (2, ntau, nbpart))
    fxtemp0 = zeros(ComplexF64, (2, ntau, nbpart))
    fytemp1 = zeros(ComplexF64, (2, ntau, nbpart))
    fytemp0 = zeros(ComplexF64, (2, ntau, nbpart))

    for istep = 1:nstep

        # preparation
        for m = 1:nbpart

            bx[m]   = 1 + 0.5 * sin(auxpx[m]) * sin(auxpy[m])
            ds[m]   = dt * bx[m]
            pl[1,m] = ds[m]
            ql[1,m] = ds[m]^2 / 2

            for i=2:ntau
                pl[i,m] = ep * 1im*(exp(-1im*ltau[i]*ds[m]/ep)-1)/ltau[i]
                ql[i,m] = ep * (ep*(1-exp(-1im*ltau[i]*ds[m]/ep))
                                -1im*ltau[i]*ds[m])/ltau[i]^2
            end

            # preparation initial data

            temp[1,1] = particles.vx[m]/bx[m]
            temp[2,1] = particles.vy[m]/bx[m]

            for n = 1:ntau
                h[1,n] = ep * (sin(tau[n]) * temp[1,1] - cos(tau[n]) * temp[2,1])
                h[2,n] = ep * (sin(tau[n]) * temp[2,1] + cos(tau[n]) * temp[1,1])
            end

            xt[1,:,m] .= auxpx[m] .+ h[1,:] .- h[1,1]
            xt[2,:,m] .= auxpy[m] .+ h[2,:] .- h[2,1]

            for n = 1:ntau

                local interv :: Float64

                interv=(1+0.5*sin(real(xt[1,n,m]))*sin(real(xt[2,n,m]))-bx[m])/ep

                temp[1,1]=((  cos(tau[n])*particles.vy[m]
                            - sin(tau[n])*particles.vx[m])
                            * interv + particles.ex[m])/bx[m]

                temp[2,1]=(( - cos(tau[n])*particles.vx[m]
                             - sin(tau[n])*particles.vy[m])
                             * interv + particles.ey[m])/bx[m]

                r[1,n] = cos(tau[n])*temp[1,1]-sin(tau[n]) * temp[2,1]
                r[2,n] = sin(tau[n])*temp[1,1]+cos(tau[n]) * temp[2,1]

            end

            tilde[1,:] .= fft(r[1,:])
            tilde[2,:] .= fft(r[2,:])

            for n = 2:ntau
                tilde[:,n] .= -1im * tilde[:,n]/ltau[n]/ntau
            end

            tilde[:,1] .= 0.0

            r[1,:] .= ifft(tilde[1,:])
            r[2,:] .= ifft(tilde[2,:])

            yt[1,:,m] .= particles.vx[m] .+ (r[1,:] .- r[1,1]) * ep
            yt[2,:,m] .= particles.vy[m] .+ (r[2,:] .- r[2,1]) * ep

        end

        for n=1:ntau
            for m=1:nbpart
                xxt1, xxt2 = real(xt[1,n,m]), real(xt[2,n,m])
                @apply_bc()
                particles.ix[m] = trunc(Int32, xxt1/dimx*nx)
                particles.dx[m] = xxt1/dx - particles.ix[m]
                particles.iy[m] = trunc(Int32, xxt2/dimy*ny)
                particles.dy[m] = xxt2/dy - particles.iy[m]
            end
            interpol_eb_m6!( particles, fields )
            Et[1,n,:] .= particles.ex
            Et[2,n,:] .= particles.ey
        end

        #  !--time iteration
        #  prediction--
        for m=1:nbpart

            for n=1:ntau

                fx[1,n]   = ( cos(tau[n]) * yt[1,n,m] + sin(tau[n])*yt[2,n,m])/bx[m]
                fx[2,n]   = (-sin(tau[n]) * yt[1,n,m] + cos(tau[n])*yt[2,n,m])/bx[m]

                interv    = (1 + 0.5*sin(real(xt[1,n,m]))*sin(real(xt[2,n,m]))-bx[m])/ep

                temp[1,n] = Et[1,n,m]+( cos(tau[n])*yt[2,n,m]-sin(tau[n])*yt[1,n,m])*interv
                temp[2,n] = Et[2,n,m]+(-cos(tau[n])*yt[1,n,m]-sin(tau[n])*yt[2,n,m])*interv

                fy[1,n]   = (cos(tau[n])*temp[1,n]-sin(tau[n])*temp[2,n])/bx[m]
                fy[2,n]   = (sin(tau[n])*temp[1,n]+cos(tau[n])*temp[2,n])/bx[m]

            end

            tilde[1,:] .= fft(fx[1,:])
            tilde[2,:] .= fft(fx[2,:])

            fxtemp0[:,:,m] .= tilde/ntau

            tilde[1,:] .= fy[1,:]
            tilde[2,:] .= fy[2,:]

            fytemp0[:,:,m] .= tilde/ntau

            tildex[1,:,m] .= xt[1,:,m] 
            tildex[2,:,m] .= xt[2,:,m] 

            for n=1:ntau
                temp[:,n] = (exp.(-1im*ltau[n]*ds[m]/ep) .* tildex[:,n,m]/ntau 
                             .+ pl[n,m] 
                             .* fxtemp0[:,n,m])
            end

            xt[1,:,m] = ifft(temp[1,:])
            xt[2,:,m] = ifft(temp[2,:])

            tildey[1,:,m] = fft(yt[1,:,m]) 
            tildey[2,:,m] = fft(yt[2,:,m]) 

            for n=1:ntau
                temp[:,n] = (exp(-1im*ltau[n]*ds[m]/ep) .* tildey[:,n,m]/ntau 
                             .+ pl[n,m]*fytemp0[:,n,m])
            end

            yt[1,:,m] = ifft(temp[1,:]) 
            yt[2,:,m] = ifft(temp[2,:]) 

        end

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
        println( poisson!( fields ) )


        for n=1:ntau
            for m=1:nbpart
                xxt1, xxt2 = real(xt[:,n,m])
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

        # correction
        for m=1:nbpart

            for n=1:ntau

                fx[1,n]=( cos(tau[n])*yt[1,n,m]+sin(tau[n])*yt[2,n,m])/bx[m]
                fx[2,n]=(-sin(tau[n])*yt[1,n,m]+cos(tau[n])*yt[2,n,m])/bx[m]

                interv=(1 + 0.5*sin(real(xt[1,n,m]))*sin(real(xt[2,n,m]))-bx[m])/ep

                temp[1,n] = Et[1,n,m]+( cos(tau[n])*yt[2,n,m]-sin(tau[n])*yt[1,n,m])*interv
                temp[2,n] = Et[2,n,m]+(-cos(tau[n])*yt[1,n,m]-sin(tau[n])*yt[2,n,m])*interv

                fy[1,n] = (cos(tau[n])*temp[1,n]-sin(tau[n])*temp[2,n])/bx[m]
                fy[2,n] = (sin(tau[n])*temp[1,n]+cos(tau[n])*temp[2,n])/bx[m]

            end

            tilde[1,:] .= fft(fx[1,:]) 
            tilde[2,:] .= fft(fx[2,:]) 

            fxtemp1[:,:,m] = tilde/ntau

            tilde[1,:] .= fft(fy[1,:]) 
            tilde[2,:] .= fft(fy[2,:]) 

            fytemp1[:,:,m] = tilde/ntau

            for n=1:ntau
                temp[:,n] = ( exp(-1im*ltau[n]*ds[m]/ep)*tildex[:,n,m]/ntau
                              .+ pl[n,m] * fxtemp0[:,n,m] .+ ql[n,m] * 
                                (fxtemp1[:,n,m] .- fxtemp0[:,n,m]) / ds[m] )
            end

            xt[1,:,m] .= ifft(temp[1,:])
            xt[2,:,m] .= ifft(temp[2,:])

            for n=1:ntau
                temp[:,n] = ( exp(-1im*ltau[n]*ds[m]/ep)*tildey[:,n,m]/ntau 
                             .+ pl[n,m]*fytemp0[:,n,m]
                             .+ ql[n,m]*(fytemp1[:,n,m]-fytemp0[:,n,m])/ds[m])
            end

            yt[1,:,m] .= ifft(temp[1,:]) 
            yt[2,:,m] .= ifft(temp[2,:]) 

        end

        for m=1:nbpart

            time = ds[m]
            @energyuse()
            @apply_bc( )
            particles.ix[m] = trunc(Int32,   xxt1/dimx*nx)
            particles.dx[m] = Float32(xxt1/dx - particles.ix[m])
            particles.iy[m] = trunc(Int32,   xxt2/dimy*ny)
            particles.dy[m] = Float32(xxt2/dy - particles.iy[m])

        end

        calcul_rho_m6!( fields, particles )

        @show poisson!( fields )

        for n=1:ntau

            for m=1:nbpart
                xxt1, xxt2 = real(xt[:,n,m])
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

        for m=1:nbpart

            tilde[1,:] .= fft(yt[1,:,m]) 
            tilde[2,:] .= fft(yt[2,:,m]) 

            temp[:,1] .= 0.

            for n=1:ntau
                temp[:,1] .+= tilde[:,n]/ntau * exp(1im*ltau[n]*ds[m]/ep)
            end

            particles.vx[m] = real(cos(ds[m]/ep)*temp[1,1]+sin(ds[m]/ep)*temp[2,1])
            particles.vy[m] = real(cos(ds[m]/ep)*temp[2,1]-sin(ds[m]/ep)*temp[1,1])

        end
    end


    true

end 

const ntau = 16

@time test_pic2d( ntau )
