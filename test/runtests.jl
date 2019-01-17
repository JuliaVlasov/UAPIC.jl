using Test
using UAPIC
using FFTW

include("test_poisson_2d.jl")

"""
UA scheme for 4d VP in Fluid-scaling with b(x)
Update b(x(tn)) every step
"""

function test_pic2d( ntau )

    nstepmax = 20000	
    kx       = 0.50
    ky       = 1.0
    dimx     = 2π/kx
    dimy     = 2π/ky 
    nx       = 128	
    ny       = 128 
    cfl      = 0.9 
    tfinal   = 1.0 
    idiag    = 10  
    bcname   = :period
    exext    = 0.	
    eyext    = 0.	
    bzext    = 0.	
    charge   = 1.0 
    masse    = 1.0  
    c        = 8.
    e0       = 1. 
    relativ  = false

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
    
    particles = plasma( mesh )

    nbpart = particles.nbpart

    poisson! = Poisson(mesh)

    calcul_rho_m6!( fields, particles )

    println("rho       :", sum(abs.(view(fields.ρ, 1:nx, 1:ny))) )

    poisson!( fields )

    interpol_eb_m6!( particles, fields )

    println(" mesh fields     : ", sum(view(fields.ex,1:nx,1:ny))*mesh.dx*mesh.dy, 
                             "\t", sum(view(fields.ey,1:nx,1:ny))*mesh.dx*mesh.dy)
    println(" particle fields : ", sum(particles.ex) , 
                             "\t", sum(particles.ey) )

    
    auxpx = zeros(Float64, nbpart)
    auxpy = zeros(Float64, nbpart)

    auxpx = (particles.dx+particles.ix) * dx
    auxpy = (particles.dy+particles.iy) * dy

    bx = zeros(Float64, nbpart)
    ds = zeros(Float64, nbpart)
    pl = zeros(ComplexF64, (ntau, nbpart))
    ql = zeros(ComplexF64, (ntau, nbpart))
    tilde = zeros(ComplexF64, (2, ntau))
    temp  = zeros(ComplexF64, (2, ntau))
    h     = zeros(ComplexF64, (2, ntau))
    r     = zeros(ComplexF64, (2, ntau))
    xt    = zeros(ComplexF64, (2, ntau, nbpart))
    yt    = zeros(ComplexF64, (2, ntau, nbpart))

    for istep = 1:nstep

        # preparation
        for m = 1:nbpart

            bx[m]   = 1 + 0.5 * sin(auxpx[m]) * sin(auxpy[m])
            ds[m]   = dt * bx[m]
            pl[1,m] = ds[m]
            ql[1,m] = ds[m]^2 / 2

            for i=2:ntau
                pl[i,m] = ep * 1im*(exp(-1im*ltau[i]*ds[m]/ep)-1)/ltau[i]
                ql[i,m] = ep * (ep*(1-exp(-1im*ltau[i]*ds[m]/ep))-1im*ltau[i]*ds[m])/ltau[i]^2
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

                temp[2,1]=(( -cos(tau[n])*particles.vx[m]
                             -sin(tau[n])*particles.vy[m])
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

#    do n=0,ntau-1
#        do m=1,nbpart
#            xxt=dreal(xt(:,n,m))
#            call apply_bc()
#            p%idx[m] = floor(xxt(1)/dimx*nx)
#            p%dpx[m] = real(xxt(1)/dx- p%idx[m], f64)
#            p%idy[m] = floor(xxt(2)/dimy*ny)
#            p%dpy[m] = real(xxt(2)/dy- p%idy[m], f64)
#        end
#        call interpol_eb_m6( f, p )
#        Et(1,n,:)=p%epx
#        Et(2,n,:)=p%epy
#    end

#    !--time iteration---
#    !--prediction--
#    do m=1,nbpart
#        do n=0,ntau-1
#            fx(1,n)=(cos(tau[n])*yt(1,n,m)+sin(tau[n])*yt(2,n,m))/bx[m]
#            fx(2,n)=(-sin(tau[n])*yt(1,n,m)+cos(tau[n])*yt(2,n,m))/bx[m]
#            interv=(1.+0.5*sin(dreal(xt(1,n,m)))*sin(dreal(xt(2,n,m)))-bx[m])/ep
#            temp(1,n)=Et(1,n,m)+(cos(tau[n])*yt(2,n,m)-sin(tau[n])*yt(1,n,m))*interv
#            temp(2,n)=Et(2,n,m)+(-cos(tau[n])*yt(1,n,m)-sin(tau[n])*yt(2,n,m))*interv
#            fy(1,n)=(cos(tau[n])*temp(1,n)-sin(tau[n])*temp(2,n))/bx[m]
#            fy(2,n)=(sin(tau[n])*temp(1,n)+cos(tau[n])*temp(2,n))/bx[m]
#        end
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(2,:), tilde(2,:))
#        fxtemp0(:,:,m)=tilde/ntau!
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(2,:), tilde(2,:))
#        fytemp0(:,:,m)=tilde/ntau
#        call sll_s_fft_exec_c2c_1d(PlnF, xt(1,:,m), tildex(1,:,m))
#        call sll_s_fft_exec_c2c_1d(PlnF, xt(2,:,m), tildex(2,:,m))
#        do n=0,ntau-1
#            temp(:,n)=exp(-1im*ltau[n]*ds[m]/ep)*tildex(:,n,m)/ntau+pl(n,m)*fxtemp0(:,n,m)
#        end
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), xt(1,:,m))!xt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), xt(2,:,m))
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(1,:,m), tildey(1,:,m))
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(2,:,m), tildey(2,:,m))
#        do n=0,ntau-1
#            temp(:,n)=exp(-1im*ltau[n]*ds[m]/ep)*tildey(:,n,m)/ntau+pl(n,m)*fytemp0(:,n,m)
#        end
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), yt(1,:,m))!yt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), yt(2,:,m))
#    end
#    do m=1,nbpart
#        time=ds[m]
#        call energyuse()
#        call apply_bc()
#        p%idx[m] = floor(xxt(1)/dimx*nx)
#        p%dpx[m] = real(xxt(1)/dx- p%idx[m], f64)
#        p%idy[m] = floor(xxt(2)/dimy*ny)
#        p%dpy[m] = real(xxt(2)/dy- p%idy[m], f64)
#    end
#    call calcul_rho_m6( p, f )
#    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
#    do n=0,ntau-1
#        do m=1,nbpart
#            xxt=dreal(xt(:,n,m))
#            call apply_bc()
#            p%idx[m] = floor(xxt(1)/dimx*nx)
#            p%dpx[m] = real(xxt(1)/dx- p%idx[m], f64)
#            p%idy[m] = floor(xxt(2)/dimy*ny)
#            p%dpy[m] = real(xxt(2)/dy- p%idy[m], f64)
#        end
#        call interpol_eb_m6( f, p )
#        Et(1,n,:)=p%epx
#        Et(2,n,:)=p%epy
#    end
#    !-correction-
#    do m=1,nbpart
#        do n=0,ntau-1
#            fx(1,n)=(cos(tau[n])*yt(1,n,m)+sin(tau[n])*yt(2,n,m))/bx[m]
#            fx(2,n)=(-sin(tau[n])*yt(1,n,m)+cos(tau[n])*yt(2,n,m))/bx[m]
#            interv=(1.+0.5*sin(dreal(xt(1,n,m)))*sin(dreal(xt(2,n,m)))-bx[m])/ep
#            temp(1,n)=Et(1,n,m)+(cos(tau[n])*yt(2,n,m)-sin(tau[n])*yt(1,n,m))*interv
#            temp(2,n)=Et(2,n,m)+(-cos(tau[n])*yt(1,n,m)-sin(tau[n])*yt(2,n,m))*interv
#            fy(1,n)=(cos(tau[n])*temp(1,n)-sin(tau[n])*temp(2,n))/bx[m]
#            fy(2,n)=(sin(tau[n])*temp(1,n)+cos(tau[n])*temp(2,n))/bx[m]
#        end
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(2,:), tilde(2,:))
#        fxtemp1(:,:,m)=tilde/ntau!
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(2,:), tilde(2,:))
#        fytemp1(:,:,m)=tilde/ntau
#        do n=0,ntau-1
#            temp(:,n)=exp(-1im*ltau[n]*ds[m]/ep)*tildex(:,n,m)/ntau+pl(n,m)*fxtemp0(:,n,m)+ql(n,m)*(fxtemp1(:,n,m)-fxtemp0(:,n,m))/ds[m]
#        end
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), xt(1,:,m))!xt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), xt(2,:,m))
#        do n=0,ntau-1
#            temp(:,n)=exp(-1im*ltau[n]*ds[m]/ep)*tildey(:,n,m)/ntau+pl(n,m)*fytemp0(:,n,m)+ql(n,m)*(fytemp1(:,n,m)-fytemp0(:,n,m))/ds[m]
#        end
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), yt(1,:,m))!yt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), yt(2,:,m))
#    end
#    do m=1,nbpart
#        time=ds[m]
#        call energyuse()
#        call apply_bc()
#        p%idx[m] = floor(xxt(1)/dimx*nx)
#        p%dpx[m] = real(xxt(1)/dx- p%idx[m], f64)
#        p%idy[m] = floor(xxt(2)/dimy*ny)
#        p%dpy[m] = real(xxt(2)/dy- p%idy[m], f64)
#    end
#    auxpx(:)=(p%dpx+p%idx)*dx
#    auxpy(:)=(p%dpy+p%idy)*dy
#    call calcul_rho_m6( p, f )
#    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
#    do n=0,ntau-1
#        do m=1,nbpart
#            xxt=dreal(xt(:,n,m))
#            call apply_bc()
#            p%idx[m] = floor(xxt(1)/dimx*nx)
#            p%dpx[m] = real(xxt(1)/dx- p%idx[m], f64)
#            p%idy[m] = floor(xxt(2)/dimy*ny)
#            p%dpy[m] = real(xxt(2)/dy- p%idy[m], f64)
#        end
#        call interpol_eb_m6( f, p )
#        Et(1,n,:)=p%epx
#        Et(2,n,:)=p%epy
#    end
#    do m=1,nbpart
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(1,:,m),tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(2,:,m),tilde(2,:))
#        temp(:,1)=0.
#        do n=0,ntau-1
#            temp(:,1)=temp(:,1)+tilde(:,n)/ntau*cdexp(1im*ltau[n]*ds[m]/ep)
#        end
#        p%vpx[m]=dreal(cos(ds[m]/ep)*temp(1,1)+sin(ds[m]/ep)*temp(2,1))
#        p%vpy[m]=dreal(cos(ds[m]/ep)*temp(2,1)-sin(ds[m]/ep)*temp(1,1))
#    end
    end

#open(unit=851,file='T5.dat')
#do i=1,nx
#do j=1,ny
#write(851,*)f%r0(i,j)
#end
#end
#close(851)
#
#
#contains
#
#subroutine apply_bc()
#do while ( xxt(1) > xmax )
#xxt(1) = xxt(1) - dimx
#end
#do while ( xxt(1) < xmin )
#xxt(1)= xxt(1) + dimx
#end
#do while ( xxt(2) > ymax )
#xxt(2)  = xxt(2)  - dimy
#end
#do while ( xxt(2)  < ymin )
#xxt(2) = xxt(2)  + dimy
#end
#end subroutine apply_bc
#
#subroutine energyuse()
#call sll_s_fft_exec_c2c_1d(PlnF, xt(1,:,m),tilde(1,:))
#call sll_s_fft_exec_c2c_1d(PlnF, xt(2,:,m),tilde(2,:))
#temp(:,1)=0.
#do n=0,ntau-1
#temp(:,1)=temp(:,1)+tilde(:,n)/ntau*cdexp(1im*ltau[n]*time/ep)
#end
#xxt=dreal(temp(:,1))
#end subroutine

    true

end 

const ntau = 16

@test test_pic2d( ntau )
