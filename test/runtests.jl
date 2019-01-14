using Test
using UAPIC

include("test_poisson_2d.jl")

"""
UA scheme for 4d VP in Fluid-scaling with b(x)
Update b(x(tn)) every step
"""

function test_pic2d( ntau )

    nstepmax = 20000	
    dimx     = 20.
    dimy     = 20.      
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

    f = MeshFields( nx, ny )
    
    time = 0.
    ep   = 0.1
    dtau = 2π / ntau
    
    m    = ntau÷2
    ltau = vcat(0:m-1, -m:-1)
    
    tau  = [ i*dtau for i=0:ntau-1 ]
    
    xmin, xmax = 0.0, dimx
    ymin, ymax = 0.0, dimy

    mesh = Mesh( xmin, xmax, nx, ymin, ymax, ny )
    
    p = plasma( mesh )

    poisson = Poisson(mesh)

    calcul_rho_m6!( f, p )

#call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
#call interpol_eb_m6( f, p )
#
#auxpx(1,:)=(p%dpx+p%idx)*dx
#auxpx(2,:)=(p%dpy+p%idy)*dy
#
#!--iteration
#do istep = 1, nstep
#    !--preparation--
#    do m=1,nbpart
#        bx(m)=1.d0+0.5d0*dsin(auxpx(1,m))*dsin(auxpx(2,m))
#        ds(m)=dt*bx(m)
#        pl(0,m)=ds(m)
#        ql(0,m)=ds(m)**2/2.0d0
#        do i=1,ntau-1
#            pl(i,m) =ep*sll_p_i1*(exp(-sll_p_i1*ltau(i)*ds(m)/ep)-1.0d0)/ltau(i)
#            ql(i,m)=ep*(ep*(1.0d0-exp(-sll_p_i1*ltau(i)*ds(m)/ep))-sll_p_i1*ltau(i)*ds(m))/ltau(i)**2
#        enddo
#        !--preparation initial data--
#        temp(1,1)=p%vpx(m)/bx(m)
#        temp(2,1)=p%vpy(m)/bx(m)
#        do n=0,ntau-1
#            h(1,n)=ep*(dsin(tau(n))*temp(1,1)-dcos(tau(n))*temp(2,1))
#            h(2,n)=ep*(dsin(tau(n))*temp(2,1)+dcos(tau(n))*temp(1,1))
#        enddo
#        xt(1,:,m)=auxpx(1,m)+h(1,:)-h(1,0)
#        xt(2,:,m)=auxpx(2,m)+h(2,:)-h(2,0)
#        do n=0,ntau-1
#            interv=(1.d0+0.5d0*dsin(dreal(xt(1,n,m)))*dsin(dreal(xt(2,n,m)))-bx(m))/ep
#            temp(1,1)=((cos(tau(n))*p%vpy(m)-sin(tau(n))*p%vpx(m))*interv+p%epx(m))/bx(m)
#            temp(2,1)=((-cos(tau(n))*p%vpx(m)-sin(tau(n))*p%vpy(m))*interv+p%epy(m))/bx(m)
#            r(1,n)=dcos(tau(n))*temp(1,1)-dsin(tau(n))*temp(2,1)
#            r(2,n)=dsin(tau(n))*temp(1,1)+dcos(tau(n))*temp(2,1)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnF, r(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, r(2,:), tilde(2,:))
#        do n=1,ntau-1
#            tilde(:,n)=-sll_p_i1*tilde(:,n)/ltau(n)/ntau
#        enddo
#        tilde(:,0)=0.0d0
#        call sll_s_fft_exec_c2c_1d(PlnB, tilde(1,:), r(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnB, tilde(2,:), r(2,:))
#        yt(1,:,m)=p%vpx(m)+(r(1,:)-r(1,0))*ep
#        yt(2,:,m)=p%vpy(m)+(r(2,:)-r(2,0))*ep
#    enddo
#    do n=0,ntau-1
#        do m=1,nbpart
#            xxt=dreal(xt(:,n,m))
#            call apply_bc()
#            p%idx(m) = floor(xxt(1)/dimx*nx)
#            p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
#            p%idy(m) = floor(xxt(2)/dimy*ny)
#            p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
#        enddo
#        call interpol_eb_m6( f, p )
#        Et(1,n,:)=p%epx
#        Et(2,n,:)=p%epy
#    enddo
#    !--time iteration---
#    !--prediction--
#    do m=1,nbpart
#        do n=0,ntau-1
#            fx(1,n)=(dcos(tau(n))*yt(1,n,m)+dsin(tau(n))*yt(2,n,m))/bx(m)
#            fx(2,n)=(-dsin(tau(n))*yt(1,n,m)+dcos(tau(n))*yt(2,n,m))/bx(m)
#            interv=(1.d0+0.5d0*dsin(dreal(xt(1,n,m)))*dsin(dreal(xt(2,n,m)))-bx(m))/ep
#            temp(1,n)=Et(1,n,m)+(cos(tau(n))*yt(2,n,m)-sin(tau(n))*yt(1,n,m))*interv
#            temp(2,n)=Et(2,n,m)+(-cos(tau(n))*yt(1,n,m)-sin(tau(n))*yt(2,n,m))*interv
#            fy(1,n)=(dcos(tau(n))*temp(1,n)-dsin(tau(n))*temp(2,n))/bx(m)
#            fy(2,n)=(dsin(tau(n))*temp(1,n)+dcos(tau(n))*temp(2,n))/bx(m)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(2,:), tilde(2,:))
#        fxtemp0(:,:,m)=tilde/ntau!
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(2,:), tilde(2,:))
#        fytemp0(:,:,m)=tilde/ntau
#        call sll_s_fft_exec_c2c_1d(PlnF, xt(1,:,m), tildex(1,:,m))
#        call sll_s_fft_exec_c2c_1d(PlnF, xt(2,:,m), tildex(2,:,m))
#        do n=0,ntau-1
#            temp(:,n)=exp(-sll_p_i1*ltau(n)*ds(m)/ep)*tildex(:,n,m)/ntau+pl(n,m)*fxtemp0(:,n,m)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), xt(1,:,m))!xt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), xt(2,:,m))
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(1,:,m), tildey(1,:,m))
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(2,:,m), tildey(2,:,m))
#        do n=0,ntau-1
#            temp(:,n)=exp(-sll_p_i1*ltau(n)*ds(m)/ep)*tildey(:,n,m)/ntau+pl(n,m)*fytemp0(:,n,m)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), yt(1,:,m))!yt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), yt(2,:,m))
#    enddo
#    do m=1,nbpart
#        time=ds(m)
#        call energyuse()
#        call apply_bc()
#        p%idx(m) = floor(xxt(1)/dimx*nx)
#        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
#        p%idy(m) = floor(xxt(2)/dimy*ny)
#        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
#    enddo
#    call calcul_rho_m6( p, f )
#    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
#    do n=0,ntau-1
#        do m=1,nbpart
#            xxt=dreal(xt(:,n,m))
#            call apply_bc()
#            p%idx(m) = floor(xxt(1)/dimx*nx)
#            p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
#            p%idy(m) = floor(xxt(2)/dimy*ny)
#            p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
#        enddo
#        call interpol_eb_m6( f, p )
#        Et(1,n,:)=p%epx
#        Et(2,n,:)=p%epy
#    enddo
#    !-correction-
#    do m=1,nbpart
#        do n=0,ntau-1
#            fx(1,n)=(dcos(tau(n))*yt(1,n,m)+dsin(tau(n))*yt(2,n,m))/bx(m)
#            fx(2,n)=(-dsin(tau(n))*yt(1,n,m)+dcos(tau(n))*yt(2,n,m))/bx(m)
#            interv=(1.d0+0.5d0*dsin(dreal(xt(1,n,m)))*dsin(dreal(xt(2,n,m)))-bx(m))/ep
#            temp(1,n)=Et(1,n,m)+(cos(tau(n))*yt(2,n,m)-sin(tau(n))*yt(1,n,m))*interv
#            temp(2,n)=Et(2,n,m)+(-cos(tau(n))*yt(1,n,m)-sin(tau(n))*yt(2,n,m))*interv
#            fy(1,n)=(dcos(tau(n))*temp(1,n)-dsin(tau(n))*temp(2,n))/bx(m)
#            fy(2,n)=(dsin(tau(n))*temp(1,n)+dcos(tau(n))*temp(2,n))/bx(m)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fx(2,:), tilde(2,:))
#        fxtemp1(:,:,m)=tilde/ntau!
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(1,:), tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, fy(2,:), tilde(2,:))
#        fytemp1(:,:,m)=tilde/ntau
#        do n=0,ntau-1
#            temp(:,n)=exp(-sll_p_i1*ltau(n)*ds(m)/ep)*tildex(:,n,m)/ntau+pl(n,m)*fxtemp0(:,n,m)+ql(n,m)*(fxtemp1(:,n,m)-fxtemp0(:,n,m))/ds(m)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), xt(1,:,m))!xt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), xt(2,:,m))
#        do n=0,ntau-1
#            temp(:,n)=exp(-sll_p_i1*ltau(n)*ds(m)/ep)*tildey(:,n,m)/ntau+pl(n,m)*fytemp0(:,n,m)+ql(n,m)*(fytemp1(:,n,m)-fytemp0(:,n,m))/ds(m)
#        enddo
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), yt(1,:,m))!yt(t1)
#        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), yt(2,:,m))
#    enddo
#    do m=1,nbpart
#        time=ds(m)
#        call energyuse()
#        call apply_bc()
#        p%idx(m) = floor(xxt(1)/dimx*nx)
#        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
#        p%idy(m) = floor(xxt(2)/dimy*ny)
#        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
#    enddo
#    auxpx(1,:)=(p%dpx+p%idx)*dx
#    auxpx(2,:)=(p%dpy+p%idy)*dy
#    call calcul_rho_m6( p, f )
#    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
#    do n=0,ntau-1
#        do m=1,nbpart
#            xxt=dreal(xt(:,n,m))
#            call apply_bc()
#            p%idx(m) = floor(xxt(1)/dimx*nx)
#            p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
#            p%idy(m) = floor(xxt(2)/dimy*ny)
#            p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
#        enddo
#        call interpol_eb_m6( f, p )
#        Et(1,n,:)=p%epx
#        Et(2,n,:)=p%epy
#    enddo
#    do m=1,nbpart
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(1,:,m),tilde(1,:))
#        call sll_s_fft_exec_c2c_1d(PlnF, yt(2,:,m),tilde(2,:))
#        temp(:,1)=0.d0
#        do n=0,ntau-1
#            temp(:,1)=temp(:,1)+tilde(:,n)/ntau*cdexp(sll_p_i1*ltau(n)*ds(m)/ep)
#        enddo
#        p%vpx(m)=dreal(dcos(ds(m)/ep)*temp(1,1)+dsin(ds(m)/ep)*temp(2,1))
#        p%vpy(m)=dreal(dcos(ds(m)/ep)*temp(2,1)-dsin(ds(m)/ep)*temp(1,1))
#    enddo
#enddo
#!--only needed for rhov--
#!do m=1,nbpart
#!time=ds(m)
#!call energyuse()
#!call apply_bc()
#!p%idx(m) = floor(xxt(1)/dimx*nx)
#!p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
#!p%idy(m) = floor(xxt(2)/dimy*ny)
#!p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
#!enddo
#!call calcul_energy( p, f, energy)
#!----
#open(unit=851,file='T5.dat')
#do i=1,nx
#do j=1,ny
#write(851,*)f%r0(i,j)
#enddo
#enddo
#close(851)
#
#call sll_s_fft_free(PlnF)
#call sll_s_fft_free(PlnB)
#
#
#contains
#
#subroutine apply_bc()
#do while ( xxt(1) > xmax )
#xxt(1) = xxt(1) - dimx
#enddo
#do while ( xxt(1) < xmin )
#xxt(1)= xxt(1) + dimx
#enddo
#do while ( xxt(2) > ymax )
#xxt(2)  = xxt(2)  - dimy
#enddo
#do while ( xxt(2)  < ymin )
#xxt(2) = xxt(2)  + dimy
#enddo
#end subroutine apply_bc
#
#subroutine energyuse()
#call sll_s_fft_exec_c2c_1d(PlnF, xt(1,:,m),tilde(1,:))
#call sll_s_fft_exec_c2c_1d(PlnF, xt(2,:,m),tilde(2,:))
#temp(:,1)=0.d0
#do n=0,ntau-1
#temp(:,1)=temp(:,1)+tilde(:,n)/ntau*cdexp(sll_p_i1*ltau(n)*time/ep)
#enddo
#xxt=dreal(temp(:,1))
#end subroutine

    true

end 

const ntau = 16

@test test_pic2d( ntau )
