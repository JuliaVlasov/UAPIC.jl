!UA scheme for 4d Vlasov in Fluid-scaling with b(x)
!External E
!finite difference solver
!with new IMEX1 and IMEX2
program efd

use, intrinsic :: iso_c_binding
use mesh_fields_m, only: mesh_fields_t, mesh_t, init_mesh, init_mesh_fields
use particles_m, only: init_particles, particles_t
use compute_rho_m, only: compute_rho_m6_real

implicit none
include "fftw3.f03" 

type(mesh_t)        :: mesh
type(mesh_fields_t) :: f
type(particles_t)   :: p

type(C_PTR) :: plan_fw
type(C_PTR) :: plan_bw

integer, parameter :: ntau = 16
integer, parameter :: npp  = 204800

real(8)    :: xxt(2)
real(8)    :: bx
real(8)    :: ds
real(8)    :: interv
real(8)    :: dtau
real(8)    :: tau(ntau)
real(8)    :: ltau(ntau)
real(8)    :: auxpx(2,npp)
real(8)    :: energy
real(8)    :: Et(2,ntau)
real(8)    :: ave(2)
real(8)    :: ave2(2)
real(8)    :: ave3(2)

complex(8) :: pl(       ntau)
complex(8) :: ql(       ntau)
complex(8) :: tildex( 2,ntau)
complex(8) :: tildey( 2,ntau)
complex(8) :: fxtemp1(2,ntau)
complex(8) :: fxtemp0(2,ntau)
complex(8) :: fytemp1(2,ntau)
complex(8) :: fytemp0(2,ntau)
complex(8) :: temp(   2,ntau)
complex(8) :: tilde(  2,ntau)
complex(8) :: h(      2,ntau)
complex(8) :: r(      2,ntau)
complex(8) :: fx(     2,ntau)
complex(8) :: fy(     2,ntau)
complex(8) :: xt(     2,ntau)
complex(8) :: yt(     2,ntau)

real(8) :: time
real(8) :: xmin
real(8) :: xmax
real(8) :: ymin
real(8) :: ymax
integer :: istep
integer :: iplot
integer :: iargc
integer :: n,m
integer :: i
integer :: j
integer :: error

real(8) :: aux1, aux2

character(len=272) :: argv

integer, parameter :: nx = 128
integer, parameter :: ny = 64
integer :: nbpart
integer :: nstep

real(8) :: alpha
real(8) :: dt
real(8) :: dx
real(8) :: dy
real(8) :: ep
real(8) :: kx
real(8) :: ky
real(8) :: dimx
real(8) :: dimy
real(8) :: tfinal
real(8) :: pi
real(8) :: poids
complex(8), parameter :: im = (0d0, 1d0)


pi    = 4d0 * atan(1d0)
alpha = 0.05d0 
kx    = 0.50_8
ky    = 1.0d0   
dimx  = 2*pi/kx
dimy  = 2*pi/ky  
poids = dimx * dimy 

dx = dimx / nx
dy = dimy / ny

dt     = pi/2.d0/(2.0d0**3)!Used to determine N in MRC, not the macro-step
tfinal = pi/2.d0
nstep  = nint(tfinal/dt)

time  = 0.d0

ep   = 0.001d0
dtau = 2.0d0*pi/ntau

plan_fw = fftw_plan_dft_1d(ntau, temp(1,:), tilde(1,:), FFTW_FORWARD, FFTW_ESTIMATE)
plan_bw = fftw_plan_dft_1d(ntau, temp(1,:), tilde(1,:), FFTW_BACKWARD,FFTW_ESTIMATE)

m = ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
do i=1,ntau
  tau(i) = (i-1)*dtau
end do

xmin = 0.0_8; xmax = dimx
ymin = 0.0_8; ymax = dimy

nbpart = npp

call init_mesh( mesh, xmin, xmax, nx, ymin, ymax, ny )
call init_mesh_fields( f, mesh )
call init_particles( p, nbpart, mesh, alpha, kx )

print"('ep = ', g15.3)", ep
print"('dt = ', g15.3)", dt

auxpx = p%x

do m=1,nbpart
    
    time  = 0.d0
    bx    = 1.d0+0.5d0*sin(auxpx(1,m))*sin(auxpx(2,m))
    ds    = dt*bx
    pl(1) = ds
    ql(1) = ds**2/2.0d0
    do i=2,ntau
        pl(i) = ep*im*(exp(-im*ltau(i)*ds/ep)-1.0d0)/ltau(i)
        ql(i) = ep*(ep*(1.0d0-exp(-im*ltau(i)*ds/ep))-im*ltau(i)*ds)/ltau(i)**2
    end do
    xt(1,:)   = auxpx(1,m)
    xt(2,:)   = auxpx(2,m)
    yt(1,:)   = p%v(1,m)
    yt(2,:)   = p%v(2,m)

    !--preparation initial data--
    temp(1,1) = p%v(1,m)/bx
    temp(2,1) = p%v(2,m)/bx
    do n=1,ntau
        h(1,n)=ep*(sin(tau(n))*temp(1,1)-cos(tau(n))*temp(2,1))
        h(2,n)=ep*(sin(tau(n))*temp(2,1)+cos(tau(n))*temp(1,1))
    end do

    xt(1,:)=auxpx(1,m)+h(1,:)-h(1,1)
    xt(2,:)=auxpx(2,m)+h(2,:)-h(2,1)!xt1st
    p%e(1,m)=(0.5d0*cos(auxpx(1,m)/2.d0)*sin(auxpx(2,m)))*(1.d0+0.5d0*sin(time))
    p%e(2,m)=(sin(auxpx(1,m)/2.d0)*cos(auxpx(2,m)))*(1.d0+0.5d0*sin(time))
    do n=1,ntau
        interv=(1.d0+0.5d0*sin(real(xt(1,n)))*sin(real(xt(2,n)))-bx)/bx
        r(1,n) =  interv*p%v(2,m)
        r(2,n) = -interv*p%v(1,m)
    end do
    call fft(r(1,:), tilde(1,:))
    call fft(r(2,:), tilde(2,:))
    ave=tilde(:,1)/ep!dot{Y}_underline^0th
    do n=2,ntau
        tilde(:,n)=-im*tilde(:,n)/ltau(n)
    end do
    tilde(:,1)=0.0d0
    call ifft( tilde(1,:), r(1,:))
    call ifft( tilde(2,:), r(2,:))
    do n=1,ntau
        r(1,n)=ep*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))/bx+r(1,n)
        r(2,n)=ep*(sin(tau(n))*p%e(2,m)-cos(tau(n))*p%e(1,m))/bx+r(2,n)
    end do
    yt(1,:)=p%v(1,m)+(r(1,:)-r(1,1))
    yt(2,:)=p%v(2,m)+(r(2,:)-r(2,1))!yt1st
    !--more preparation--
    do n=1,ntau
        temp(1,n)=ep*(cos(tau(n))*yt(1,n)+sin(tau(n))*yt(2,n))/bx
        temp(2,n)=ep*(cos(tau(n))*yt(2,n)-sin(tau(n))*yt(1,n))/bx
    end do
    call fft(temp(1,:), tilde(1,:))
    call fft(temp(2,:), tilde(2,:))
    do n=2,ntau
        tilde(:,n)=-im*tilde(:,n)/ltau(n)
    end do
    tilde(:,1)=0.0d0
    call ifft(tilde(1,:), h(1,:))
    call ifft(tilde(2,:), h(2,:))
    do n=1,ntau
        h(1,n)=h(1,n)-ep**2/bx*(-cos(tau(n))*ave(1)-sin(tau(n))*ave(2))
        h(2,n)=h(2,n)-ep**2/bx*(-cos(tau(n))*ave(2)+sin(tau(n))*ave(1))!h2nd
    end do
    xt(1,:)=auxpx(1,m)+h(1,:)-h(1,1)
    xt(2,:)=auxpx(2,m)+h(2,:)-h(2,1)!x2nd,residue O(eps^3)
    p%e(1,m)=(0.5d0*cos(auxpx(1,m)/2.d0)*sin(auxpx(2,m)))*0.5d0*cos(time)
    p%e(2,m)=(sin(auxpx(1,m)/2.d0)*cos(auxpx(2,m)))*0.5d0*cos(time)!partial_tE
    do n=1,ntau
        interv=(1.d0+0.5d0*sin(real(xt(1,n)))*sin(real(xt(2,n)))-bx)/bx
        fx(1,n)=interv*ave(2)
        fx(2,n)=-interv*ave(1)
        fy(1,n)=ep/bx*(sin(tau(n))*ave(1)-cos(tau(n))*ave(2))
        fy(2,n)=ep/bx*(cos(tau(n))*ave(1)+sin(tau(n))*ave(2))!partial_sh^1st
        interv=cos(auxpx(1,m))*sin(auxpx(2,m))*fy(1,n)+sin(auxpx(1,m))*cos(auxpx(2,m))*fy(2,n)
        fy(1,n)=interv/bx/2.d0*p%v(2,m)+fx(1,n)
        fy(2,n)=-interv/bx/2.d0*p%v(1,m)+fx(2,n)
        fx(1,n)=ep/bx**2*(-sin(tau(n))*p%e(2,m)+cos(tau(n))*p%e(1,m))
        fx(2,n)=ep/bx**2*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))
    end do
    temp=fy+fx
    call fft(temp(1,:), tilde(1,:))
    call fft(temp(2,:), tilde(2,:))
    do n=2,ntau
        fx(:,n)=-im*tilde(:,n)/ltau(n)
        tilde(:,n)=-tilde(:,n)/ltau(n)**2
    end do
    fx(:,1)    = 0.0d0
    tilde(:,1) = 0.0d0
    call ifft(tilde(1,:), temp(1,:))
    call ifft(tilde(2,:), temp(2,:))!L^-1partial_sr^1st
    r=-ep*temp
    call ifft(fx(1,:), fy(1,:))
    call ifft(fx(2,:), fy(2,:))!partial_sr^1st
    fxtemp0=fy!partial_sr^1st!!!
    do n=1,ntau
        tilde(1,n)=(cos(tau(n))*fy(1,n)+sin(tau(n))*fy(2,n))/bx
        tilde(2,n)=(cos(tau(n))*fy(2,n)-sin(tau(n))*fy(1,n))/bx
    end do
    call fft(tilde(1,:), temp(1,:))
    call fft(tilde(2,:), temp(2,:))
    ave2=temp(:,1)/ntau!ddot{X_underline}!!!!
    do n=1,ntau
        p%e(1,m)=(0.5d0*cos(real(xt(1,n))/2.d0)*sin(real(xt(2,n))))*(1.d0+0.5d0*sin(time))
        p%e(2,m)=(sin(real(xt(1,n))/2.d0)*cos(real(xt(2,n))))*(1.d0+0.5d0*sin(time))
        interv=(1.d0+0.5d0*sin(real(xt(1,n)))*sin(real(xt(2,n)))-bx)/bx
        temp(1,n)=interv*yt(2,n)+ep/bx*(-sin(tau(n))*p%e(2,m)+cos(tau(n))*p%e(1,m))
        temp(2,n)=-interv*yt(1,n)+ep/bx*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))
    end do
    call fft( temp(1,:), tilde(1,:))
    call fft( temp(2,:), tilde(2,:))
    tildex(:,1)=tilde(:,1)/ep!dot{Y_underline}
    do n=2,ntau
        tilde(:,n)=-im*tilde(:,n)/ltau(n)
    end do
    tilde(:,1)=0.0d0
    call ifft(tilde(1,:), temp(1,:))
    call ifft(tilde(2,:), temp(2,:))
    r=r+temp
    yt(1,:)=p%v(1,m)+r(1,:)-r(1,1)
    yt(2,:)=p%v(2,m)+r(2,:)-r(2,1)
    !----end more preparation
    !--even more preparation--
    do n=1,ntau
        temp(1,n)=(cos(tau(n))*r(1,n)+sin(tau(n))*r(2,n))/bx
        temp(2,n)=(cos(tau(n))*r(2,n)-sin(tau(n))*r(1,n))/bx
    end do
    call fft( temp(1,:), tilde(1,:))
    call fft( temp(2,:), tilde(2,:))
    temp(:,1)=tilde(:,1)
    ave3=temp(:,1)!dot{X_underline}!!!
    interv=cos(auxpx(1,m))*sin(auxpx(2,m))*temp(1,1)+sin(auxpx(1,m))*cos(auxpx(2,m))*temp(2,1)
    fx(1,1)=interv/ep/bx*p%v(2,m)/2.d0
    fx(2,1)=-interv/ep/bx*p%v(1,m)/2.d0
    do n=1,ntau
        temp(1,n)=(1.d0+0.5d0*sin(real(xt(1,n)))*sin(real(xt(2,n)))-bx)/bx
    end do
    call fft(temp(1,:), tilde(1,:))
    fx(1,1)=fx(1,1)+tilde(1,1)/ep*ave(2)
    fx(2,1)=fx(2,1)-tilde(1,1)/ep*ave(1)!ddot{Y}_underline^0th!!!!
    do n=1,ntau
        tildey(:,n)=tildex(:,1)+fy(:,n)
        temp(1,n)=(cos(tau(n))*tildey(1,n)+sin(tau(n))*tildey(2,n))
        temp(2,n)=(cos(tau(n))*tildey(2,n)-sin(tau(n))*tildey(1,n))
    end do
    call fft( temp(1,:), tilde(1,:))
    call fft( temp(2,:), tilde(2,:))
    do n=2,ntau
        tilde(:,n)=-im*tilde(:,n)/ltau(n)
    end do
    tilde(:,1)=0.0d0
    call ifft( tilde(1,:), temp(1,:))
    call ifft( tilde(2,:), temp(2,:))
    do n=1,ntau
        fy(1,n)=temp(1,n)*ep/bx-ep**2/bx*(-cos(tau(n))*fx(1,1)-sin(tau(n))*fx(2,1))
        fy(2,n)=temp(2,n)*ep/bx-ep**2/bx*(-cos(tau(n))*fx(2,1)+sin(tau(n))*fx(1,1))!partial_sh2!!!
    end do
    call fft( fy(1,:), tilde(1,:))
    call fft( fy(2,:), tilde(2,:))
    do n=2,ntau
        tilde(:,n)=-im*tilde(:,n)/ltau(n)
    end do
    tilde(:,1)=0.0d0
    call ifft( tilde(1,:), temp(1,:))
    call ifft( tilde(2,:), temp(2,:))
    h=-ep*temp
    do n=1,ntau
        temp(1,n)=ep*(cos(tau(n))*yt(1,n)+sin(tau(n))*yt(2,n))/bx
        temp(2,n)=ep*(cos(tau(n))*yt(2,n)-sin(tau(n))*yt(1,n))/bx
    end do
    call fft( temp(1,:), tilde(1,:))
    call fft( temp(2,:), tilde(2,:))
    do n=2,ntau
        tilde(:,n)=-im*tilde(:,n)/ltau(n)
    end do
    tilde(:,1)=0.0d0
    call ifft( tilde(1,:), temp(1,:))
    call ifft( tilde(2,:), temp(2,:))
    h=h+temp
    xt(1,:)=auxpx(1,m)+h(1,:)-h(1,1)
    xt(2,:)=auxpx(2,m)+h(2,:)-h(2,1)
    
    !--end even more
    !--iteration
    do istep = 1, nstep

        !---imex2 New---
        call compute_fy( xt, yt)
        fytemp0=yt+ds/2.d0*fy
        call fft( fytemp0(1,:), fy(1,:))
        call fft( fytemp0(2,:), fy(2,:))
        do n=1,ntau
            fy(:,n)=fy(:,n)/(1.0d0+im*ds/2.d0*ltau(n)/ep)
        end do
        call ifft( fy(1,:), tildey(1,:))!yt(tn+1/2)
        call ifft( fy(2,:), tildey(2,:))
        do n=1,ntau
            fx(1,n)=(cos(tau(n))*tildey(1,n)+sin(tau(n))*tildey(2,n))/bx
            fx(2,n)=(cos(tau(n))*tildey(2,n)-sin(tau(n))*tildey(1,n))/bx
        end do
        fxtemp0=xt+ds/2.d0*fx
        call fft( fxtemp0(1,:), fx(1,:))
        call fft( fxtemp0(2,:), fx(2,:))
        do n=1,ntau
            fx(:,n)=fx(:,n)/(1.0d0+im*ds/2.d0*ltau(n)/ep)
        end do
        call ifft( fx(1,:), tildex(1,:))!xt(tn+1/2)
        call ifft( fx(2,:), tildex(2,:))
        time=time+dt/2.d0

        call compute_fy( tildex, tildey)

        call fft( fy(1,:), fytemp0(1,:))
        call fft( fy(2,:), fytemp0(2,:))
        call fft( yt(1,:), tildey(1,:))
        call fft( yt(2,:), tildey(2,:))
        do n=1,ntau
            fy(:,n)=(tildey(:,n)*(1.0d0-im*ds/ep/2.0d0*ltau(n)) &
&            +ds*fytemp0(:,n))/(1.0d0+im*ds/2.0d0*ltau(n)/ep)
        end do
        tildey=yt
        call ifft( fy(1,:), yt(1,:))!yt(tn+1)
        call ifft( fy(2,:), yt(2,:))
        tildey=(yt+tildey)/2.d0
        do n=1,ntau
            fx(1,n)=(cos(tau(n))*tildey(1,n)+sin(tau(n))*tildey(2,n))/bx
            fx(2,n)=(cos(tau(n))*tildey(2,n)-sin(tau(n))*tildey(1,n))/bx
        end do
        call fft( fx(1,:), fxtemp0(1,:))
        call fft( fx(2,:), fxtemp0(2,:))
        call fft( xt(1,:), tildex(1,:))
        call fft( xt(2,:), tildex(2,:))
        do n=1,ntau
            fx(:,n)=(tildex(:,n)*(1.0d0-im*ds/ep/2.0d0*ltau(n)) &
        &       +ds*fxtemp0(:,n))/(1.0d0+im*ds/2.0d0*ltau(n)/ep)
        end do
        call ifft( fx(1,:), xt(1,:))!xt(tn+1)
        call ifft( fx(2,:), xt(2,:))
        time=time+dt/2.d0
        !---end imex2 New---

    end do
    call fft( xt(1,:),tilde(1,:))
    call fft( xt(2,:),tilde(2,:))
    temp(:,1)=0.d0
    do n=1,ntau
        temp(:,1)=temp(:,1)+tilde(:,n)*exp(im*ltau(n)*tfinal*bx/ep)
    end do
    xxt=real(temp(:,1))
    call apply_bc()
    p%x(1,m) = xxt(1)
    p%x(2,m) = xxt(2)
    call fft(yt(1,:),tilde(1,:))
    call fft(yt(2,:),tilde(2,:))
    temp(:,1)=0.d0
    do n=1,ntau
        temp(:,1)=temp(:,1)+tilde(:,n)*exp(im*ltau(n)*tfinal*bx/ep)
    end do
    p%v(1,m)=real(cos(tfinal*bx/ep)*temp(1,1)+sin(tfinal*bx/ep)*temp(2,1))
    p%v(2,m)=real(cos(tfinal*bx/ep)*temp(2,1)-sin(tfinal*bx/ep)*temp(1,1))

end do
print*, sum(p%v(1,:))+857.95049281063064, sum(p%v(2,:))+593.40700170710875 


call compute_rho_m6_real(f, p)
stop
open(unit=851,file='UA3ep0001T1.dat')
do i=1,nx
do j=1,ny
write(851,*)i*dx, j*dy, f%rho(i,j)
end do
write(851,*)
end do
close(851)
!call calcul_energy( p, f,energy )
!open(unit=851,file='UA3ep0001T1v.dat')
!do i=1,nx
!do j=1,ny
!write(851,*)i*dx, j*dy, f%rho(i,j)
!end do
!write(851,*)
!end do
!close(851)

contains


subroutine compute_fy( xt, yt )

    complex(8) :: xt(:,:)
    complex(8) :: yt(:,:)

    do n=1,ntau
        Et(1,n)=(0.5d0*cos(real(xt(1,n))/2.d0)*sin(real(xt(2,n))))*(1.d0+0.5d0*sin(time))
        Et(2,n)=(cos(real(xt(2,n)))*sin(real(xt(1,n))/2.d0))*(1.d0+0.5d0*sin(time))
        interv=(1.d0+0.5d0*sin(real(xt(1,n)))*sin(real(xt(2,n)))-bx)/bx/ep
        temp(1,n)=(cos(tau(n))*Et(1,n)-sin(tau(n))*Et(2,n))/bx
        temp(2,n)=(cos(tau(n))*Et(2,n)+sin(tau(n))*Et(1,n))/bx
        fy(1,n)=temp(1,n)+interv*yt(2,n)
        fy(2,n)=temp(2,n)-interv*yt(1,n)
    end do

end subroutine compute_fy

subroutine fft( source, destination)

    complex(8), intent(inout)  :: source(:)
    complex(8), intent(out)    :: destination(:)

    call fftw_execute_dft(plan_fw, source, destination)

    destination = destination / ntau

end subroutine fft

subroutine ifft( source, destination)

    complex(8), intent(inout) :: source(:)
    complex(8), intent(out)   :: destination(:)

    call fftw_execute_dft(plan_bw, source, destination)

end subroutine ifft



subroutine apply_bc()

    do while ( xxt(1) > xmax )
        xxt(1) = xxt(1) - dimx
    end do

    do while ( xxt(1) < xmin )
        xxt(1)= xxt(1) + dimx
    end do

    do while ( xxt(2) > ymax )
        xxt(2)  = xxt(2)  - dimy
    end do

    do while ( xxt(2)  < ymin )
        xxt(2) = xxt(2)  + dimy
    end do

end subroutine apply_bc


end program efd



