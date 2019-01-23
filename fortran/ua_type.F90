module ua_m

use, intrinsic :: iso_c_binding
use mesh_fields_m
use particles_m

implicit none

include 'fftw3.f03'

type :: ua_t
  
    integer                                :: ntau 
    real(8)                                :: eps 
    real(8), allocatable                   :: tau(:)
    real(8), allocatable                   :: ltau(:)

    complex(8), allocatable                :: pl(:,:)
    complex(8), allocatable                :: ql(:,:)

    type(C_PTR)                            :: fw
    type(C_PTR)                            :: bw
    complex(C_DOUBLE_COMPLEX), allocatable :: ft(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: rt(:,:)
    complex(C_DOUBLE_COMPLEX), allocatable :: rf(:,:)

end type ua_t

complex(8), parameter :: j = (0d0, 1d0)
real(8)               :: pi

integer, private      :: m
integer, private      :: n

contains

subroutine init_ua( ua, ntau, eps, nbpart )

    type(ua_t)          :: ua
    integer, intent(in) :: ntau
    real(8), intent(in) :: eps
    integer, intent(in) :: nbpart

    real(8) :: dtau


    pi = 4d0 * atan(1d0)

    ua%ntau = ntau
    ua%eps  = eps

    dtau = 2d0 * pi / ntau
    
    allocate(ua%ltau(ntau))

    do n = 1, ntau/2
        ua%ltau(n) = n-1
    end do
    do n = ntau/2+1, ntau
        ua%ltau(n) = n - ntau - 1 
    end do
    
    allocate(ua%tau(ntau))

    do n = 1, ntau
        ua%tau(n)  = (n-1) * dtau
    end do

    allocate(ua%pl(ntau, nbpart))
    allocate(ua%ql(ntau, nbpart))

    allocate(ua%ft(ntau))
    allocate(ua%rt(ntau,2))
    allocate(ua%rf(ntau,2))

    ua%fw = fftw_plan_dft_1d(ntau, ua%rt(:,1), ua%rf(:,1), &
    &         FFTW_FORWARD, FFTW_ESTIMATE)
    ua%bw = fftw_plan_dft_1d(ntau, ua%rf(:,1), ua%rt(:,1), &
    &         FFTW_BACKWARD,FFTW_ESTIMATE)

end subroutine init_ua

subroutine preparation( ua, dt, particles, xt, yt) 

    type(ua_t)         , intent(inout) :: ua 
    real(8)            , intent(in)    :: dt
    type(particles_t)  , intent(inout) :: particles 
    complex(8)         , intent(out)   :: xt(:,:,:)         
    complex(8)         , intent(out)   :: yt(:,:,:)         

    real(8)    :: x1
    real(8)    :: x2
    complex(8) :: elt
    real(8)    :: t
    real(8)    :: b
    real(8)    :: h1
    real(8)    :: h2
    real(8)    :: ex
    real(8)    :: ey
    real(8)    :: vx
    real(8)    :: vy
    real(8)    :: vxb
    real(8)    :: vyb
    real(8)    :: xt1
    real(8)    :: xt2
    integer    :: ntau
    integer    :: nbpart
    real(8)    :: eps
    real(8)    :: interv
    real(8)    :: exb
    real(8)    :: eyb

    eps    = ua%eps
    ntau   = ua%ntau
    nbpart = particles%nbpart

    do m = 1,nbpart

        x1 = particles%x(1,m)
        x2 = particles%x(2,m)

        t = dt * particles%b(m)
        b = 1d0 + 0.5d0 * sin(x1) * sin(x2)

        particles%b(m) = b
        particles%t(m) = t

        ua%pl(1,m) = particles%t(m)
        ua%ql(1,m) = particles%t(m)**2 / 2


        do n=2,ua%ntau
            elt = exp(-j*ua%ltau(n)*t/eps) 
            ua%pl(n,m) = eps * j*(elt-1)/ua%ltau(n)
            ua%ql(n,m) = eps * (eps*(1-elt) -j*ua%ltau(n)*t)/ua%ltau(n)**2
        end do

        ex  = particles%e(1,m)
        ey  = particles%e(2,m)
        vx  = particles%v(1,m)
        vy  = particles%v(2,m)
        vxb = vx/b
        vyb = vy/b

        do n = 1,ntau

            h1 = eps * (sin(ua%tau(n)) * vxb - cos(ua%tau(n)) * vyb)
            h2 = eps * (sin(ua%tau(n)) * vyb + cos(ua%tau(n)) * vxb)

            xt1 = x1 + h1 + eps * vyb
            xt2 = x2 + h2 - eps * vxb

            xt(n,1,m) = xt1
            xt(n,2,m) = xt2

            interv=(1+0.5*sin(xt1)*sin(xt2)-b)/eps

            exb = ((  cos(ua%tau(n))*vy - sin(ua%tau(n))*vx) * interv + ex)/b
            eyb = (( -cos(ua%tau(n))*vx - sin(ua%tau(n))*vy) * interv + ey)/b

            ua%rt(n,1) = cos(ua%tau(n))* exb - sin(ua%tau(n)) * eyb
            ua%rt(n,2) = sin(ua%tau(n))* exb + cos(ua%tau(n)) * eyb

        end do

        call fftw_execute_dft(ua%fw, ua%rt(:,1), ua%rf(:,1)) 
        call fftw_execute_dft(ua%fw, ua%rt(:,2), ua%rf(:,2)) 

        do n = 2,ntau
            ua%rf(n,1) = -j * ua%rf(n,1)/ua%ltau(n)
            ua%rf(n,2) = -j * ua%rf(n,2)/ua%ltau(n)
        end do

        call fftw_execute_dft(ua%bw, ua%rf(:,1), ua%rt(:,1)) 
        call fftw_execute_dft(ua%bw, ua%rf(:,2), ua%rt(:,2)) 

        do n = 1,ntau
            yt(n,1,m) = vx + (ua%rt(n,1) - ua%rt(1,1)) * eps
            yt(n,2,m) = vy + (ua%rt(n,2) - ua%rt(1,2)) * eps
        end do

    end do

end subroutine preparation

subroutine update_particles_e( particles, et, fields, ua, xt ) 

    type(particles_t),   intent(out) :: particles
    real(8),             intent(out) :: et(:,:,:)
    type(mesh_fields_t), intent(in)  :: fields
    type(ua_t),          intent(in)  :: ua
    complex(8),          intent(in)  :: xt(:,:,:)

    call interpol_eb_m6_complex( et, fields, xt, particles%nbpart, ua%ntau) 

end subroutine update_particles_e

subroutine update_particles_x( particles, fields, ua, xt) 

    type(particles_t)   :: particles
    type(mesh_fields_t) :: fields
    type(ua_t)          :: ua
    complex(8)          :: xt(:,:,:)

    call calcul_rho_m6_complex( fields, particles, xt, ua )

end

subroutine compute_f( fx, fy, ua, particles, xt, yt, et )

    complex(8)        :: fx(:,:,:)
    complex(8)        :: fy(:,:,:)
    type(ua_t)        :: ua 
    type(particles_t) :: particles
    complex(8)        :: xt(:,:,:)
    complex(8)        :: yt(:,:,:)
    real(8)           :: et(:,:,:)

    real(8)    :: interv
    real(8)    :: tau
    real(8)    :: tmp1
    real(8)    :: tmp2
    real(8)    :: xt1
    real(8)    :: xt2
    complex(8) :: yt1
    complex(8) :: yt2
    real(8)    :: b

    do m=1,particles%nbpart
    
        b = particles%b(m)

        do n=1,ua%ntau
    
            xt1 = real(xt(n,1,m))
            xt2 = real(xt(n,2,m))
    
            yt1 = yt(n,1,m) 
            yt2 = yt(n,2,m) 
    
            tau = ua%tau(n)
    
            fx(n,1,m) = ( cos(tau) * yt1 + sin(tau) * yt2)/b
            fx(n,2,m) = (-sin(tau) * yt1 + cos(tau) * yt2)/b
    
            interv = (1 + 0.5*sin(xt1)*sin(xt2)-b)/ua%eps
    
            tmp1 = et(n,1,m)+(  cos(tau)*yt2 - sin(tau)*yt1)*interv
            tmp2 = et(n,2,m)+(- cos(tau)*yt1 - sin(tau)*yt2)*interv
    
            fy(n,1,m) = (cos(tau)*tmp1-sin(tau)*tmp2)/b
            fy(n,2,m) = (sin(tau)*tmp1+cos(tau)*tmp2)/b
    
        end do

        call fftw_execute_dft(ua%fw, fx(:,1,m), fx(:,1,m)) 
        call fftw_execute_dft(ua%fw, fy(:,1,m), fy(:,1,m)) 
    
    end do


end subroutine compute_f

subroutine ua_step( xt, xf, ua, particles, fx )

    complex(8)        , intent(out) :: xt(:,:,:)
    complex(8)        , intent(in)  :: xf(:,:,:)
    type(ua_t)        , intent(in)  :: ua
    type(particles_t) , intent(in)  :: particles 
    complex(8)        , intent(in)  :: fx(:,:,:)

    real(8)           :: t
    integer           :: m
    integer           :: n
    complex(8)        :: elt

    do m=1,particles%nbpart

        t = particles%t(m)

        do n=1,ua%ntau

            elt = exp(-j*ua%ltau(n)*t/ua%eps) 
            xt(n,1,m) = elt * xf(n,1,m) + ua%pl(n,m) * fx(n,1,m)
            xt(n,2,m) = elt * xf(n,2,m) + ua%pl(n,m) * fx(n,2,m)

        end do

    end do

end subroutine ua_step

end module ua_m
