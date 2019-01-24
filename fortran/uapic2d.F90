program uapic_2d

    use mesh_fields_m
    use particles_m
    use poisson_m
    use interpolation_m
    use compute_rho_m
    use ua_steps_m

    implicit none


    real(8), parameter :: alpha    = 0.05d0 
    real(8), parameter :: kx       = 0.5d0
    real(8), parameter :: ky       = 1d0
    integer, parameter :: nx       = 128
    integer, parameter :: ny       = 64 
    real(8), parameter :: eps      = 0.1d0
    integer, parameter :: nbpart   = 204800
    integer, parameter :: ntau     = 16
   
    real(8)    :: t = 0d0
    complex(8) :: elt

    real(8) :: tfinal   = 1.0d0 

    type(mesh_t)        :: mesh
    type(mesh_fields_t) :: fields
    type(particles_t)   :: particles
    type(poisson_t)     :: poisson
    type(ua_t)          :: ua

    complex(8), allocatable :: xt(:,:,:)
    complex(8), allocatable :: xf(:,:,:)
    complex(8), allocatable :: yt(:,:,:)
    complex(8), allocatable :: yf(:,:,:)
    complex(8), allocatable :: fx(:,:,:)
    complex(8), allocatable :: fy(:,:,:)
    complex(8), allocatable :: gx(:,:,:)
    complex(8), allocatable :: gy(:,:,:)

    real(8), allocatable :: et(:,:,:)

    real(8) :: dimx
    real(8) :: dimy
    integer :: istep
    integer :: nstep
    integer :: m
    integer :: n

    real(8) :: px
    real(8) :: py
    real(8) :: dx
    real(8) :: dy
    real(8) :: dt

    pi    = 4d0 * atan(1d0)
    dimx  = 2d0*pi/kx
    dimy  = 2d0*pi/ky 

    call init_mesh( mesh, 0d0, dimx, nx, 0d0, dimy, ny )

    dx = mesh%dx
    dy = mesh%dy

    dt = pi / 2d0 / (2d0**3) 
    tfinal = pi / 2d0

    nstep  = floor(tfinal/dt)

    call init_mesh_fields(fields, mesh )
    
    call init_particles( particles, nbpart, alpha, kx, dimx, dimy )

    call init_poisson( poisson, mesh )

    call init_ua( ua, ntau, eps, nbpart )

    allocate(et(ntau, 2, nbpart))

    allocate(xt(ntau, 2, nbpart))
    allocate(xf(ntau, 2, nbpart))
    allocate(yt(ntau, 2, nbpart))
    allocate(yf(ntau, 2, nbpart))

    allocate(fx(ntau, 2, nbpart))
    allocate(fy(ntau, 2, nbpart))
    allocate(gx(ntau, 2, nbpart))
    allocate(gy(ntau, 2, nbpart))

    call compute_rho_m6_real( fields, particles )

    call solve_poisson( poisson, fields )

    call interpolate_eb_m6_real( particles, fields )

    do istep = 1,2

        call preparation( ua, dt, particles, xt, yt) 


        call interpolation( particles, et, fields, ua, xt)

        print*, " xt : ", sum(xt)
        print*, " yt : ", sum(yt)
        print*, " et : ", sum(et)

        call compute_f( fx, fy, ua, particles, xt, yt, et )

        print*, " fx :", sum(fx)
        print*, " fy :", sum(fy)

        print*, " xt :", sum(xt)
        print*, " yt :", sum(yt)

        call ua_step( xt, xf, ua, particles, fx )

        call ua_step( yt, yf, ua, particles, fy )

        print*, " xt :", sum(xt)
        print*, " yt :", sum(yt)
        stop


        call deposition( particles, fields, ua, xt)

        call solve_poisson( poisson, fields ) 

        call interpolation( particles, et, fields, ua, xt)

        call compute_f( gx, gy, ua, particles, xt, yt, et )

        call ua_step( xt, xf, ua, particles, fx ) 

        call ua_step( yt, yf, ua, particles, fy )

        do m=1,nbpart

            t = particles%t(m)

            do n=1,ntau

                xt(n,1,m) = xt(n,1,m) + ua%ql(n,m) * (gx(n,1,m) - fx(n,1,m)) / t
                xt(n,2,m) = xt(n,2,m) + ua%ql(n,m) * (gx(n,2,m) - fx(n,2,m)) / t
                yt(n,1,m) = yt(n,1,m) + ua%ql(n,m) * (gy(n,1,m) - fy(n,1,m)) / t
                yt(n,2,m) = yt(n,2,m) + ua%ql(n,m) * (gy(n,2,m) - fy(n,2,m)) / t

            end do

        end do

        do m = 1, nbpart
            call dfftw_execute_dft( ua%bw, xt(:,1,m), xt(:,1,m))
            call dfftw_execute_dft( ua%bw, xt(:,2,m), xt(:,2,m))
        end do

        xt = xt / real(ntau, kind=8)

        call deposition( particles, fields, ua, xt)

        call solve_poisson( poisson, fields )

        call interpolation( particles, et, fields, ua, xt)

        do m=1,nbpart

            t = particles%t(m)

            px = 0d0
            py = 0d0

            do n = 1,ntau
                elt = exp(cmplx(0d0,1d0,kind=8)*ua%ltau(n)*t/eps) 
                px = px + real(yt(n,1,m)/ntau * elt, kind=8)
                py = py + real(yt(n,2,m)/ntau * elt, kind=8)
            end do

            particles%v(1,m) = cos(t/eps)*px+sin(t/eps)*py
            particles%v(2,m) = cos(t/eps)*py-sin(t/eps)*px

        end do


    end do

    print*, sum(particles%v(1,:)), sum(particles%v(2,:))

end program uapic_2d


