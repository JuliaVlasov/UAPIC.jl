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

    complex(8) :: px
    complex(8) :: py
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

    print*, " xp : ", sum(particles%x(1,:))
    print*, " yp : ", sum(particles%x(2,:))

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

    do istep = 1,nstep

        call preparation( ua, dt, particles, xt, yt) 

        print*, " xt : ", sum(xt)
        print*, " yt : ", sum(yt)

        call interpolation( particles, et, fields, ua, xt)

        print*, " et : ", sum(et)

        call compute_f( fx, fy, ua, particles, xt, yt, et )

        print*, " fx :", sum(fx)
        print*, " fy :", sum(fy)

        call ua_step1( xt, xf, ua, particles, fx )
        call ua_step1( yt, yf, ua, particles, fy )
        print*, " xt :", sum(xt)
        print*, " yt :", sum(yt)
        print*, " xf :", sum(xf)
        print*, " yf :", sum(yf)

        call deposition( particles, fields, ua, xt)

        call solve_poisson( poisson, fields ) 
        print*," nrj = ", sum(fields%e(1,:,:)**2+fields%e(2,:,:)**2)*dx*dy


        print*, " ex :", sum(fields%e(1,:,:))
        print*, " ey :", sum(fields%e(2,:,:))

        call interpolation( particles, et, fields, ua, xt)
        
        print*, " et : ", sum(et)

        call compute_f( gx, gy, ua, particles, xt, yt, et )

        print*, " gx :", sum(gx)
        print*, " gy :", sum(gy)

        call ua_step2( xt, xf, ua, particles, fx, gx ) 
        call ua_step2( yt, yf, ua, particles, fy, gy )

        print*, " xt :", sum(xt)
        print*, " yt :", sum(yt)

        call deposition( particles, fields, ua, xt)

        call solve_poisson( poisson, fields )

        call interpolation( particles, et, fields, ua, xt)

        print*, " et : ", sum(et)

        do m=1,nbpart

            t = particles%t(m)

            call fftw_execute_dft(ua%fw, yt(:,1,m), yf(:,1,m)) 
            call fftw_execute_dft(ua%fw, yt(:,2,m), yf(:,2,m)) 

            px = (0d0, 0d0)
            py = (0d0, 0d0)

            do n = 1,ntau
                elt = exp(cmplx(0d0,1d0,kind=8)*ua%ltau(n)*t/eps) 
                px = px + yf(n,1,m)/real(ntau,kind=8) * elt
                py = py + yf(n,2,m)/real(ntau,kind=8) * elt
            end do

            particles%v(1,m) = real(cos(t/eps)*px+sin(t/eps)*py)
            particles%v(2,m) = real(cos(t/eps)*py-sin(t/eps)*px)

        end do
        print*, " vx : ", sum(particles%v(1,:))
        print*, " vy : ", sum(particles%v(2,:))


    end do


end program uapic_2d


