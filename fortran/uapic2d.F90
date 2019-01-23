program uapic_2d

    use mesh_fields_m
    use particles_m
    use poisson_m
    use interpolation_m
    use compute_rho

    real(8), parameter :: kx       = 0.50
    real(8), parameter :: ky       = 1.0
    real(8), parameter :: pi       = 3.141592654d0
    integer, parameter :: nx       = 128	
    integer, parameter :: ny       = 64 
    real(8), parameter :: tfinal   = 1.0d0 
    real(8), parameter :: eps      = 0.1d0

    real(8) :: t = 0.

    real(8) :: xmin = 0d0
    real(8) :: xmax = dimx
    real(8) :: ymin = 0d0
    real(8) :: ymax = dimy

    type(mesh_t)        :: mesh
    type(mesh_fields_t) :: fields
    type(particles_t)   :: particles
    type(poisson_t)     :: particles

    complex(8), allocatable :: xt(:,:,:)
    complex(8), allocatable :: xf(:,:,:)
    complex(8), allocatable :: yt(:,:,:)
    complex(8), allocatable :: yf(:,:,:)
    complex(8), allocatable :: fx(:,:,:)
    complex(8), allocatable :: fy(:,:,:)
    complex(8), allocatable :: gx(:,:,:)
    complex(8), allocatable :: gy(:,:,:)

    real(8), allocatable :: et(:,:,:)

    real(8) :: pi 
    real(8) :: dimx
    real(8) :: dimy


    pi       = 4d0 * atan(1d0)
    dimx     = 2*pi/kx
    dimy     = 2*pi/ky 

    call init_mesh( mesh, xmin, xmax, nx, ymin, ymax, ny )

    dx = mesh%dx
    dy = mesh%dy

    dt = pi / 2 / (2**3) 
    tfinal = pi / 2

    nstep  = floor(tfinal/dt)

    call init_mesh_fields(fields, mesh )
    
    call init_particles( particles, mesh )

    nbpart = particles%nbpart

    call init_poisson( poisson, mesh )

    call init_ua( ua, ntau, eps, nbpart )

    tau  = ua.tau
    ltau = ua.ltau

    et  = zeros(Float64, (ntau, 2, nbpart))

    allocate(xt(ntau, 2, nbpart))
    allocate(xf(ntau, 2, nbpart))
    allocate(yt(ntau, 2, nbpart))
    allocate(yf(ntau, 2, nbpart))

    ftau = plan_fft(xt,  1)

    allocate(fx(ntau, 2, nbpart))
    allocate(fy(ntau, 2, nbpart))
    allocate(gx(ntau, 2, nbpart))
    allocate(gy(ntau, 2, nbpart))

    call calcul_rho_m6( fields, particles )

    call poisson( fields )

    call interpolate_eb_m6_real( particles, fields )

    do istep = 1,2

        call preparation( ua, dt, particles, xt, yt) 

        call update_particles_e( particles, et, fields, ua, xt)

        call compute_f( fx, fy, ua, particles, xt, yt, et )

        do m = 1, nbpart
            call dfftw_execute_dft( ua%fw, xt(:,1,m), xf(:,1,m))
            call dfftw_execute_dft( ua%fw, xt(:,2,m), xf(:,2,m))
        end do

        call ua_step( xt, xf, ua, particles, fx )

        do m = 1, nbpart
            call dfftw_execute_dft( ua%bw, xt(:,1,m), xt(:,1,m))
            call dfftw_execute_dft( ua%bw, xt(:,2,m), xt(:,2,m))
        end do

        do m = 1, nbpart
            call dfftw_execute_dft( ua%fw, yt(:,1,m), yf(:,1,m))
            call dfftw_execute_dft( ua%fw, yt(:,2,m), yf(:,2,m))
        end do

        call ua_step( yt, yf, ua, particles, fy )

        do m = 1, nbpart
            call dfftw_execute_dft( ua%bw, yt(:,1,m), yt(:,1,m))
            call dfftw_execute_dft( ua%bw, yt(:,2,m), yt(:,2,m))
        end do

        call update_particles_x( particles, fields, ua, xt)

        call poisson( fields ) 

        call update_particles_e( particles, et, fields, ua, xt)

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

        call update_particles_x( particles, fields, ua, xt)

        call nrj = poisson( fields )

        call update_particles_e( particles, et, fields, ua, xt)

        do m=1,nbpart

            t = particles%t(m)

            px = 0d0
            py = 0d0

            do n = 1,ntau
                elt = exp(1im*ltau(n)*t/eps) 
                px = px + yt(n,1,m)/ntau * elt
                py = py + yt(n,2,m)/ntau * elt
            end do

            particles%v(1,m) = real(cos(t/eps)*px+sin(t/eps)*py)
            particles%v(2,m) = real(cos(t/eps)*py-sin(t/eps)*px)

        end do


    end do

    print*, sum(particles%v(1,:)), sum(particles%v(2,:))

end program uapic_2d


