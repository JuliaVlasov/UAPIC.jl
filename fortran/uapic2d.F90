program uapic_2d

    use mesh_fields_m
    use particles_m
    use poisson_m
    use interpolation_m
    use compute_rho_m

    real(8), parameter :: kx       = 0.50
    real(8), parameter :: ky       = 1.0
    real(8), parameter :: pi       = 3.141592654d0
    real(8), parameter :: dimx     = 2*pi/kx
    real(8), parameter :: dimy     = 2*pi/ky 
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

    ua = UA( ntau, ε, nbpart )

    tau  = ua.tau
    ltau = ua.ltau

    et  = zeros(Float64, (ntau, 2, nbpart))

    xt  = zeros(ComplexF64, (ntau, 2, nbpart))
    x̃t  = zeros(ComplexF64, (ntau, 2, nbpart))
    yt  = zeros(ComplexF64, (ntau, 2, nbpart))
    ỹt  = zeros(ComplexF64, (ntau, 2, nbpart))

    ftau = plan_fft(xt,  1)

    fx = zeros(ComplexF64, (ntau, 2, nbpart))
    fy = zeros(ComplexF64, (ntau, 2, nbpart))

    gx = zeros(ComplexF64, (ntau, 2, nbpart))
    gy = zeros(ComplexF64, (ntau, 2, nbpart))

    calcul_rho_m6!( fields, particles )

    nrj =  poisson!( fields )

    interpol_eb_m6!( particles, fields )

    for istep = 1:2

        preparation!( ua, dt, particles, xt, yt) 

        update_particles_e!( particles, et, fields, ua, xt)

        #  prediction

        compute_f!( fx, fy, ua, particles, xt, yt, et )

        mul!(x̃t, ftau, xt) 

        ua_step!( xt, x̃t, ua, particles, fx )

        ifft!(xt,1)

        mul!(ỹt, ftau, yt) 

        ua_step!( yt, ỹt, ua, particles, fy )

        ifft!(yt,1) 

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields ) 

        update_particles_e!( particles, et, fields, ua, xt)

        # correction

        compute_f!( gx, gy, ua, particles, xt, yt, et )

        ua_step!( xt, x̃t, ua, particles, fx ) 

        ua_step!( yt, ỹt, ua, particles, fy )

        for m=1:nbpart

            t = particles.t[m]

            for n=1:ntau

                xt[n,1,m] += ua.ql[n,m] * (gx[n,1,m] - fx[n,1,m]) / t
                xt[n,2,m] += ua.ql[n,m] * (gx[n,2,m] - fx[n,2,m]) / t
                yt[n,1,m] += ua.ql[n,m] * (gy[n,1,m] - fy[n,1,m]) / t
                yt[n,2,m] += ua.ql[n,m] * (gy[n,2,m] - fy[n,2,m]) / t

            end

        end

        ifft!(xt,1)

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields )

        update_particles_e!( particles, et, fields, ua, xt)

        for m=1:nbpart

            t = particles.t[m]

            px, py = 0.0, 0.0
            for n = 1:ntau
                elt = exp(1im*ltau[n]*t/ε) 
                px += yt[n,1,m]/ntau * elt
                py += yt[n,2,m]/ntau * elt
            end

            particles.v[1,m] = real(cos(t/ε)*px+sin(t/ε)*py)
            particles.v[2,m] = real(cos(t/ε)*py-sin(t/ε)*px)

        end


    end

    @show @views sum(particles.v[1,:]), sum(particles.v[2,:])

    true

end program uapic_2d


