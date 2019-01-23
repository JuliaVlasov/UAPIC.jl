module particles_m

implicit none

type :: particles_t

   real(8), allocatable :: x(:,:)
   real(8), allocatable :: v(:,:)
   real(8), allocatable :: e(:,:)
   real(8), allocatable :: b(:)
   real(8), allocatable :: t(:)
   real(8)              :: w

end type particles_t

contains

subroutine init_particles( self, nbpart, alpha, kx, dimx, dimy )

    type(particles_t)  :: self
    integer            :: nbpart
    real(8)            :: alpha
    real(8)            :: kx
    real(8)            :: dimx
    real(8)            :: dimy

    integer              :: k
    real(8), parameter   :: eps = 1.d-12
    real(8)              :: xi, yi, zi
    real(8)              :: temm
    integer              :: nseed
    integer, allocatable :: seed(:)

    allocate(self%x(2,nbpart))
    allocate(self%v(2,nbpart))
    allocate(self%e(2,nbpart))
    allocate(self%b(nbpart))
    allocate(self%t(nbpart))

    self%w = dimx * dimy / nbpart

    nseed = 33
    allocate(seed(nseed))
    
    seed = [-1584649339, -1457681104, 1579121008,   -819547200, &
              249798090, -517237887,   177452147,   -981503238, &
             1418301473,  1989625004, 2065424384,   -296364178, &
             1658790794, -435188152, -1643185032,   1461389312, &
             1869073641,  1321930686,  483734018,   1269936416, &
            -1999561453,  906251506,   782514880,    428753705, &
            -2031262823,  263953581,  1026600222,  -1118515860, &
             1633712916, -464192498, -1860714528,   1436611533, 0]
    
    call random_seed(put  = seed)


    k = 1
    do while (k<=nbpart)

        call random_number(xi)
        xi   = dimx * xi
        call random_number(yi)
        yi   = dimy * yi
        call random_number(zi)
        zi   = (2d0+alpha)*zi
        temm = 1d0+sin(yi)+alpha*cos(kx*xi)
        if (temm>=zi) then
            self%x(1,k) = xi
            self%x(2,k) = yi
            k = k + 1
        end if

    end do
    
    k = 1
    do while (k<=nbpart)

        call random_number(xi)
        xi = (xi-0.5d0) * 10d0
        call random_number(yi)
        yi = (yi-0.5d0) * 10d0
        call random_number(zi)

        temm = (exp(-((xi-2d0)**2+yi**2)/2d0) &
        &      +exp(-((xi+2d0)**2+yi**2)/2d0))/2d0

        if (temm>=zi) then
            self%v(1,k) = xi
            self%v(2,k) = yi
            k = k + 1
        end if

    end do

end subroutine init_particles

end module particles_m
