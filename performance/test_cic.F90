program test_cic

implicit none

integer, parameter :: nstep  = 1000
integer, parameter :: nbpart = 409600
integer, parameter :: nx     = 128
integer, parameter :: ny     = 128

integer    :: i, j, k
integer    :: istep, n
real(8)    :: dt, pi
real(8)    :: xmin, xmax, ymin, ymax
real(8)    :: dpx, dpy, dx , dy, xn, yn
real(8)    :: a1, a2, a3, a4
real(8)    :: vpx, vpy

real(8), allocatable :: xp(:), yp(:)
real(8), allocatable :: vx(:,:), vy(:,:)
character(len=9) :: filename

allocate(xp(nbpart))
allocate(yp(nbpart))
allocate(vx(nx+1,ny+1))
allocate(vy(nx+1,ny+1))

xmin = -5d0; xmax = 5d0
ymin = -5d0; ymax = 5d0
dx = (xmax-xmin)/real(nx,kind=8)
dy = (ymax-ymin)/real(ny,kind=8)
dt = 0.1d0

call random_number(xp)
call random_number(yp)

do j = 1, ny+1
    do i = 1, nx+1
        vx(i,j) =   real(j-1)*dy - 5.0
        vy(i,j) = - real(i-1)*dx + 5.0
    end do
end do

print*, nstep, nbpart

do istep = 1, nstep

    ! i,j+1_________i+1,j+1
    !  |     |        |
    !  | a2     a1    |
    !  |____ P _______|
    !  |              |
    !  | a3  |  a4    |
    !  |     |        |
    ! i,j____|______i+1,j
    
    do k=1,nbpart
    
       xp(k) = xmin + modulo(xp(k), xmax-xmin)
       yp(k) = ymin + modulo(yp(k), ymax-ymin)
       xn  = (xp(k)-xmin)/dx
       yn  = (yp(k)-ymin)/dy
       i   = floor(xn)
       j   = floor(yn)
       dpx = xn - real(i, kind=8)
       dpy = yn - real(j, kind=8)
       i   = i + 1
       j   = j + 1
    
       a1 = (1.0-dpx) * (1.0-dpy)
       a2 = (    dpx) * (1.0-dpy)
       a3 = (    dpx) * (    dpy)
       a4 = (1.0-dpx) * (    dpy)
    
       vpx = a1 * vx(i  ,j  ) + a2 * vx(i+1,j  ) &
         & + a3 * vx(i+1,j+1) + a4 * vx(i  ,j+1)
       vpy = a1 * vy(i  ,j  ) + a2 * vy(i+1,j  ) &
         & + a3 * vy(i+1,j+1) + a4 * vy(i  ,j+1)

       xp(k) = xp(k) + vpx * dt
       yp(k) = yp(k) + vpy * dt

    end do

end do

print*, sum(xp)/nbpart, sum(yp)/nbpart

end program test_cic
