module cic_m

implicit none

contains

subroutine push( xp, yp, vx, vy, nbpart, mesh )

    real(8), intent(inout) :: xp(nbpart)
    real(8), intent(inout) :: yp(nbpart)
    real(8), intent(in)    :: vx(:,:)
    real(8), intent(in)    :: vy(:,:)
    integer, intent(in)    :: nbpart
    type(mesh_t)           :: mesh

    ! i,j+1_________i+1,j+1
    !  |     |        |
    !  | a2     a1    |
    !  |____ P _______|
    !  |              |
    !  | a3  |  a4    |
    !  |     |        |
    ! i,j____|______i+1,j
    
    do k=1,nbpart
    
       xp(k) = xmin + modulo(xp(k)-xmin, xmax-xmin)
       yp(k) = ymin + modulo(yp(k)-xmin, ymax-ymin)
    
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
    
end subroutine push

end module cic_m
