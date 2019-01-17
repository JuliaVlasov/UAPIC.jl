export interpol_eb_m6!

function interpol_eb_m6!( particles :: Particles, fields :: MeshFields )

    nx = fields.mesh.nx
    ny = fields.mesh.ny

    for k=1:particles.nbpart
    
       i = particles.ix[k]
       j = particles.iy[k]
    
       im3 = mod(i-3,nx) + 1
       im2 = mod(i-2,nx) + 1
       im1 = mod(i-1,nx) + 1
       ip1 = mod(i+1,nx) + 1
       ip2 = mod(i+2,nx) + 1
       ip3 = mod(i+3,nx) + 1
       jm3 = mod(j-3,ny) + 1
       jm2 = mod(j-2,ny) + 1
       jm1 = mod(j-1,ny) + 1
       jp1 = mod(j+1,ny) + 1
       jp2 = mod(j+2,ny) + 1
       jp3 = mod(j+3,ny) + 1

       i = i + 1
       j = i + 1
    
       dpx = particles.dx[k]
       dpy = particles.dy[k]
    
       cm3x = f_m6(3+dpx)
       cp3x = f_m6(3-dpx)
       cm2x = f_m6(2+dpx)
       cp2x = f_m6(2-dpx)
       cm1x = f_m6(1+dpx)
       cp1x = f_m6(1-dpx)
       cx   = f_m6(  dpx)
       cy   = f_m6(  dpy)
       cm3y = f_m6(3+dpy)
       cp3y = f_m6(3-dpy)
       cm2y = f_m6(2+dpy)
       cp2y = f_m6(2-dpy)
       cm1y = f_m6(1+dpy)
       cp1y = f_m6(1-dpy)
    
       particles.ex[k] =  (
                    + cm3x * cm3y * fields.ex[im3,jm3]   
                    + cm3x * cm2y * fields.ex[im3,jm2]   
                    + cm3x * cm1y * fields.ex[im3,jm1]   
                    + cm3x * cy   * fields.ex[im3,j  ]   
                    + cm3x * cp1y * fields.ex[im3,jp1]   
                    + cm3x * cp2y * fields.ex[im3,jp2]   
                    + cm3x * cp3y * fields.ex[im3,jp3]   
                    + cm2x * cm3y * fields.ex[im2,jm3]   
                    + cm2x * cm2y * fields.ex[im2,jm2]   
                    + cm2x * cm1y * fields.ex[im2,jm1]   
                    + cm2x * cy   * fields.ex[im2,j  ]   
                    + cm2x * cp1y * fields.ex[im2,jp1]   
                    + cm2x * cp2y * fields.ex[im2,jp2]   
                    + cm2x * cp3y * fields.ex[im2,jp3]   
                    + cm1x * cm3y * fields.ex[im1,jm3]   
                    + cm1x * cm2y * fields.ex[im1,jm2]   
                    + cm1x * cm1y * fields.ex[im1,jm1]   
                    + cm1x * cy   * fields.ex[im1,j  ]   
                    + cm1x * cp1y * fields.ex[im1,jp1]   
                    + cm1x * cp2y * fields.ex[im1,jp2]   
                    + cm1x * cp3y * fields.ex[im1,jp3]   
                    + cx   * cm3y * fields.ex[i  ,jm3]   
                    + cx   * cm2y * fields.ex[i  ,jm2]   
                    + cx   * cm1y * fields.ex[i  ,jm1]   
                    + cx   * cy   * fields.ex[i  ,j  ]   
                    + cx   * cp1y * fields.ex[i  ,jp1]   
                    + cx   * cp2y * fields.ex[i  ,jp2]   
                    + cx   * cp3y * fields.ex[i  ,jp3]   
                    + cp1x * cm3y * fields.ex[ip1,jm3]   
                    + cp1x * cm2y * fields.ex[ip1,jm2]   
                    + cp1x * cm1y * fields.ex[ip1,jm1]   
                    + cp1x * cy   * fields.ex[ip1,j  ]   
                    + cp1x * cp1y * fields.ex[ip1,jp1]   
                    + cp1x * cp2y * fields.ex[ip1,jp2]   
                    + cp1x * cp3y * fields.ex[ip1,jp3]   
                    + cp2x * cm3y * fields.ex[ip2,jm3]   
                    + cp2x * cm2y * fields.ex[ip2,jm2]   
                    + cp2x * cm1y * fields.ex[ip2,jm1]   
                    + cp2x * cy   * fields.ex[ip2,j  ]   
                    + cp2x * cp1y * fields.ex[ip2,jp1]   
                    + cp2x * cp2y * fields.ex[ip2,jp2]   
                    + cp2x * cp3y * fields.ex[ip2,jp3]   
                    + cp3x * cm3y * fields.ex[ip3,jm3]   
                    + cp3x * cm2y * fields.ex[ip3,jm2]   
                    + cp3x * cm1y * fields.ex[ip3,jm1]   
                    + cp3x * cy   * fields.ex[ip3,j  ]   
                    + cp3x * cp1y * fields.ex[ip3,jp1]   
                    + cp3x * cp2y * fields.ex[ip3,jp2]   
                    + cp3x * cp3y * fields.ex[ip3,jp3]   )
    
       particles.ey[k] =  (
                    + cm3x * cm3y * fields.ey[im3,jm3]   
                    + cm3x * cm2y * fields.ey[im3,jm2]   
                    + cm3x * cm1y * fields.ey[im3,jm1]   
                    + cm3x * cy   * fields.ey[im3,j  ]   
                    + cm3x * cp1y * fields.ey[im3,jp1]   
                    + cm3x * cp2y * fields.ey[im3,jp2]   
                    + cm3x * cp3y * fields.ey[im3,jp3]   
                    + cm2x * cm3y * fields.ey[im2,jm3]   
                    + cm2x * cm2y * fields.ey[im2,jm2]   
                    + cm2x * cm1y * fields.ey[im2,jm1]   
                    + cm2x * cy   * fields.ey[im2,j  ]   
                    + cm2x * cp1y * fields.ey[im2,jp1]   
                    + cm2x * cp2y * fields.ey[im2,jp2]   
                    + cm2x * cp3y * fields.ey[im2,jp3]   
                    + cm1x * cm3y * fields.ey[im1,jm3]   
                    + cm1x * cm2y * fields.ey[im1,jm2]   
                    + cm1x * cm1y * fields.ey[im1,jm1]   
                    + cm1x * cy   * fields.ey[im1,j  ]   
                    + cm1x * cp1y * fields.ey[im1,jp1]   
                    + cm1x * cp2y * fields.ey[im1,jp2]   
                    + cm1x * cp3y * fields.ey[im1,jp3]   
                    + cx   * cm3y * fields.ey[i  ,jm3]   
                    + cx   * cm2y * fields.ey[i  ,jm2]   
                    + cx   * cm1y * fields.ey[i  ,jm1]   
                    + cx   * cy   * fields.ey[i  ,j  ]   
                    + cx   * cp1y * fields.ey[i  ,jp1]   
                    + cx   * cp2y * fields.ey[i  ,jp2]   
                    + cx   * cp3y * fields.ey[i  ,jp3]   
                    + cp1x * cm3y * fields.ey[ip1,jm3]   
                    + cp1x * cm2y * fields.ey[ip1,jm2]   
                    + cp1x * cm1y * fields.ey[ip1,jm1]   
                    + cp1x * cy   * fields.ey[ip1,j  ]   
                    + cp1x * cp1y * fields.ey[ip1,jp1]   
                    + cp1x * cp2y * fields.ey[ip1,jp2]   
                    + cp1x * cp3y * fields.ey[ip1,jp3]   
                    + cp2x * cm3y * fields.ey[ip2,jm3]   
                    + cp2x * cm2y * fields.ey[ip2,jm2]   
                    + cp2x * cm1y * fields.ey[ip2,jm1]   
                    + cp2x * cy   * fields.ey[ip2,j  ]   
                    + cp2x * cp1y * fields.ey[ip2,jp1]   
                    + cp2x * cp2y * fields.ey[ip2,jp2]   
                    + cp2x * cp3y * fields.ey[ip2,jp3]   
                    + cp3x * cm3y * fields.ey[ip3,jm3]   
                    + cp3x * cm2y * fields.ey[ip3,jm2]   
                    + cp3x * cm1y * fields.ey[ip3,jm1]   
                    + cp3x * cy   * fields.ey[ip3,j  ]   
                    + cp3x * cp1y * fields.ey[ip3,jp1]   
                    + cp3x * cp2y * fields.ey[ip3,jp2]   
                    + cp3x * cp3y * fields.ey[ip3,jp3] ) 
    
    end
    

end 
