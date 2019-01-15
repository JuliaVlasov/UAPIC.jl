"""
 M6 function 
 Quintic spline
               |  (3-q)^5-6(2-q)^5+15(1-q)^5  for 0 <= q < 1
 M6(x) = 1/120 |  (3-q)^5-6(2-q)^5            for 1 <= q < 2
               |  (3-q)^5                     for 2 <= q < 3
               |  0                           for q >= 3
"""

function f_m6( q :: Float32 )

    f_m6 = 0.0

    if ( q < 1.0 ) 
        f_m6 = (3-q)^5-6*(2-q)^5+15*(1-q)^5
    elseif ( q >= 1.0 && q < 2.0 )
        f_m6 = (3-q)^5-6*(2-q)^5
    elseif ( q >= 2 && q < 3 )
        f_m6 = (3-q)^5
    end
    
    f_m6 / 120

end 

export calcul_rho_m6!

function calcul_rho_m6!( fields  :: MeshFields, particles :: Particles) 

    fill!(fields.ρ , 0.0)
    nx = fields.mesh.nx
    ny = fields.mesh.ny

    for k = 1:particles.nbpart
    
	    i      = particles.ix[k]
        j      = particles.iy[k]
        dpx    = particles.dx[k]
        dpy    = particles.dy[k]
        weight = particles.p
      
        im3 = mod(i-3,nx)+1
        im2 = mod(i-2,nx)+1
        im1 = mod(i-1,nx)+1
        ip1 = mod(i+1,nx)+1
        ip2 = mod(i+2,nx)+1
        ip3 = mod(i+3,nx)+1
        jm3 = mod(j-3,ny)+1
        jm2 = mod(j-2,ny)+1
        jm1 = mod(j-1,ny)+1
        jp1 = mod(j+1,ny)+1
        jp2 = mod(j+2,ny)+1
        jp3 = mod(j+3,ny)+1

        i = i+1
        j = j+1
      
        cm3x = f_m6(3+dpx)
        cp3x = f_m6(3-dpx)
        cm2x = f_m6(2+dpx)
        cp2x = f_m6(2-dpx)
        cm1x = f_m6(1+dpx)
        cp1x = f_m6(1-dpx)
        cx   = f_m6(dpx)
        cy   = f_m6(dpy)
        cp1y = f_m6(1-dpy)
        cm1y = f_m6(1+dpy)
        cp2y = f_m6(2-dpy)
        cm2y = f_m6(2+dpy)
        cp3y = f_m6(3-dpy)
        cm3y = f_m6(3+dpy)
      
	    fields.ρ[im3,jm3] += cm3x * cm3y * weight
        fields.ρ[im3,jm2] += cm3x * cm2y * weight
        fields.ρ[im3,jm1] += cm3x * cm1y * weight
        fields.ρ[im3,j  ] += cm3x * cy   * weight
        fields.ρ[im3,jp1] += cm3x * cp1y * weight
        fields.ρ[im3,jp2] += cm3x * cp2y * weight
        fields.ρ[im3,jp3] += cm3x * cp3y * weight
      
	    fields.ρ[im2,jm3] += cm2x * cm3y * weight
        fields.ρ[im2,jm2] += cm2x * cm2y * weight
        fields.ρ[im2,jm1] += cm2x * cm1y * weight
        fields.ρ[im2,j  ] += cm2x * cy   * weight
        fields.ρ[im2,jp1] += cm2x * cp1y * weight
        fields.ρ[im2,jp2] += cm2x * cp2y * weight
        fields.ρ[im2,jp3] += cm2x * cp3y * weight
      
	    fields.ρ[im1,jm3] += cm1x * cm3y * weight
        fields.ρ[im1,jm2] += cm1x * cm2y * weight
        fields.ρ[im1,jm1] += cm1x * cm1y * weight
        fields.ρ[im1,j  ] += cm1x * cy   * weight
        fields.ρ[im1,jp1] += cm1x * cp1y * weight
        fields.ρ[im1,jp2] += cm1x * cp2y * weight
        fields.ρ[im1,jp3] += cm1x * cp3y * weight
      
	    fields.ρ[i  ,jm3] += cx   * cm3y * weight
        fields.ρ[i  ,jm2] += cx   * cm2y * weight
        fields.ρ[i  ,jm1] += cx   * cm1y * weight
        fields.ρ[i  ,j  ] += cx   * cy   * weight
        fields.ρ[i  ,jp1] += cx   * cp1y * weight
        fields.ρ[i  ,jp2] += cx   * cp2y * weight
        fields.ρ[i  ,jp3] += cx   * cp3y * weight
      
	    fields.ρ[ip1,jm3] += cp1x * cm3y * weight
        fields.ρ[ip1,jm2] += cp1x * cm2y * weight
        fields.ρ[ip1,jm1] += cp1x * cm1y * weight
        fields.ρ[ip1,j  ] += cp1x * cy   * weight
        fields.ρ[ip1,jp1] += cp1x * cp1y * weight
        fields.ρ[ip1,jp2] += cp1x * cp2y * weight
        fields.ρ[ip1,jp3] += cp1x * cp3y * weight
      
	    fields.ρ[ip2,jm3] += cp2x * cm3y * weight
        fields.ρ[ip2,jm2] += cp2x * cm2y * weight
        fields.ρ[ip2,jm1] += cp2x * cm1y * weight
        fields.ρ[ip2,j  ] += cp2x * cy   * weight
        fields.ρ[ip2,jp1] += cp2x * cp1y * weight
        fields.ρ[ip2,jp2] += cp2x * cp2y * weight
        fields.ρ[ip2,jp3] += cp2x * cp3y * weight
      
	    fields.ρ[ip3,jm3] += cp3x * cm3y * weight
        fields.ρ[ip3,jm2] += cp3x * cm2y * weight
        fields.ρ[ip3,jm1] += cp3x * cm1y * weight
        fields.ρ[ip3,j  ] += cp3x * cy   * weight
        fields.ρ[ip3,jp1] += cp3x * cp1y * weight
        fields.ρ[ip3,jp2] += cp3x * cp2y * weight
        fields.ρ[ip3,jp3] += cp3x * cp3y * weight

    end
    

    dimx = fields.mesh.xmax - fields.mesh.xmin
    dimy = fields.mesh.ymax - fields.mesh.ymin

    dx = fields.mesh.dx
    dy = fields.mesh.dy

    fields.ρ[1:nx,ny+1] .= view(fields.ρ,1:nx,1)
    fields.ρ[nx+1,1:ny] .= view(fields.ρ,1,1:ny)
    fields.ρ[nx+1,ny+1]  = fields.ρ[1,1]
    
    fields.ρ .= fields.ρ / (dx*dy)
    
    rho_total = sum(view(fields.ρ,1:nx,1:ny)) * dx * dy

    println( " rho_total = $rho_total ")
    
    fields.ρ .= fields.ρ .- rho_total/dimx/dimy

end 

