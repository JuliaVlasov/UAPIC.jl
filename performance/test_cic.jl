function test_cic( nstep, nbpart)

    @show nstep, nbpart

    nx     = 128
    ny     = 128 

    xp  = zeros(Float64, nbpart)
    yp  = zeros(Float64, nbpart)
    vx  = zeros(Float64, (nx+1,ny+1))
    vy  = zeros(Float64, (nx+1,ny+1))

    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0
    dx = (xmax-xmin)/nx
    dy = (ymax-ymin)/ny
    dt = 0.1

    xp .= rand(nbpart)
    yp .= rand(nbpart)

    for j = 1:ny+1, i = 1:nx+1
        vx[i,j] =   (j-1)*dy + xmin
        vy[i,j] = - (i-1)*dx + ymin
    end

    for istep = 1:nstep 

        for k=1:nbpart
        
           xp[k] = xmin + mod(xp[k]-xmin,xmax - xmin)
           yp[k] = ymin + mod(yp[k]-ymin,ymax - ymin)

           xc = (xp[k] - xmin)/dx
           yc = (yp[k] - ymin)/dy

           i = trunc(Int64,xc)
           j = trunc(Int64,yc)

           dpx = xc - i
           dpy = yc - j 

           i += 1
           j += 1
        
           a1 = (1.0-dpx) * (1.0-dpy)
           a2 = (    dpx) * (1.0-dpy)
           a3 = (    dpx) * (    dpy)
           a4 = (1.0-dpx) * (    dpy)
        
           vpx = ( a1 * vx[i  ,j  ] + a2 * vx[i+1,j  ] 
                 + a3 * vx[i+1,j+1] + a4 * vx[i  ,j+1] )
           vpy = ( a1 * vy[i  ,j  ] + a2 * vy[i+1,j  ] 
                 + a3 * vy[i+1,j+1] + a4 * vy[i  ,j+1] )
        
           xp[k] += vpx * dt
           yp[k] += vpy * dt

        end

        #println(sum(xp)/nbpart, " - ", sum(yp)/nbpart)
    
    end
    

end 

test_cic( 1, 1) #trigger comoilation
@time test_cic( 1000, 409600)
