using UAPIC
using Test

@testset " Test particles input file read " begin


mesh = Mesh( 0, 4π, 128, 0, 2π, 64 )

particles = read_particles("particles.dat", mesh)

@test particles.nbpart == 204800


end

@testset " Test Particles-MeshFields interaction " begin

     xmin, xmax = 0.0, 20.0
     ymin, ymax = 0.0, 20.0
     nx,   ny   = 20, 20

     mesh = Mesh( xmin, xmax, nx, ymin, ymax, ny )

     dx, dy = mesh.dx, mesh.dy

     fields = MeshFields( mesh )

     nbpart = 121
     p      = 1/nbpart

     particles = Particles( nbpart, p )

     k = 1
     for i = 5:nx-5, j = 5:ny-5

         particles.ix[k] = Int32(i-1)
         particles.iy[k] = Int32(j-1)
         particles.dx[k] = Float32(0.5)
         particles.dy[k] = Float32(0.5)

         k += 1

     end

     calcul_rho_m6!( fields, particles )

     @test integrate( fields.ρ, mesh ) ≈ 0.0 atol = 1e-4

     for i = 1:nx+1, j = 1:ny+1

         fields.ex[i,j] = (i-1)*dx
         fields.ey[i,j] = (j-1)*dy

     end
     
     interpol_eb_m6!( particles, fields )

     err_x, err_y = 0.0, 0.0

     for k = 1:nbpart
         
         ix = particles.ix[k] 
         iy = particles.iy[k] 

         xp = (ix + particles.dx[k]) * dx
         yp = (iy + particles.dy[k]) * dy
         ex = particles.ex[k]
         ey = particles.ey[k]


         err_x += abs( particles.ex[k] - xp ) 
         err_y += abs( particles.ey[k] - yp ) 

     end

     err_x /= nbpart
     err_y /= nbpart

     @test err_x ≈ 0.0 atol = 1e-6
     @test err_y ≈ 0.0 atol = 1e-6

end
