@testset "Poisson 2D on rectangular grid" begin

    alpha = 0.05
    kx    = 0.5
    ky    = 1.0

    mesh = Mesh( 0, 2π/kx, 64, 0, 2π/ky, 128)

    fields = MeshFields( mesh )

    x = range(mesh.xmin, stop=mesh.xmax, length=mesh.nx+1) |> collect
    y = range(mesh.ymin, stop=mesh.ymax, length=mesh.ny+1) |> collect
    
    fields.ρ  .= - 8 * sin.(2*x) .* cos.(2*y')

    poisson! = Poisson( mesh )

    poisson!( fields )

    err_ex = maximum(abs.( fields.ex .- 2 * (cos.(2 .* x) .* cos.(2 .* y')))) 
    err_ey = maximum(abs.( fields.ey .+ 2 * (sin.(2 .* x) .* sin.(2 .* y')))) 

    println(" errors : $err_ex $err_ey ")

    @test err_ex ≈ 0.0 atol = 1e-14
    @test err_ey ≈ 0.0 atol = 1e-14

    nx    = mesh.nx
    ny    = mesh.ny
    dx    = mesh.dx
    dy    = mesh.dy
    for i=1:nx+1
	aux2 = alpha * cos(kx*(i-1)*dx)
        for j=1:ny+1
	    fields.ρ[i,j] = aux2+sin(ky*(j-1)*dy)
        end
    end

    poisson!( fields )

    err_x = 0.0
    err_y = 0.0
    for i=1:nx+1
	aux1 = alpha/kx * sin(kx*(i-1)*dx)
        for j=1:ny+1
	    err_x += abs(fields.ex[i,j] - aux1)
	    err_y += abs(fields.ey[i,j] - cos(ky*(j-1)*dy))
        end
    end

    println( err_x, " ", err_y )

end
