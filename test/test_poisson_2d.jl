@testset "Poisson 2D on rectangular grid" begin

    alpha = 0.05
    kx    = 0.5
    ky    = 1.0

    mesh = Mesh( 0, 2π/kx, 64, 0, 2π/ky, 128)

    fields = MeshFields( mesh )
    solutions = MeshFields( mesh)

    x = range(mesh.xmin, stop=mesh.xmax, length=mesh.nx+1) |> collect
    y = range(mesh.ymin, stop=mesh.ymax, length=mesh.ny+1) |> collect
    
    fields.ρ  .= - 8 * sin.(2*x) .* cos.(2*y')

    solutions.ex .=   2 * (cos.(2 .* x) .* cos.(2 .* y'))
    solutions.ey .= - 2 * (sin.(2 .* x) .* sin.(2 .* y'))

    poisson! = Poisson( mesh )

    poisson!( fields )

    err_ex, err_ey = errors( fields, solutions )

    println(" errors : $err_ex $err_ey ")

    @test err_ex ≈ 0.0 atol = 1e-14
    @test err_ey ≈ 0.0 atol = 1e-14

    fields.ρ .= - 4 * ( sin.(2*x) .+ cos.(2*y') )

    poisson!( fields )


    nx, ny = mesh.nx, mesh.ny
    err_x, err_y = 0., 0.

    for j in 1:ny+1, i in 1:nx+1
        solutions.ex[i,j] =   2*cos(2*x[i])
        solutions.ey[i,j] = - 2*sin(2*y[j])
    end

    gnuplot("test2.dat", fields)
    gnuplot("solu2.dat", solutions)

    println(" errors : $err_ex $err_ey ")
    err_ex, err_ey = errors( fields, solutions )


    @test err_ex ≈ 0.0 atol = 1e-14
    @test err_ey ≈ 0.0 atol = 1e-14

end
