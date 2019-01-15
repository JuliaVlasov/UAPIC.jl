@testset "Poisson 2D on rectangular grid" begin

    mesh = Mesh( 0, 2π, 64, 0, 2π, 128)

    fields = MeshFields( mesh )

    x = range(mesh.xmin, stop=mesh.xmax, length=mesh.nx+1) |> collect
    y = range(mesh.ymin, stop=mesh.ymax, length=mesh.ny+1) |> collect
    
    fields.ρ  .= - 2 * sin.(x) .* cos.(y')

    poisson! = Poisson( mesh )

    poisson!( fields )

    err_ex = maximum(abs.( real(fields.ex) .- (cos.(x) .* cos.(y')))) 
    err_ey = maximum(abs.( real(fields.ey) .+ (sin.(x) .* sin.(y')))) 

    println(" errors : $err_ex $err_ey ")

    @test err_ex ≈ 0.0 atol = 1e-15
    @test err_ey ≈ 0.0 atol = 1e-15

end
