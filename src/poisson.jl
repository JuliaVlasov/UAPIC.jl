using FFTW

export Poisson

"""

   poisson!(ρ, mesh, ex, ey)

Solve the equation Δ Φ = - ρ

 ex = ∂ Φ / ∂ x
 ey = ∂ Φ / ∂ y 

WARNING: the ρ array is destroyed

"""
struct Poisson

    nx :: Int64
    ny :: Int64
    kx :: Vector{Float64}
    ky :: Vector{Float64}

    function Poisson( mesh :: Mesh )

	nx = mesh.nx
	ny = mesh.ny
        kx0 = 2π / (mesh.xmax - mesh.xmin)
        ky0 = 2π / (mesh.ymax - mesh.ymin)

	kx = zeros(Float64, nx)
	ky = zeros(Float64, ny)

        kx .= kx0 * vcat(0:nx÷2-1,-nx÷2:-1)
        kx[1] = 1.0
        ky .= ky0 * vcat(0:ny÷2-1,-ny÷2:-1)
        kx[1] = 1.0

	new( nx, ny, kx, ky)

    end	    

end


function ( p :: Poisson )( ρ  :: Array{Complex{Float64},2}, 
		           ex :: Array{Complex{Float64},2}, 
		           ey :: Array{Complex{Float64},2} )


    fft!(ρ,[1,2])

    for i = 1:p.nx
       kx2 = p.kx[i] * p.kx[i]
       for j =  1:p.ny
          k2 = kx2 + p.ky[j] * p.ky[j]
          ex[i,j] = -1im * p.kx[i]/k2 * ρ[i,j]
	  ey[i,j] = -1im * p.ky[j]/k2 * ρ[i,j]
       end 
    end

    ifft!(ex,[1,2])
    ifft!(ey,[1,2])

end
