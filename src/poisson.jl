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

    ρ̃  :: Array{ComplexF64, 2}

    fft_ex :: Array{ComplexF64, 2}
    fft_ey :: Array{ComplexF64, 2}


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

        ρ̃ = zeros(ComplexF64,(nx,ny))
        fft_ex = similar(ρ̃)
        fft_ey = similar(ρ̃)

	    new( nx, ny, kx, ky, ρ̃, fft_ex, fft_ey)

    end	    

end


function ( p :: Poisson )( fields :: MeshFields )

    nx, ny = p.nx, p.ny

    p.ρ̃ .= fields.ρ[1:nx,1:ny]

    fft!(p.ρ̃,[1,2])

    for i = 1:nx
       kx2 = p.kx[i] * p.kx[i]
       for j =  1:ny
          k2 = kx2 + p.ky[j] * p.ky[j]
          p.fft_ex[i,j] = -1im * p.kx[i]/k2 * p.ρ̃[i,j]
	      p.fft_ey[i,j] = -1im * p.ky[j]/k2 * p.ρ̃[i,j]
       end 
    end

    ifft!(p.fft_ex,[1,2])
    ifft!(p.fft_ey,[1,2])

    fields.ex[1:nx,1:ny] .= real(p.fft_ex)
    fields.ey[1:nx,1:ny] .= real(p.fft_ey)

    fields.ex[nx+1,:] = view(fields.ex,1,:)
    fields.ex[:,ny+1] = view(fields.ex,:,1)
    fields.ey[nx+1,:] = view(fields.ey,1,:)
    fields.ey[:,ny+1] = view(fields.ey,:,1)
    

end
