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

    mesh   :: Mesh
    kx     :: Array{Float64, 2}
    ky     :: Array{Float64, 2}
    ρ̃      :: Array{ComplexF64, 2}
    fft_ex :: Array{ComplexF64, 2}
    fft_ey :: Array{ComplexF64, 2}


    function Poisson( mesh :: Mesh )

        nx = mesh.nx
        ny = mesh.ny

        kx0 = 2π / (mesh.xmax - mesh.xmin)
        ky0 = 2π / (mesh.ymax - mesh.ymin)

        kx = zeros(Float64, (nx÷2+1,ny))
        ky = zeros(Float64, (nx÷2+1,ny))

        for ik=1:nx÷2+1
           kx1 = (ik-1)*kx0
           for jk = 1:ny÷2
              kx[ik,jk] = kx1
              ky[ik,jk] = (jk-1)*ky0
           end
           for jk = ny÷2+1:ny
              kx[ik,jk] = kx1
              ky[ik,jk] = (jk-1-ny)*ky0
           end 
        end

        kx[1,1] = 1.0
        k2  = kx .* kx .+ ky .* ky
        kx .= kx ./ k2
        ky .= ky ./ k2

        ρ̃ = zeros(ComplexF64,(nx÷2+1,ny))

        fft_ex = similar(ρ̃)
        fft_ey = similar(ρ̃)

        new( mesh, kx, ky, ρ̃, fft_ex, fft_ey)

    end	    

end


function ( p :: Poisson )( fields :: MeshFields )

    nx, ny = p.mesh.nx, p.mesh.ny
    dx, dy = p.mesh.dx, p.mesh.dy

    p.ρ̃ .= rfft(fields.ρ[1:nx,1:ny])

    p.fft_ex .= -1im .* p.kx .* p.ρ̃
    p.fft_ey .= -1im .* p.ky .* p.ρ̃

    fields.ex[1:nx,1:ny] .= irfft(p.fft_ex, nx)
    fields.ey[1:nx,1:ny] .= irfft(p.fft_ey, nx)

    fields.ex[nx+1,:] .= fields.ex[1,:]
    fields.ex[:,ny+1] .= fields.ex[:,1]
    fields.ey[nx+1,:] .= fields.ey[1,:]
    fields.ey[:,ny+1] .= fields.ey[:,1]


    sum(fields.ex[1:nx,1:ny] .^2 .+ fields.ey[1:nx,1:ny].^2)*dx*dy

end
