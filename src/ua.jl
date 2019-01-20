import FFTW

export UA

struct UA
  
    nτ :: Int64
    ε  :: Float64
    τ  :: Vector{Float64}
    lτ :: Vector{Float64}
    pτ :: FFTW.cFFTWPlan{Complex{Float64},-1,false,1}

    function UA( nτ, ε )

        dτ = 2π / nτ
        
        lτ  = zeros(Float64, nτ)
        lτ .= vcat(0:nτ÷2-1, -nτ÷2:-1) 
        
        τ   = zeros(Float64, nτ)
        τ  .= [ i*dτ for i=0:nτ-1 ]

        pτ  = FFTW.plan_fft(τ)
        new( nτ, ε, τ, lτ, pτ )

    end

end
