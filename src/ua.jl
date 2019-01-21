import FFTW

export UA

struct UA
  
    ntau :: Int64
    ε    :: Float64
    tau  :: Vector{Float64}
    ltau :: Vector{Float64}
    ftau :: Vector{Float64}
    ptau :: FFTW.cFFTWPlan{ComplexF64,-1,false,1}

    function UA( ntau, ε )

        dtau = 2π / ntau
        
        ltau  = zeros(Float64, ntau)
        ltau .= vcat(0:ntau÷2-1, -ntau÷2:-1) 
        
        tau   = zeros(Float64, ntau)
        tau  .= [ i*dtau for i=0:ntau-1 ]

        ftau  = zeros(ComplexF64, ntau)

        ptau  = FFTW.plan_fft(ftau)

        new( ntau, ε, tau, ltau, ftau, ptau )

    end

end
