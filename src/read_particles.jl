export read_particles

function read_particles( filename, mesh :: Mesh )

    nbpart = countlines(filename)
    println( " nbpart   : ", nbpart )
    println( " filename : ", filename )

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx   = mesh.xmax - mesh.xmin
    dimy   = mesh.ymax - mesh.ymin

    weight    = (dimx * dimy) / nbpart
    particles = Particles( nbpart, weight)

    open(filename) do f

        for (k,line) in enumerate(eachline(f))
            str_values = split(line)
            particles.ix[k] = parse(Int32,   str_values[1])
            particles.iy[k] = parse(Int32,   str_values[2])
            particles.dx[k] = parse(Float32, str_values[3])
            particles.dy[k] = parse(Float32, str_values[4])
            particles.vx[k] = parse(Float64, str_values[5])
            particles.vy[k] = parse(Float64, str_values[6])
        end

    end

    particles

end

