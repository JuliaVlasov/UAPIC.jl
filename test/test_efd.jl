"""
UA scheme for 4d Vlasov in Fluid-scaling with b(x)
External E
finite difference solver
with new IMEX1 and IMEX2

"""
function efd( ntau, nbpart )

tau   = zeros(Float64, ntau)
ltau  = zeros(Float64, ntau)
et    = zeros(Float64, (2,ntau))

pl    = zeros(ComplexF64,  ntau)
ql    = zeros(ComplexF64,  ntau)
xf    = zeros(ComplexF64, (ntau,2))
yf    = zeros(ComplexF64, (ntau,2))
gx    = zeros(ComplexF64, (ntau,2))
gy    = zeros(ComplexF64, (ntau,2))
temp  = zeros(ComplexF64, (ntau,2))
tilde = zeros(ComplexF64, (ntau,2))
h     = zeros(ComplexF64, (ntau,2))
r     = zeros(ComplexF64, (ntau,2))
fx    = zeros(ComplexF64, (ntau,2))
fy    = zeros(ComplexF64, (ntau,2))
xt    = zeros(ComplexF64, (ntau,2))
yt    = zeros(ComplexF64, (ntau,2))

nx = 128
ny = 64

alpha = 0.05 
kx    = 0.50
ky    = 1.0   
dimx  = 2pi/kx
dimy  = 2pi/ky  
poids = dimx * dimy 
dx    = dimx / nx
dy    = dimy / ny

dt     = pi/2./(2.0^3)
tfinal = pi/2.
nstep  = trunc(Int64,tfinal/dt)
time   = 0.
ep     = 0.001
dtau   = 2pi/ntau

call init_fft( ntau, temp[:,1), tilde[:,1) )

m    = ntau÷2
ltau = [0:m-1; -m:-1]
for i=1:ntau
  tau[i) = [i-1)*dtau
end

xmin = 0.0_8; xmax = dimx
ymin = 0.0_8; ymax = dimy

mesh = init_mesh( xmin, xmax, nx, ymin, ymax, ny )
f    = init_mesh_fields( mesh )
p    = init_particles( nbpart, mesh, alpha, kx )

for m=1:nbpart
    
    x1 = p.x[1,m]
    x2 = p.x[2,m]
    v1 = p.v[1,m]
    v2 = p.v[2,m]

    time  = 0.
    bx    = 1+0.5*sin(x1)*sin(x2)
    ds    = dt*bx

    pl[1] = ds
    ql[1] = ds^2/2.0

    for i=2,ntau
        pl[i) = ep*im*(exp(-im*ltau[i)*ds/ep)-1.0)/ltau[i)
        ql[i) = ep*(ep*(1.0-exp(-im*ltau[i)*ds/ep))-im*ltau[i)*ds)/ltau[i)^2
    end


    xt[:,1] = x1
    xt[:,2] = x2
    yt[:,1] = v1
    yt[:,2] = v2

    #--preparation initial data--
    temp[1,1) = v1/bx
    temp[1,2) = v2/bx
    for n=1:ntau
        h[n,1)=ep*(sin(tau[n))*temp[1,1)-cos(tau[n))*temp(1,2))
        h[n,2)=ep*(sin(tau[n))*temp[1,2)+cos(tau[n))*temp[1,1))
    end

    xt[:,1)=x1+h[:,1)-h[1,1)
    xt[:,2)=x2+h[:,2)-h[1,2)!xt1st

    p%e[1,m)=(0.5*cos(x1/2.)*sin(x2))*(1.+0.5*sin(time))
    p%e[2,m)=(sin(x1/2.)*cos(x2))*(1.+0.5*sin(time))
    for n=1:ntau
        interv=(1.+0.5*sin(real(xt[n,1)))*sin(real(xt[n,2)))-bx)/bx
        r[n,1) =  interv*v2
        r[n,2) = -interv*v1
    end

    call fft(r, tilde)
    ave = tilde[1,:)/ep!dot{Y}_underline^0th
    for n=2,ntau
        tilde[n,:)=-im*tilde[n,:)/ltau[n)
    end

    tilde[1,:)=0.0
    call ifft( tilde, r)

    for n=1:ntau
        r[n,1)=ep*(sin(tau[n))*p%e[1,m)+cos(tau[n))*p%e[2,m))/bx+r[n,1)
        r[n,2)=ep*(sin(tau[n))*p%e[2,m)-cos(tau[n))*p%e[1,m))/bx+r[n,2)
    end

    yt[:,1)=v1+(r[:,1)-r[1,1))
    yt[:,2)=v2+(r[:,2)-r[1,2))!yt1st

    #--more preparation--

    for n=1:ntau
        temp[n,1)=ep*(cos(tau[n))*yt[n,1)+sin(tau[n))*yt[n,2))/bx
        temp[n,2)=ep*(cos(tau[n))*yt[n,2)-sin(tau[n))*yt[n,1))/bx
    end

    call fft(temp, tilde)

    for n=2,ntau
        tilde[n,:)=-im*tilde[n,:)/ltau[n)
    end

    tilde[1,:)=0.0
    call ifft(tilde, h)

    for n=1:ntau
        h[n,1)=h[n,1)-ep^2/bx*[-cos(tau[n))*ave[1)-sin(tau[n))*ave[2))
        h[n,2)=h[n,2)-ep^2/bx*(-cos(tau[n))*ave[2)+sin(tau(n))*ave(1))!h2nd
    end

    xt(:,1)=x1+h(:,1)-h(1,1)
    xt(:,2)=x2+h(:,2)-h(1,2)!x2nd,residue O(eps^3)

    p%e(1,m)=(0.5*cos(x1/2.)*sin(x2))*0.5*cos(time)
    p%e(2,m)=(sin(x1/2.)*cos(x2))*0.5*cos(time)!partial_tE

    for n=1:ntau

        interv=(1.+0.5*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx

        fx(n,1)=interv*ave(2)
        fx(n,2)=-interv*ave(1)

        fy(n,1)=ep/bx*(sin(tau(n))*ave(1)-cos(tau(n))*ave(2))
        fy(n,2)=ep/bx*(cos(tau(n))*ave(1)+sin(tau(n))*ave(2))!partial_sh^1st

        interv=cos(x1)*sin(x2)*fy(n,1)+sin(x1)*cos(x2)*fy(n,2)

        fy(n,1)=interv/bx/2.*v2+fx(n,1)
        fy(n,2)=-interv/bx/2.*v1+fx(n,2)

        fx(n,1)=ep/bx^2*(-sin(tau(n))*p%e(2,m)+cos(tau(n))*p%e(1,m))
        fx(n,2)=ep/bx^2*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))

    end

    temp=fy+fx

    call fft(temp, tilde)

    for n=2,ntau
        fx(n,:)=-im*tilde(n,:)/ltau(n)
        tilde(n,:)=-tilde(n,:)/ltau(n)^2
    end

    fx(1,:)    = 0.0
    tilde(1,:) = 0.0

    call ifft(tilde, temp)

    r = - ep * temp

    call ifft(fx, fy)

    gx=fy
    for n=1:ntau
        tilde(n,1)=(cos(tau(n))*fy(n,1)+sin(tau(n))*fy(n,2))/bx
        tilde(n,2)=(cos(tau(n))*fy(n,2)-sin(tau(n))*fy(n,1))/bx
    end

    call fft(tilde, temp)

    ave2=temp(1,:)/ntau
    for n=1:ntau

        p%e(1,m)=(0.5*cos(real(xt(n,1))/2.)*sin(real(xt(n,2))))*(1.+0.5*sin(time))
        p%e(2,m)=(sin(real(xt(n,1))/2.)*cos(real(xt(n,2))))*(1.+0.5*sin(time))

        interv=(1.+0.5*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx

        temp(n,1)=interv*yt(n,2)+ep/bx*(-sin(tau(n))*p%e(2,m)+cos(tau(n))*p%e(1,m))
        temp(n,2)=-interv*yt(n,1)+ep/bx*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))
    end

    call fft( temp, tilde)

    xf(1,:)=tilde(1,:)/ep!dot{Y_underline}
    for n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end
    tilde(1,:)=0.0
    call ifft(tilde, temp)
    r=r+temp

    yt(:,1)=v1+r(:,1)-r(1,1)
    yt(:,2)=v2+r(:,2)-r(1,2)

    #----end more preparation
    #--even more preparation--

    for n=1:ntau
        temp(n,1)=(cos(tau(n))*r(n,1)+sin(tau(n))*r(n,2))/bx
        temp(n,2)=(cos(tau(n))*r(n,2)-sin(tau(n))*r(n,1))/bx
    end

    call fft( temp, tilde)

    temp(1,:)=tilde(1,:)
    ave3=temp(1,:)!dot{X_underline}!!!
    interv=cos(x1)*sin(x2)*temp(1,1)+sin(x1)*cos(x2)*temp(1,2)

    fx(1,1)=interv/ep/bx*v2/2.
    fx(1,2)=-interv/ep/bx*v1/2.

    for n=1:ntau
        temp(n,1)=(1.+0.5*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx
    end

    call fft(temp(:,1), tilde(:,1))

    fx(1,1)=fx(1,1)+tilde(1,1)/ep*ave(2)
    fx(1,2)=fx(1,2)-tilde(1,1)/ep*ave(1)!ddot{Y}_underline^0th!!!!

    for n=1:ntau
        yf(n,:)=xf(1,:)+fy(n,:)
        temp(n,1)=(cos(tau(n))*yf(n,1)+sin(tau(n))*yf(n,2))
        temp(n,2)=(cos(tau(n))*yf(n,2)-sin(tau(n))*yf(n,1))
    end

    call fft( temp, tilde)

    for n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end

    tilde(1,:)=0.0

    call ifft( tilde, temp)

    for n=1:ntau
        fy(n,1)=temp(n,1)*ep/bx-ep^2/bx*(-cos(tau(n))*fx(1,1)-sin(tau(n))*fx(1,2))
        fy(n,2)=temp(n,2)*ep/bx-ep^2/bx*(-cos(tau(n))*fx(1,2)+sin(tau(n))*fx(1,1))!partial_sh2!!!
    end

    call fft( fy, tilde)

    for n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end
    tilde(1,:)=0.0

    call ifft( tilde, temp)

    h=-ep*temp

    for n=1:ntau
        temp(n,1)=ep*(cos(tau(n))*yt(n,1)+sin(tau(n))*yt(n,2))/bx
        temp(n,2)=ep*(cos(tau(n))*yt(n,2)-sin(tau(n))*yt(n,1))/bx
    end

    call fft( temp, tilde)

    for n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end

    tilde(1,:)=0.0

    call ifft( tilde, temp)

    h = h + temp

    xt(:,1)=x1+h(:,1)-h(1,1)
    xt(:,2)=x2+h(:,2)-h(1,2)
    
    #--end even more
    #--iteration
    for istep = 1, nstep

        #---imex2 New---
        call compute_fy( xt, yt)

        gy=yt+ds/2.*fy

        call fft( gy, fy)

        for n=1:ntau
            fy(n,:)=fy(n,:)/(1.0+im*ds/2.*ltau(n)/ep)
        end

        call ifft( fy, yf)!yt(tn+1/2)

        for n=1:ntau
            fx(n,1)=(cos(tau(n))*yf(n,1)+sin(tau(n))*yf(n,2))/bx
            fx(n,2)=(cos(tau(n))*yf(n,2)-sin(tau(n))*yf(n,1))/bx
        end

        gx = xt + ds/2. * fx

        call fft( gx, fx)

        for n=1:ntau
            fx(n,:)=fx(n,:)/(1.0+im*ds/2.*ltau(n)/ep)
        end

        call ifft( fx, xf)!xt(tn+1/2)

        time = time + dt/2.

        call compute_fy( xf, yf)

        call fft( fy, gy)
        call fft( yt, yf)

        for n=1:ntau
            fy(n,:)=(yf(n,:)*(1.0-im*ds/ep/2.0*ltau(n)) &
&            +ds*gy(n,:))/(1.0+im*ds/2.0*ltau(n)/ep)
        end

        yf=yt

        call ifft( fy, yt)!yt(tn+1)

        yf=(yt+yf)/2.
        for n=1:ntau
            fx(n,1)=(cos(tau(n))*yf(n,1)+sin(tau(n))*yf(n,2))/bx
            fx(n,2)=(cos(tau(n))*yf(n,2)-sin(tau(n))*yf(n,1))/bx
        end

        call fft( fx, gx)
        call fft( xt, xf)

        for n=1:ntau
            fx(n,:)=(xf(n,:)*(1.0-im*ds/ep/2.0*ltau(n)) &
        &       +ds*gx(n,:))/(1.0+im*ds/2.0*ltau(n)/ep)
        end

        call ifft( fx, xt)!xt(tn+1)

        time=time+dt/2.

        #---end imex2 New---

    end

    call fft( xt,tilde)

    temp(1,:)=0.
    for n=1:ntau
        temp(1,:)=temp(1,:)+tilde(n,:)*exp(im*ltau(n)*tfinal*bx/ep)
    end

    xxt=real(temp(1,:))

    call apply_bc()

    p%x(1,m) = xxt(1)
    p%x(2,m) = xxt(2)

    call fft(yt,tilde)

    temp(1,:)=0.
    for n=1:ntau
        temp(1,:)=temp(1,:)+tilde(n,:)*exp(im*ltau(n)*tfinal*bx/ep)
    end

    p%v(1,m)=real(cos(tfinal*bx/ep)*temp(1,1)+sin(tfinal*bx/ep)*temp(1,2))
    p%v(2,m)=real(cos(tfinal*bx/ep)*temp(1,2)-sin(tfinal*bx/ep)*temp(1,1))

end
print*, sum(p%v(1,:))+857.95049281063064, sum(p%v(2,:))+593.40700170710875 


compute_rho_m6_real!(f, p)

#open(unit=851,file='UA3ep0001T1.dat')
#for i=1,nx
#for j=1,ny
#write(851,*)i*dx, j*dy, f%rho(i,j)
#end
#write(851,*)
#end
#close(851)

function compute_fy( xt, yt )

    complex(8) :: xt(:,:)
    complex(8) :: yt(:,:)

    for n=1:ntau
        et(1,n)=(0.5*cos(real(xt(n,1))/2.)*sin(real(xt(n,2))))*(1.+0.5*sin(time))
        et(2,n)=(cos(real(xt(n,2)))*sin(real(xt(n,1))/2.))*(1.+0.5*sin(time))
        interv=(1.+0.5*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx/ep
        temp(n,1)=(cos(tau(n))*et(1,n)-sin(tau(n))*et(2,n))/bx
        temp(n,2)=(cos(tau(n))*et(2,n)+sin(tau(n))*et(1,n))/bx
        fy(n,1)=temp(n,1)+interv*yt(n,2)
        fy(n,2)=temp(n,2)-interv*yt(n,1)
    end

end 

function apply_bc()

    while ( xxt(1) > xmax )
        xxt(1) = xxt(1) - dimx
    end

    while ( xxt(1) < xmin )
        xxt(1)= xxt(1) + dimx
    end

    while ( xxt(2) > ymax )
        xxt(2)  = xxt(2)  - dimy
    end

    while ( xxt(2)  < ymin )
        xxt(2) = xxt(2)  + dimy
    end

end 

end 

test_efd( 16, 1 ) # trigger build

@time test_efd( 16, 204800 )
