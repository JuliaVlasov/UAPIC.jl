module poisson_m

use, intrinsic :: iso_c_binding 

implicit none

include 'fftw3.f03'

type :: poisson_t

    type(mesh_t)            :: mesh
    real(8),    allocatable :: kx(:,:)
    real(8),    allocatable :: ky(:,:)
    complex(8), allocatable :: rho(:,:)
    complex(8), allocatable :: ex(:,:)
    complex(8), allocatable :: ey(:,:)
    integer(8)              :: fw
    integer(8)              :: bw


end type poisson_t

contains

    subroutine init_poisson( self, mesh )
        
        type(poisson_t)      :: self
        type(mesh_t)         :: mesh

        real(8)              :: kx0
        real(8)              :: ky0
        real(8)              :: pi
        real(8), allocatable :: tmp(:,:)

        pi = 4d0 * atan(1d0)

        kx0 = 2d0 * pi / (mesh%xmax - mesh%xmin)
        ky0 = 2d0 * pi / (mesh%ymax - mesh%ymin)

        allocate(self%kx(nx÷2+1,ny))
        allocate(self%ky(nx÷2+1,ny))

        do ik=1,mesh%nx/2+1
           kx1 = (ik-1)*kx0
           do jk = 1,mesh%ny/2
              self%kx(ik,jk) = kx1
              self%ky(ik,jk) = (jk-1)*ky0
           end do
           do jk = mesh%ny/2+1,mesh%ny
              self%kx(ik,jk) = kx1
              self%ky(ik,jk) = (jk-1-ny)*ky0
           end do
        end do

        self%kx(1,1) = 1.0
        k2 = self%kx * self%kx + self%ky * self%ky
        self%kx = self%kx / k2
        self%ky = self%ky / k2

        allocate(self%rho(nx/2+1,ny))
        allocate(self%ex(nx/2+1,ny))
        allocate(self%ey(nx/2+1,ny))
        allocate(tmp(mesh%nx, mesh%ny))

        call dfftw_plan_dft_r2c_2d(self%fw,nx,ny,tmp,self%rho,FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(self%bw,nx,ny,self%rho,tmp,FFTW_ESTIMATE)

    end subroutine init_poisson	    

    subroutine solve_poisson ( self, fields )

        type(poisson_t)       :: self
        type(mesh_fields_t)   :: fields
        complex(8), parameter :: j = (0d0, 1d0)

        nx = fields.mesh.nx
        ny = fields.mesh.ny

        call dfftw_execute_dft( self.fw, fields.rho(1:nx,1:ny), self.rho )

        self%ex = -j * self%kx * self%rho
        self%ey = -j * self%ky * self%rho

        call dfftw_execute_dft( self.bw, self%ex, fields%e(1,1:nx,1:ny) )
        call dfftw_execute_dft( self.fw, self.ey, fields%e(2,1:nx,1:ny) )

        fields%e(1,nx+1,:) = fields%e(1,1,:)
        fields%e(1,:,ny+1) = fields%e(1,:,1)
        fields%e(2,nx+1,:) = fields%e(2,1,:)
        fields%e(2,:,ny+1) = fields%e(2,:,1)

    end subroutine solve_poisson

end module poisson_m
