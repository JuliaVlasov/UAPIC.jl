module mesh_fields_m

implicit none

type :: mesh_t

    real(8) :: xmin 
    real(8) :: xmax 
    integer :: nx   
    real(8) :: dx   
    real(8) :: ymin 
    real(8) :: ymax 
    integer :: ny   
    real(8) :: dy   

end type mesh_t

type :: mesh_fields_t

    type(mesh_t)         :: mesh 
    real(8), allocatable :: e(:,:,:)
    real(8), allocatable :: rho(:,:)

end type mesh_fields_t

contains

    subroutine init_mesh( self, xmin, xmax, nx, ymin, ymax, ny )

        type(mesh_t) :: self

        real(8), intent(in) :: xmin 
        real(8), intent(in) :: xmax 
        integer, intent(in) :: nx   
        real(8), intent(in) :: ymin 
        real(8), intent(in) :: ymax 
        integer, intent(in) :: ny   

        self%xmin = xmin
        self%xmax = xmax
        self%nx   = nx
        self%ymin = ymin
        self%ymax = ymax
        self%ny   = ny
        self%dx   = (xmax - xmin) / real(nx,8)
        self%dy   = (ymax - ymin) / real(ny,8)

    end

    subroutine init_mesh_fields( self, mesh )

         type(mesh_fields_t) :: self
         type(mesh_t)        :: mesh
	    
         self%mesh = mesh

         allocate(self%e(2, mesh%nx+1, mesh%ny+1))
         allocate(self%rho(mesh%nx+1, mesh%ny+1))

    end

end module mesh_fields_m
