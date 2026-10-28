!> Routines to import a restart file written out by a routine in [[export_restart]].
module mod_import_grid
implicit none

character(len=20), parameter :: rst_file_ind_fmt(2) = (/'(a,i6.6)', '(a,i5.5)'/)

contains
!> Imports a restart file written out by the routine export_restart.

subroutine import_grid(node_list, element_list, filename, format_rst, ierr)

  use tr_module
  use data_structure
  use phys_module
  use pellet_module
  use equil_info
  use mod_boundary
  use basis_at_gaussian

  implicit none
  
  ! --- Routine parameters
  type(type_node_list),         intent(inout)           :: node_list
  type(type_element_list),      intent(inout)           :: element_list
  character*(*)          ,      intent(in)              :: filename
  integer,                      intent(out)             :: ierr
  integer,                      intent(in)              :: format_rst  ! format of restart file 
 
  ! --- Local parameters
  type (type_bnd_element_list)           :: bnd_elm_list    
  type (type_bnd_node_list)              :: bnd_node_list 

  ! Initialise basis functions before element tree is populated
  call initialise_basis()
  

  write(*,*) " Restart from HDF5 file " // trim(filename) // '.h5'
  call import_hdf5_grid(node_list, element_list, &
      trim(filename)//'.h5', format_rst, ierr)
      
  
  ! --- Required initializations to update equilibrium state
  call initialise_basis
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  call update_equil_state(0,node_list, element_list, bnd_elm_list, xpoint, xcase)
  write(*,*) " "
  write(*,*) " The equilibrium state has been updated "
  
end subroutine import_grid


!
! Import an HDF5 restart file
subroutine import_hdf5_grid(node_list, element_list, filename, format_rst, error)

#include "version.h"

  use tr_module 
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only: import_HDF5_restart_vacuum, current_FB_fact
  use mod_element_rtree, only: populate_element_rtree
  use mod_neighbours, only: update_neighbours
#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
  use mod_parameters 
#endif
  
  implicit none
  
  ! --- Routine parameters
  type(type_node_list),         intent(inout)           :: node_list
  type(type_element_list),      intent(inout)           :: element_list
  character(len=*),             intent(in)              :: filename
  integer,                      intent(in)              :: format_rst  ! format of restart file
  integer,                      intent(out)             :: error
  
  ! --- Perturbation-Import variables

  ! --- Local variables
  integer              :: i, j, m, k, n_tor_tmp, n_coord_tor_tmp, jorek_model_tmp, n_var_tmp, n_order_tmp, n_period_tmp, rst_hdf5_version_tmp, i_p, p_begin
  integer              :: n_plane_tmp, n_vertex_max_tmp, n_nodes_max_tmp, n_elements_max_tmp,n_boundary_max_tmp, n_nodes_tmp, n_dof_tmp
  integer              :: n_pieces_max_tmp, n_degrees_tmp, nref_max_tmp, n_ref_list_tmp, n_new_modes
  real*8               :: growth_mag, growth_kin, amplitude
  character*50         :: version_control, version_control_tmp, t_treat_axis
  logical              :: kept, modes_changed, import_3xx_4xx
  
#ifdef USE_HDF5
  integer(HID_T)     :: file_id, datatype, dataset
  integer            :: ind, n_spi_check, n_inj_check
  character          :: t_current_prof_initialized
  integer            :: var_rank       !< Rank of an array in the restart file
  integer(HSIZE_T),dimension(:),allocatable :: var_dims       !< Dimension of an array in the restart file
  
  real(RKIND), allocatable :: t_x(:,:,:,:)

  integer,     allocatable :: t_index(:,:)
  integer,     allocatable :: t_boundary(:)
  character,   allocatable :: t_axis_node(:)     
  integer,     allocatable :: t_axis_dof(:)

  integer,     allocatable :: t_vertex(:,:)
  real(RKIND), allocatable :: t_size(:,:,:)

#endif
  error = 0

  
#ifdef USE_HDF5

  ! ->  Reading HDF5 file
  write(*,*) 'Importing HDF5 restart file "', trim(filename), '".' 
  
  ! -> Open HDF5 file
  call HDF5_open(trim(filename),file_id,error)
  if ( error /= 0 ) then
    write(*,*) '...failed!'
    return
  end if

    
  !write(*,'(2(A,i5))') ' Importing ',n_tor_tmp,' harmonics with n_period=', n_period_tmp 
  call HDF5_integer_reading(file_id,n_nodes_tmp,"n_nodes")
  call HDF5_integer_reading(file_id,n_dof_tmp,"n_dof")
  n_degrees_tmp =n_degrees
  write(*,*) "n_dof", n_dof_tmp, n_degrees
  node_list%n_nodes =  n_nodes_tmp
  write(*,*) "NODES", n_nodes_tmp, node_list%n_nodes
  call HDF5_integer_reading(file_id,element_list%n_elements,"n_elements")

  ! initialise and allocate node_list
  call init_node_list(node_list, n_nodes_tmp, n_dof_tmp, n_var)

  ! -> Allocate temporary arrays 
  call tr_allocate(t_x,     1,node_list%n_nodes,1,n_coord_tor,1,n_degrees,1,n_dim,         "node_list%x",     CAT_UNKNOWN)
  call tr_allocate(t_index,      1,node_list%n_nodes,1,n_degrees_tmp,"index",      CAT_UNKNOWN)
  call tr_allocate(t_boundary,   1,node_list%n_nodes,                "boundary",   CAT_UNKNOWN)

  ! type_element, element_list%n_elements
  call tr_allocate(t_vertex,      1,element_list%n_elements,1,n_vertex_max,                 "vertex",CAT_UNKNOWN)
  call tr_allocate(t_size,        1,element_list%n_elements,1,n_vertex_max,1,n_degrees_tmp, "size",CAT_UNKNOWN)
  write(*,*) "x"
  call HDF5_array3D_reading(file_id,t_x(:,1,:,:),        'x')
  write(*,*) "index"
  call HDF5_array2D_reading_int (file_id,t_index,       'index')
  write(*,*) "boundary"
  call HDF5_array1D_reading_int (file_id,t_boundary,    'boundary')
  write(*,*) "axis_node"
  call tr_allocate(t_axis_node,  1,node_list%n_nodes,                "axis_node",  CAT_UNKNOWN)
  call HDF5_array1D_reading_char(file_id,t_axis_node,   'axis_node')
  write(*,*) "axis_dof"


  
  do i=1,node_list%n_nodes
    do j=1,n_degrees_tmp
      node_list%node(i)%x(:,j,:)  = t_x(i,:,j,:) 
    enddo
    write(12,*) t_x(i,1,1,1), t_x(i,1,1,2)
    write(13,*) t_index(i,:)
    write(14,*) t_boundary(i)
    node_list%node(i)%values = 0.d0 
    node_list%node(i)%deltas = 0.d0 

    node_list%node(i)%index(1:n_degrees_tmp) = t_index(i,1:n_degrees_tmp)
    node_list%node(i)%boundary = t_boundary(i)
    if (t_axis_node(i) == 'T') then
       node_list%node(i)%axis_node = .true.
    else
       node_list%node(i)%axis_node = .false.
    end if
    node_list%node(i)%axis_dof = 0
    node_list%node(i)%parents = 0
    node_list%node(i)%parent_elem = 0
    node_list%node(i)%ref_lambda = 0
    node_list%node(i)%ref_mu = 0
    node_list%node(i)%constrained = .false.
  end do

  call HDF5_array2D_reading_int(file_id,t_vertex,      'vertex')
  call HDF5_array3D_reading    (file_id,t_size,        'size')


  do i=1,element_list%n_elements
    write(15,*) t_vertex(i,:)
    write(16,*) t_size(i,:,:)
    element_list%element(i)%vertex                  = t_vertex(i,:)
    element_list%element(i)%size(:,1:n_degrees_tmp) = t_size(i,:,1:n_degrees_tmp)
    element_list%element(i)%father                  = 0
    element_list%element(i)%n_sons                  = 0
    element_list%element(i)%n_gen                   = 0
    element_list%element(i)%sons                    = 0
    element_list%element(i)%contain_node            = 0
    element_list%element(i)%nref                    = 0
  end do

!  call  update_neighbours(node_list,element_list, force_rtree_initialize=.true.)

    
  treat_axis = .false.
  
  call HDF5_close(file_id)
  
  call tr_deallocate(t_x,"t_x",CAT_UNKNOWN)
 
  call tr_deallocate(t_index,"index",CAT_UNKNOWN)
  call tr_deallocate(t_boundary,"boundary",CAT_UNKNOWN)
  call tr_deallocate(t_vertex,"t_vertex",CAT_UNKNOWN)
  call tr_deallocate(t_size,"t_size",CAT_UNKNOWN)

  write(*,*) ' restart complete '

#endif
  

  return
end subroutine import_hdf5_grid







!< Checks if a restart file exists in the current directory
!< Returns -1 if not found, and the digit format index if found (1 for 6 digits), (2 for 5 digits)
integer function restart_file_exists(i_step)

  use phys_module, only : rst_hdf5

  implicit none

  integer, intent(in) :: i_step
  integer             :: i_fmt
  character(len=64)   :: file_name, extension
  logical             :: file_exists

  restart_file_exists = -1

  ! Determine the file extension
  extension = '.rst'
  if (rst_hdf5 .ne. 0) extension = '.h5'

  ! Check each possible format
  do i_fmt = 1, size(rst_file_ind_fmt)
    write(file_name, rst_file_ind_fmt(i_fmt)) 'jorek', i_step
    inquire(file=trim(file_name) // extension, exist=file_exists)

    if (file_exists) then
      restart_file_exists = i_fmt
      return
    end if
  end do

end function restart_file_exists


end module mod_import_grid
