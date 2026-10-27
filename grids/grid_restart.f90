subroutine grid_restart(node_list, element_list)


use constants
use tr_module 
use data_structure
use mod_neighbours, only: update_neighbours
use mod_interp
use phys_module, only: force_central_node, write_ps, fix_axis_nodes, treat_axis,rst_format, xcase
use mod_grid_conversions
use mod_poiss
use mod_node_indices
use equil_info, only: find_xpoint
use mod_import_grid

implicit none

! --- Routine parameters
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list

! --- local variables
type (type_surface_list) :: flux_list

type (type_node_list),    pointer :: newnode_list
type (type_element_list), pointer :: newelement_list

! --- Unused (just for call to Poisson for psi-projection)
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

real*8, allocatable :: s_values(:), theta_sep(:), R_sep(:), Z_sep(:), R_max(:), Z_max(:), R_min(:), Z_min(:),s_tmp(:)
real*8              :: psi_axis, R_axis, Z_axis, s_axis, t_axis, R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), psi_xpoint(2)
real*8              :: s_find(8), t_find(8), st_find(8), tht_x, theta, delta, ss, tmp1, tmp2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8,allocatable  :: R_polar(:,:,:),Z_polar(:,:,:),xout(:),xp(:),yp(:)
real*8              :: R_cub1d(4), Z_cub1d(4), dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z
real*8, allocatable :: RR_new(:,:),ZZ_new(:,:),s_flux(:,:),t_flux(:,:),t_tht(:,:)
integer,allocatable :: ielm_flux(:,:), keep(:,:,:), k_cross(:,:)
integer             :: i, j, k, l, m, n_psi, n_flux_2, n_open_2, n_tht_2, n_psi_2, i2, j2
integer             :: n_private_2, i_surf, n_pieces, ierr
integer             :: i_elm_axis, i_elm_xpoint(2), i_elm_find(8), i_sep, i_max, i_find, npl, ifail
integer             :: node, index, node_start, index_xpoint, n_xpoint, j_start, j_end
integer             :: iv, ivp, node_iv, node_ivp, i_elm, n_leg_2, n_tht_3
integer             :: my_id, ielm_out
real*8              :: Rmid, Zmid, R0,Z0, RP,ZP, dR0, dZ0, dRP, dZP, size_0, size_p, denom
real*8              :: R1, Z1, s_out, t_out, R_out, Z_out
real*8              :: EJAC, RX, RY, SX, SY, CRR, CZZ, CRZ, alpha1, alpha2, alpha_max, alpha_min
real*8              :: RL1, RL2, RL3, RL4, RL5, RL6, RL7, RL8, RL9, ZL1, ZL2, ZL3, ZL4, ZL5, ZL6, ZL7, ZL8, ZL9
real*8              :: angle_L1, angle_L8, angle_L9, rr1, ss1
logical             :: xpoint
real*8,external     :: root
character*4         :: label
integer             :: i_elm1, i_vertex1, i_node1, i_node_save
integer             :: i_elm2, i_vertex2, i_node2
integer             :: node_indices( (n_order+1)/2, (n_order+1)/2 ), ii, jj

xpoint = .true.
my_id  = 0

write(*,*) '*************************************'
write(*,*) '*          X-point grid             *'
write(*,*) '*************************************'


!***********************************************************************
!*     define the new nodes and finite elements    (nodes first)       *
!***********************************************************************

! --- Allocate data structures for new nodes and elements and initialize them.
allocate(newnode_list)
call init_node_list(newnode_list, n_nodes_max, newnode_list%n_dof, n_var)
call tr_register_mem(sizeof(newnode_list),"newnode_list",CAT_GRID)
allocate(newelement_list)
call tr_register_mem(sizeof(newelement_list),"newelement_list",CAT_GRID)

newnode_list%n_nodes = 0
newnode_list%n_dof   = 0
do i = 1, n_nodes_max
  newnode_list%node(i)%x           = 0.d0
  newnode_list%node(i)%values      = 0.d0
  newnode_list%node(i)%deltas      = 0.d0
  newnode_list%node(i)%index       = 0
  newnode_list%node(i)%boundary    = 0
  newnode_list%node(i)%parents     = 0
  newnode_list%node(i)%parent_elem = 0
  newnode_list%node(i)%ref_lambda  = 0.d0
  newnode_list%node(i)%ref_mu      = 0.d0
  newnode_list%node(i)%constrained = .false.
end do
newelement_list%n_elements = 0
do i = 1, n_elements_max
  newelement_list%element(i)%vertex       = 0
  newelement_list%element(i)%neighbours   = 0
  newelement_list%element(i)%size         = 0.d0
  newelement_list%element(i)%father       = 0
  newelement_list%element(i)%n_sons       = 0
  newelement_list%element(i)%n_gen        = 0
  newelement_list%element(i)%sons         = 0
  newelement_list%element(i)%contain_node = 0
  newelement_list%element(i)%nref         = 0
end do

call import_grid(newnode_list, newelement_list, "jorek000000_edited_3", rst_format, ierr)

do i=1,newnode_list%n_nodes

  R1 = newnode_list%node(i)%x(1,1,1)
  Z1 = newnode_list%node(i)%x(1,1,2)

  call find_RZ(node_list,element_list,R1,Z1,R_out,Z_out,ielm_out,s_out,t_out,ifail)

  call interp_RZ(node_list,element_list,ielm_out,s_out,t_out,   &
                 RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

  call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

  RZ_jac  = dRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
  PSI_R  = (   dZZg1_ds * dPSg1_dr - dZZg1_dr * dPSg1_ds ) / RZ_jac
  PSI_Z  = ( - dRRg1_ds * dPSg1_dr + dRRg1_dr * dPSg1_ds ) / RZ_jac

  newnode_list%node(i)%values(1,1,1) = PSg1
  newnode_list%node(i)%values(1,2,1) = PSI_R * newnode_list%node(i)%x(1,2,1) + PSI_Z * newnode_list%node(i)%x(1,2,2)
  newnode_list%node(i)%values(1,3,1) = PSI_R * newnode_list%node(i)%x(1,3,1) + PSI_Z * newnode_list%node(i)%x(1,3,2)
  newnode_list%node(i)%values(1,4,1) = PSI_R * newnode_list%node(i)%x(1,4,1) + PSI_Z * newnode_list%node(i)%x(1,4,2)

  if (newnode_list%node(i)%boundary .eq. 2) newnode_list%node(i)%values(1,3,1) = 0.d0

enddo



!----------------------------- empty old nodes/elements
do i=1,node_list%n_nodes
  node_list%node(i)%x        = 0.d0
  node_list%node(i)%values   = 0.d0
  node_list%node(i)%index    = 0
  node_list%node(i)%boundary = 0
enddo
node_list%n_nodes = 0

do i=1,element_list%n_elements
  element_list%element(i)%vertex     = 0
  element_list%element(i)%size       = 0.d0
  element_list%element(i)%neighbours = 0
enddo

!---------------------------- copy new grid into nodes/elements
node_list%n_nodes = newnode_list%n_nodes
node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)

element_list%n_elements = newelement_list%n_elements
element_list%element(1:element_list%n_elements) = newelement_list%element(1:element_list%n_elements)


!----temporary, needs to be completed, neighbour information
do i=1, element_list%n_elements
  element_list%element(i)%father = 0
  element_list%element(i)%n_sons = 0
  element_list%element(i)%sons(:) = 0
enddo

call tr_unregister_mem(sizeof(newnode_list),"newnode_list",CAT_GRID)
call dealloc_node_list(newnode_list) ! deallocates all the node values in newnode_list
deallocate(newnode_list)
call tr_unregister_mem(sizeof(newelement_list),"newelement_list",CAT_GRID)
deallocate(newelement_list)

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_restart
