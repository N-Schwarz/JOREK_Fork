!> mod_bremsstrahlung_light_dist_vertices_test contains all variables
!> and procedures for testing the bremsstrahlung light vertices
module mod_bremsstrahlung_light_dist_vertices_test
use fruit
use mod_particle_types,                    only: particle_kinetic_relativistic_id
use mod_particle_sim,                      only: particle_sim
use mod_spectra_monte_carlo,               only: spectrum_rng_uniform
use mod_bremsstrahlung_light_dist_vertices,only: bremsstrahlung_light_dist
implicit none

private
public :: run_fruit_bremsstrahlung_light_dist_vertices

!> Variables ---------------------------------------------------------
!> general parameters
logical,parameter                           :: use_xor_time_pid=.true.
real*8,parameter                            :: tol_real8=2.5d-11
real*8,parameter                            :: tol2_real8=3.5d-8
real*8,parameter                            :: mass_RE=5.48579909065d-4 !< electron mass [AMU]
real*8,parameter                            :: Zeff_sol=1.7d0
!> parameters for generating bremsstrahlung lights (single particle type: kinetic relativistic)
integer,parameter :: n_mhd_sol=1 !< local electron density
integer,parameter :: n_x=3
integer,parameter :: n_properties=8
integer,parameter :: fill_type_base=1 !< use cylindrical initialisation
integer,parameter :: n_times_sol=2
integer,parameter :: n_particle_types_check_sol=1
integer,dimension(n_times_sol),parameter      :: n_groups_per_sim=(/2,1/)
integer,parameter                             :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/97,158,211,0/),shape(n_particles_per_group))
integer,parameter                             :: n_particles_max=maxval(n_particles_per_group)
real*8,parameter                              :: survival_threshold=0.33
!> parameters for generating spectra (photon energy [J], ~20keV-500keV HXR range)
integer,parameter                             :: n_spectra=2
integer,parameter                             :: n_lines_per_spectrum=11
real*8,parameter                              :: eV2J=1.602176565d-19
real*8,dimension(n_spectra),parameter         :: min_k=(/2.0d4*eV2J,5.0d4*eV2J/)
real*8,dimension(n_spectra),parameter         :: max_k=(/1.5d5*eV2J,5.0d5*eV2J/)
!> parameters for generating shadowed points
integer,parameter                             :: n_shadowed_per_particle=5
real*8,parameter                              :: maximum_cos_half_angle_sol=1d0
real*8,dimension(2),parameter                 :: length_shadowed=(/2.d-1,7.d0/)
!> variables for generating bremsstrahlung lights
type(bremsstrahlung_light_dist)               :: vertex_sol
type(particle_sim),dimension(n_times_sol)     :: sims_particles
integer,dimension(n_particle_types_check_sol) :: particle_types_check_sol
integer,dimension(n_times_sol)                :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol)   :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)                 :: time_vector_sol
real*8,dimension(n_particles_max*n_groups_max,n_times_sol)     :: weight_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol) :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_sol
!> variables for generating spectra
type(spectrum_rng_uniform)                    :: spectrum
!> variables for generating shadowed points
real*8,dimension(:,:,:,:),allocatable         :: x_shadowed

!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
subroutine run_fruit_bremsstrahlung_light_dist_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: bremsstrahlung light vertices tests"
  call setup
  write(*,'(/A)') "  ... running: bremsstrahlung light vertices tests"
  call run_test_case(test_setup_bremsstrahlung_radiation_class,&
  'test_setup_bremsstrahlung_radiation_class')
  call run_test_case(test_compute_bremsstrahlung_mhd_fields,&
  'test_compute_bremsstrahlung_mhd_fields')
  call run_test_case(test_compute_bremsstrahlung_light_properties,&
  'test_compute_bremsstrahlung_light_properties')
  call run_test_case(test_bremsstrahlung_irradiance_directional_func,&
  'test_bremsstrahlung_irradiance_directional_func')
  call run_test_case(test_compute_particle_from_bremsstrahlung_light,&
  'test_compute_particle_from_bremsstrahlung_light')
  write(*,'(/A)') "  ... tearing-down: bremsstrahlung light vertices tests"
  call teardown
end subroutine run_fruit_bremsstrahlung_light_dist_vertices

!> Set-up and tear-down procedures------------------------------------
subroutine setup()
  use mod_rng,                        only: type_rng
  use mod_pcg32_rng,                  only: pcg32_rng
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_common_test_tools,          only: omp_initialize_rngs
  use mod_particle_common_test_tools, only: sim_time_interval
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: fill_groups,fill_mass_RE
  use mod_particle_common_test_tools, only: fill_particles_tokamak
  use mod_particle_common_test_tools, only: invalidate_particles
  use mod_particle_common_test_tools, only: obtain_active_particle_ids
  use mod_fields_linear,              only: jorek_fields_interp_linear
  !$ use omp_lib
  implicit none
  !> variables
  integer :: ii,ifail,n_particles_RE_max,n_threads
  integer,dimension(n_groups_max,n_times_sol) :: particle_types_sol
  class(type_rng),dimension(:),allocatable :: rngs
  !> initialisation
  particle_types_check_sol = (/particle_kinetic_relativistic_id/)
  particle_types_sol = particle_kinetic_relativistic_id
  vertex_sol%n_property_vertex = n_properties; vertex_sol%Zeff = Zeff_sol; ifail = 0;
  n_active_particles_sol = 0; n_threads = 1;
  !$ n_threads = omp_get_max_threads()
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)

  !> allocate and initialise particle list
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(jorek_fields_interp_linear::sims_particles(ii)%fields)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_groups(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_mass_RE(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_particles_tokamak(n_groups_per_sim(ii),sims_particles(ii)%groups,fill_type_base)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
  enddo
  n_particles_RE_max = maxval(n_active_vertices_sol)
  !> initialise positions and properties tables
  call compute_brems_x_properties_ana()

  !> initialise monte-carlo spectra
  allocate(pcg32_rng::rngs(n_threads))
  call omp_initialize_rngs(n_lines_per_spectrum,n_threads,rngs,use_xor_time_pid_in=use_xor_time_pid)
  spectrum = spectrum_rng_uniform(n_lines_per_spectrum,n_spectra,min_k,max_k)
  call spectrum%generate_spectrum(rngs)
  deallocate(rngs)

  !> generate shadowed points positions
  allocate(x_shadowed(n_x,n_shadowed_per_particle,n_particles_RE_max,n_times_sol))
  call compute_x_shadowed_particles()
end subroutine setup

subroutine teardown()
  implicit none
  call vertex_sol%deallocate_vertices; call spectrum%deallocate_spectrum;
  deallocate(x_shadowed);
end subroutine teardown

!> Tests -------------------------------------------------------------
!> Test setup bremsstrahlung radiation class
subroutine test_setup_bremsstrahlung_radiation_class()
  implicit none
  call vertex_sol%setup_light_class
  call assert_equals(n_properties,vertex_sol%n_property_vertex,&
  "Error check setup bremsstrahlung light class: wrong size of the vertex properties array!")
  call assert_equals(n_mhd_sol,vertex_sol%n_mhd,&
  "Error check setup bremsstrahlung light class: wrong size of the mhd array!")
  call assert_equals(n_particle_types_check_sol,vertex_sol%n_particle_types,&
  "Error check setup bremsstrahlung light class: wrong size of the particle types array!")
  call assert_equals(particle_types_check_sol,vertex_sol%particle_types,&
  n_particle_types_check_sol,"Error check setup bremsstrahlung light class: wrong particle types list!")
end subroutine test_setup_bremsstrahlung_radiation_class

!> test the calculation of the local electron density
subroutine test_compute_bremsstrahlung_mhd_fields()
  use mod_particle_common_test_tools, only: compute_test_n_e
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  real*8,dimension(n_mhd_sol,n_particles_max*n_groups_max) :: mhd_fields_sol
  real*8,dimension(n_mhd_sol,n_particles_max*n_groups_max) :: mhd_fields_test
  do kk=1,n_times_sol
    counter = 0; mhd_fields_sol = 0d0; mhd_fields_test = 0d0;
    do jj=1,n_groups_per_sim(kk)
      do ii=1,n_particles_per_group(jj,kk)
        if(sims_particles(kk)%groups(jj)%particles(ii)%i_elm.le.0) cycle
        counter = counter + 1
        call compute_test_n_e(sims_particles(kk)%groups(jj)%particles(ii)%x,mhd_fields_sol(1,counter))
        call vertex_sol%compute_mhd_fields(sims_particles(kk)%fields,&
        sims_particles(kk)%groups(jj)%particles(ii),kk,&
        sims_particles(kk)%groups(jj)%mass,mhd_fields_test(:,counter))
      enddo
    enddo
    call assert_equals(mhd_fields_sol,mhd_fields_test,n_mhd_sol,&
    n_particles_max*n_groups_max,tol_real8,&
    "Error bremsstrahlung light compute mhd fields: too large errors!")
  enddo
end subroutine test_compute_bremsstrahlung_mhd_fields

!> test the property function of bremsstrahlung light properties
subroutine test_compute_bremsstrahlung_light_properties()
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  real*8,dimension(n_mhd_sol) :: mhd_fields
  real*8,dimension(n_properties,n_particles_max*n_groups_max) :: error,zeros
  zeros = 0.d0
  do kk=1,n_times_sol
    counter = 0; error = 0.d0;
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is(particle_kinetic_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.le.0) cycle
          counter = counter + 1
          call compute_brems_mhd_ana(p_list(ii),mhd_fields)
          call vertex_sol%compute_light_properties(counter,kk,p_list(ii),&
          sims_particles(kk)%groups(jj)%mass,mhd_fields)
        enddo
      end select
    enddo
    where(properties_sol(:,:,kk).ne.0d0) error = abs((properties_sol(:,:,kk)-&
    vertex_sol%properties(:,:,kk))/properties_sol(:,:,kk))
    call assert_equals(zeros,error,n_properties,n_particles_max*n_groups_max,&
    tol_real8,"Error bremsstrahlung light compute properties: too large errors!")
  enddo
end subroutine test_compute_bremsstrahlung_light_properties

!> test the bremsstrahlung light irradiance and directional functions
subroutine test_bremsstrahlung_irradiance_directional_func()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer :: ii,jj,kk
  integer,dimension(n_times_sol) :: n_particles_time
  real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol)          :: x_cart_loc
  real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_loc
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shadowed_per_particle) :: dir_fun,irradiance
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shadowed_per_particle) :: dir_fun_sol,irradiance_sol
  x_cart_loc = 0.d0; properties_loc = 0.d0; n_particles_time = 0;
  do ii=1,n_times_sol
    n_particles_time(ii) = n_active_vertices_sol(ii)
    x_cart_loc(:,1:n_particles_time(ii),ii) = x_cart_sol(:,1:n_particles_time(ii),ii)
    properties_loc(:,1:n_particles_time(ii),ii) = properties_sol(:,1:n_particles_time(ii),ii)
  enddo
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles)

  do kk=1,n_times_sol
    do jj=1,n_particles_time(kk)
      dir_fun_sol = 0d0; dir_fun = 0d0; irradiance = 0d0;
      do ii=1,n_shadowed_per_particle
        call compute_brems_directionality_irradiance(x_shadowed(:,ii,jj,kk),&
        x_cart_loc(:,jj,kk),properties_loc(:,jj,kk),dir_fun_sol(:,:,ii),irradiance_sol(:,:,ii))
        call vertex_sol%directionality_funct(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),dir_fun(:,:,ii))
        call vertex_sol%spectral_irradiance(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),irradiance(:,:,ii))
      enddo
      call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
      n_shadowed_per_particle,dir_fun_sol,dir_fun,tol2_real8,&
      "Error bremsstrahlung directionality function: directionality function mismatch!")
      call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
      n_shadowed_per_particle,irradiance_sol,irradiance,tol2_real8,&
      "Error bremsstrahlung irradiance: irradiance mismatch!")
    enddo
  enddo
end subroutine test_bremsstrahlung_irradiance_directional_func

!> test the reconstruction of the emitting electron from bremsstrahlung lights
subroutine test_compute_particle_from_bremsstrahlung_light()
  use mod_particle_types,        only: particle_base,particle_kinetic_relativistic
  implicit none
  !> variables
  class(particle_base),dimension(:),allocatable :: particle_test
  integer                                       :: ii,jj,kk,counter
  real*8,dimension(n_mhd_sol)                   :: mhd_fields
  real*8                                         :: err_p,err_w
  do kk=1,n_times_sol
    counter = 0;
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is(particle_kinetic_relativistic)
        allocate(particle_kinetic_relativistic::particle_test(n_particles_per_group(jj,kk)))
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.le.0) cycle
          counter = counter + 1
          call compute_brems_mhd_ana(p_list(ii),mhd_fields)
          call vertex_sol%compute_light_properties(counter,kk,p_list(ii),&
          sims_particles(kk)%groups(jj)%mass,mhd_fields)
          call vertex_sol%compute_particle_from_light(sims_particles(kk)%fields,&
          counter,kk,sims_particles(kk)%groups(jj)%mass,particle_test(ii))
          select type (p_test=>particle_test(ii))
          type is (particle_kinetic_relativistic)
            err_p = norm2(p_test%p-p_list(ii)%p)/norm2(p_list(ii)%p)
            err_w = abs((p_test%weight-p_list(ii)%weight)/p_list(ii)%weight)
            call assert_true(err_p.lt.tol2_real8,&
            "Error reconstruct particle from bremsstrahlung light: momentum mismatch!")
            call assert_true(err_w.lt.tol2_real8,&
            "Error reconstruct particle from bremsstrahlung light: weight mismatch!")
            call assert_equals(int(-1,kind=1),p_test%q,&
            "Error reconstruct particle from bremsstrahlung light: wrong charge!")
          end select
        enddo
        deallocate(particle_test)
      end select
    enddo
  enddo
end subroutine test_compute_particle_from_bremsstrahlung_light

!> Tools -------------------------------------------------------------
!> compute the positions of points shadowed by particle lights,
!> taken within the ~1/gamma emission cone of each relativistic electron
subroutine compute_x_shadowed_particles()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  !> variables
  integer               :: ii,jj,kk,n_particles_time
  real*8                 :: cos_half_angle
  real*8,dimension(n_x)  :: p_dir,x_part,rng
  do kk=1,n_times_sol
    n_particles_time = n_active_vertices_sol(kk)
    do jj=1,n_particles_time
      x_part = x_cart_sol(:,jj,kk)
      p_dir  = properties_sol(1:3,jj,kk)
      cos_half_angle = cos(1.d0/(properties_sol(4,jj,kk)))
      do ii=1,n_shadowed_per_particle
        call random_number(rng)
        x_shadowed(:,ii,jj,kk) = sample_uniform_cone([cos_half_angle,&
        maximum_cos_half_angle_sol],rng,p_dir,x_part,length_shadowed)
      enddo
    enddo
  enddo
end subroutine compute_x_shadowed_particles

!> compute and fill particle positions and properties for RE
subroutine compute_brems_x_properties_ana()
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_particle_types,        only: particle_kinetic_relativistic
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  real*8,dimension(n_mhd_sol) :: mhd_fields
  x_cart_sol = 0.d0; properties_sol = 0.d0;
  do kk=1,n_times_sol
    counter = 0
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
      type is(particle_kinetic_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.le.0) cycle
          counter = counter + 1
          x_cart_sol(:,counter,kk) = cylindrical_to_cartesian(p_list(ii)%x)
          call compute_brems_mhd_ana(p_list(ii),mhd_fields)
          call compute_brems_properties_ana_1p(&
          sims_particles(kk)%groups(jj)%mass,p_list(ii),mhd_fields,&
          weight_sol(counter,kk),properties_sol(:,counter,kk))
        enddo
      end select
    enddo
  enddo
end subroutine compute_brems_x_properties_ana

!> compute the analytical (test) local electron density at a particle location
subroutine compute_brems_mhd_ana(particle,mhd_fields)
  use mod_particle_types,             only: particle_base
  use mod_particle_common_test_tools, only: compute_test_n_e
  implicit none
  class(particle_base),intent(in) :: particle
  real*8,dimension(n_mhd_sol),intent(out) :: mhd_fields
  call compute_test_n_e(particle%x,mhd_fields(1))
end subroutine compute_brems_mhd_ana

!> compute the bremsstrahlung lights directionality function and irradiance
!> using an independent re-implementation of Hoppe et al 2018 Eq.9
subroutine compute_brems_directionality_irradiance(x_illum,&
x_light,property,dir_func,irradiance)
  implicit none
  !> inputs
  real*8,dimension(n_x),intent(in)          :: x_illum,x_light
  real*8,dimension(n_properties),intent(in) :: property
  !> outputs
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(out) :: dir_func,irradiance
  !> variables
  integer               :: ii,jj
  integer,dimension(0)  :: int_param
  real*8                :: costheta
  real*8                :: gamma,p_norm,rest_mass_energy,omega_cone,pf,beta_c_over_omega
  real*8                :: k_dim,gamma_p,p_p,eps,eps_p,L
  real*8                :: TT1,TT2,TT3,TT4,TT5,LT1,LT2,LTF,LT3,LT4,LT5,bracket
  !> initialisations
  dir_func = 0d0; irradiance = 0d0;
  !> check is the shaded point is in the ~1/gamma emission cone
  costheta = dot_product(x_illum-x_light,property(1:3))/norm2(x_illum-x_light)
  if((costheta.lt.0).or.((sqrt(1.d0-costheta**2)*property(4)).gt.1.d0)) return
  !> local copies
  gamma = property(4); p_norm = property(5); rest_mass_energy = property(6);
  omega_cone = property(7); pf = property(8);
  beta_c_over_omega = ((p_norm/gamma)*299792458.d0)/(rest_mass_energy*omega_cone)
  !> compute directionality function and irradiance
  do ii=1,spectrum%n_spectra
    do jj=1,spectrum%n_points
      k_dim = spectrum%points(jj,ii)/rest_mass_energy
      if((k_dim.le.0.d0).or.(k_dim.ge.(gamma-1.d0))) cycle
      gamma_p = gamma-k_dim; p_p = sqrt(gamma_p**2-1.d0);
      eps   = 2.d0*log(gamma  +p_norm)
      eps_p = 2.d0*log(gamma_p+p_p)
      L     = 2.d0*log((gamma*gamma_p-1.d0+p_norm*p_p)/k_dim)
      TT1 = 4.d0/3.d0
      TT2 =-2.d0*gamma*gamma_p*(p_norm**2+p_p**2)/((p_norm**2)*(p_p**2))
      TT3 = eps*gamma_p/(p_norm**3)
      TT4 = eps_p*gamma/(p_p**3)
      TT5 =-eps*eps_p/(p_norm*p_p)
      LT1 = 8.d0*gamma*gamma_p/(3.d0*p_norm*p_p)
      LT2 = (k_dim**2)*((gamma**2)*(gamma_p**2)+(p_norm**2)*(p_p**2))/((p_norm**3)*(p_p**3))
      LTF = k_dim/(2.d0*p_norm*p_p)
      LT3 = eps*(gamma*gamma_p+p_norm**2)/(p_norm**3)
      LT4 =-eps_p*(gamma*gamma_p+p_p**2)/(p_p**3)
      LT5 = 2.d0*k_dim*gamma*gamma_p/((p_norm**2)*(p_p**2))
      bracket = TT1+TT2+TT3+TT4+TT5 + L*(LT1+LT2+LTF*(LT3+LT4+LT5))
      dir_func(jj,ii)   = beta_c_over_omega*(p_p/(k_dim*p_norm))*bracket
      irradiance(jj,ii) = dir_func(jj,ii)*pf
    enddo
  enddo
end subroutine compute_brems_directionality_irradiance

!> compute bremsstrahlung electron properties for one particle
subroutine compute_brems_properties_ana_1p(mass,particle,mhd_fields,weight,property)
  use constants,           only: ATOMIC_MASS_UNIT,SPEED_OF_LIGHT,TWOPI
  use constants,           only: EL_RAD,ALPHA_FINE_STRUCTURE
  use mod_particle_types,  only: particle_kinetic_relativistic
  implicit none
  !> inputs-outputs
  type(particle_kinetic_relativistic),intent(inout) :: particle
  !> inputs
  real*8,intent(in)                         :: mass
  real*8,dimension(n_mhd_sol),intent(in)    :: mhd_fields
  !> outputs
  real*8,intent(out) :: weight
  real*8,dimension(n_properties),intent(out) :: property
  !> variables
  real*8 :: velocity,beta_gamma,gamma,beta,p_norm
  velocity   = norm2(particle%p)
  beta_gamma = velocity/SPEED_OF_LIGHT
  gamma      = sqrt(1.d0+(beta_gamma/mass)**2)
  beta       = beta_gamma/(mass*gamma)
  p_norm     = gamma*beta
  weight = particle%weight
  property(1:3) = particle%p/velocity
  property(4)   = gamma
  property(5)   = p_norm
  property(6)   = mass*ATOMIC_MASS_UNIT*(SPEED_OF_LIGHT**2)
  property(7)   = TWOPI/(gamma*(gamma+p_norm))
  property(8)   = weight*mhd_fields(1)*Zeff_sol*(EL_RAD**2)*ALPHA_FINE_STRUCTURE
end subroutine compute_brems_properties_ana_1p

!>--------------------------------------------------------------------
end module mod_bremsstrahlung_light_dist_vertices_test
