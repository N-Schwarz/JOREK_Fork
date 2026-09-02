!> the mod_bremsstrahlung_light_dist_vertices implements
!> variables and procedures defining the in-flight (thin-target)
!> bremsstrahlung hard x-ray light distribution of relativistic
!> kinetic particle light sources. The photon is assumed to be
!> emitted, like synchrotron radiation, within a forward cone of
!> half-angle ~1/gamma around the electron velocity direction
!> (hence this type extends synchrotron_light and reuses its
!> emission-cone check procedures unchanged). The photon number
!> spectrum emitted within that cone is obtained by taking the
!> angle-integrated Born-approximation (Bethe-Heitler) bremsstrahlung
!> photon spectrum and spreading it uniformly over the cone's solid
!> angle. The model used is:
!> M. Hoppe, O. Embreus, C. Paz-Soldan, R.A. Moyer and T. Fulop,
!> Nucl. Fusion, vol.58, 082001, 2018 (their Eq.9), cross-checked
!> term-by-term against the open-source implementation in the
!> ConeBremsstrahlungEmission class of the SOFT2 code
!> (https://github.com/hoppe93/SOFT2)
module mod_bremsstrahlung_light_dist_vertices
use mod_synchrotron_light_vertices, only: synchrotron_light
implicit none

private
public :: bremsstrahlung_light_dist

!> Variables ---------------------------------------
type,extends(synchrotron_light) :: bremsstrahlung_light_dist
  real*8 :: Zeff=1.d0 !< plasma effective charge for bremsstrahlung: n_e*Zeff = sum_i n_i*Z_i**2
  contains
  procedure,pass(light_vert) :: directionality_funct => &
                                bremsstrahlung_directionality_funct
  procedure,pass(light_vert) :: spectral_irradiance => &
                                bremsstrahlung_spectral_irradiance
  procedure,pass(light_vert) :: compute_mhd_fields => &
                                compute_bremsstrahlung_mhd_fields
  procedure,pass(light_vert) :: compute_light_properties => &
                                compute_bremsstrahlung_light_properties
  procedure,pass(light_vert) :: setup_light_class => &
                                setup_bremsstrahlung_light_class
  procedure,pass(light_vert) :: compute_particle_from_light => &
                                compute_particle_from_bremsstrahlung_light
end type bremsstrahlung_light_dist
!> Interfaces --------------------------------------

contains

!> Procedures --------------------------------------

!> bremsstrahlung_directionality_funct computes the bremsstrahlung
!> photon-number spectral distribution per unit of the light's
!> "pf" prefactor (light_vert%properties(8)), i.e. the quantity
!> that when multiplied by properties(8) gives the number of
!> photons emitted per unit time, per unit solid angle and per
!> unit photon energy towards the shaded point x_shaded. Outside
!> of the ~1/gamma emission cone (reusing the synchrotron cone
!> check) or for kinematically forbidden photon energies
!> (k>=gamma-1) the returned value is zero.
!> inputs:
!>   light_vert: (bremsstrahlung_light_dist) bremsstrahlung light sources
!>   spectra:    (spectrum_base) spectral intervals (photon energy [J]) and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs:
!>   light_vert: (bremsstrahlung_light_dist) bremsstrahlung light sources
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   light_dstb: (real8)(n_points,n_spectra) bremsstrahlung photon spectral-angular
!>               distribution per unit of the light's pf prefactor towards x_shaded
subroutine bremsstrahlung_directionality_funct(light_vert,spectra,time_id,&
light_id,x_shaded,light_dstb)
  use constants,    only: SPEED_OF_LIGHT
  use mod_spectra,  only: spectrum_base
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  class(bremsstrahlung_light_dist),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)             :: spectra
  !> inputs:
  integer,intent(in)                          :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  !> variables
  logical :: in_parallel
  integer :: ii,jj
  integer,dimension(0) :: int_param
  real*8  :: gamma,p_norm,rest_mass_energy,beta_c_over_omega
  real*8,dimension(light_vert%n_property_vertex) :: light_properties

  !> initialisations
  in_parallel = .false.; light_dstb = 0d0;
  !$ in_parallel = omp_in_parallel()
  light_properties = light_vert%properties(:,light_id,time_id)
  !> check if the shaded point is in the ~1/gamma emission cone (reuses
  !> the synchrotron cone check: real_param(1:3)=direction, (4)=gamma)
  if(.not.light_vert%check_x_shaded_in_emission_zone(light_vert%n_x,x_shaded,&
  light_vert%x(:,light_id,time_id),0,4,int_param,light_properties(1:4))) return
  !> local copies of the electron state
  gamma = light_properties(4); p_norm = light_properties(5);
  rest_mass_energy = light_properties(6)
  !> beta*c / (rest_mass_energy * Omega_cone) prefactor, common to all spectral points
  beta_c_over_omega = ((p_norm/gamma)*SPEED_OF_LIGHT)/(rest_mass_energy*light_properties(7))
  if(in_parallel) then
#ifdef USE_TASKLOOP
    !$omp taskloop default(shared) private(ii,jj) &
    !$omp firstprivate(gamma,p_norm,rest_mass_energy,beta_c_over_omega) collapse(2)
#endif
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        call compute_bremsstrahlung_directionality_funct(spectra%points(jj,ii),&
        gamma,p_norm,rest_mass_energy,beta_c_over_omega,light_dstb(jj,ii))
      enddo
    enddo
#ifdef USE_TASKLOOP
    !$omp end taskloop
#endif
  else
    !$omp parallel do default(shared) private(ii,jj) &
    !$omp firstprivate(gamma,p_norm,rest_mass_energy,beta_c_over_omega) collapse(2)
    do ii=1,spectra%n_spectra
      do jj=1,spectra%n_points
        call compute_bremsstrahlung_directionality_funct(spectra%points(jj,ii),&
        gamma,p_norm,rest_mass_energy,beta_c_over_omega,light_dstb(jj,ii))
      enddo
    enddo
    !$omp end parallel do
  endif
end subroutine bremsstrahlung_directionality_funct

!> bremsstrahlung_spectral_irradiance computes the bremsstrahlung photon
!> spectral-angular distribution (photons per unit time, solid angle and
!> photon energy) emitted towards the shaded point x_shaded
!> inputs:
!>   light_vert: (bremsstrahlung_light_dist) bremsstrahlung light sources
!>   spectra:    (spectrum_base) spectral intervals and integrators
!>   time_id:    (integer) the time index
!>   light_id:   (integer) the light index
!>   x_shaded:   (real8)(3) shaded point position in cartesian coord
!> outputs:
!>   light_vert:            (bremsstrahlung_light_dist) bremsstrahlung light sources
!>   spectra:               (spectrum_base) spectral intervals and integrators
!>   light_spec_irradiance: (real8)(n_points,n_spectra) bremsstrahlung photon
!>                          spectral-angular distribution [1/(s.sr.J)] at x_shaded
subroutine bremsstrahlung_spectral_irradiance(light_vert,spectra,time_id,&
light_id,x_shaded,light_spec_irradiance)
  use mod_spectra,  only: spectrum_base
  implicit none
  !> inputs-outputs:
  class(bremsstrahlung_light_dist),intent(inout) :: light_vert
  class(spectrum_base),intent(inout)             :: spectra
  !> inputs:
  integer,intent(in)                          :: time_id,light_id
  real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
  !> outputs:
  real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance

  !> compute the directionality function
  call light_vert%directionality_funct(spectra,time_id,light_id,x_shaded,light_spec_irradiance)
  !> multiply by the light's pf prefactor (weight, local n_e, Zeff and QED constants)
  light_spec_irradiance = light_spec_irradiance*light_vert%properties(8,light_id,time_id)
end subroutine bremsstrahlung_spectral_irradiance

!> interpolate the JOREK background electron density required for
!> computing the bremsstrahlung radiation properties
!> inputs:
!>   light_vert:  (bremsstrahlung_light_dist) empty bremsstrahlung lights
!>   fields:      (fields_base) JOREK MHD fields
!>   particle_in: (particle_base) JOREK particle base structure
!>   time_id:     (integer) index of simulation time
!>   mass:        (real8) particle mass
!> outputs:
!>   mhd_fields: (real8)(n_mhd) background plasma fields at the particle location
!>               1 -> local electron density n_e [m^-3]
subroutine compute_bremsstrahlung_mhd_fields(light_vert,fields,&
particle_in,time_id,mass,mhd_fields)
  use mod_fields,         only: fields_base
  use mod_particle_types, only: particle_base
  !> used only for unit testing but required for compilation
  use mod_particle_common_test_tools, only: compute_test_n_e
  implicit none
  !> Inputs:
  class(bremsstrahlung_light_dist),intent(in) :: light_vert
  class(fields_base),intent(in)               :: fields
  class(particle_base),intent(in)             :: particle_in
  integer,intent(in)                          :: time_id
  real*8,intent(in)                           :: mass
  !> Outputs:
  real*8,dimension(light_vert%n_mhd),intent(out) :: mhd_fields
  !> Variables:
  real*8 :: T_e
#ifndef UNIT_TESTS_AFIELDS
  !> interpolate the local JOREK electron density
  call fields%calc_NeTe(light_vert%times(time_id),particle_in%i_elm,&
  particle_in%st,particle_in%x(3),mhd_fields(1),T_e)
#else
  !> analytical density only for unit testing
  call compute_test_n_e(particle_in%x,mhd_fields(1))
#endif
end subroutine compute_bremsstrahlung_mhd_fields

!> compute_bremsstrahlung_light_properties computes the
!> bremsstrahlung radiation properties from a kinetic
!> relativistic particle.
!> inputs:
!>   light_vert:  (bremsstrahlung_light_dist) empty bremsstrahlung lights
!>   property_id: (integer) index of the property to be initialised
!>   time_id:     (integer) time index
!>   particle_in: (particle_kinetic_relativistic) jorek particle
!>   mass:        (real8) mass of the particle [AMU]
!>   mhd_fields:  (real8)(n_mhd) background plasma fields
!>                1 -> local electron density n_e [m^-3]
!> outputs:
!>   light_vert: (bremsstrahlung_light_dist) bremsstrahlung lights with
!>               initialised properties. First dimension of the properties are:
!>                 1:3 -> electron velocity direction (cartesian) -> T = v/||v||
!>                 4   -> relativistic factor gamma = sqrt(1+(p/(mass*c))**2)
!>                 5   -> normalised electron momentum p = gamma*beta
!>                 6   -> electron rest mass energy [J] = mass*ATOMIC_MASS_UNIT*c**2
!>                 7   -> solid angle of the ~1/gamma emission cone [sr]
!>                 8   -> pf prefactor = weight*n_e*Zeff*EL_RAD**2*ALPHA_FINE_STRUCTURE [1/m]
subroutine compute_bremsstrahlung_light_properties(light_vert,&
property_id,time_id,particle_in,mass,mhd_fields)
  use constants,           only: ATOMIC_MASS_UNIT,SPEED_OF_LIGHT,TWOPI
  use constants,           only: EL_RAD,ALPHA_FINE_STRUCTURE
  use mod_particle_types,  only: particle_base,particle_kinetic_relativistic
  implicit none
  !> inputs-outputs
  class(bremsstrahlung_light_dist),intent(inout) :: light_vert
  !> inputs
  class(particle_base),intent(in)                :: particle_in
  integer,intent(in)                             :: property_id,time_id
  real*8,intent(in)                               :: mass
  real*8,dimension(light_vert%n_mhd),intent(in)   :: mhd_fields
  !> variables
  real*8 :: velocity,beta_gamma,gamma,beta,p_norm

  select type(p_in=>particle_in)
    type is (particle_kinetic_relativistic)
    !> compute velocity direction, relativistic factor and normalised momentum
    velocity   = norm2(p_in%p)
    beta_gamma = velocity/SPEED_OF_LIGHT
    gamma      = sqrt(1.d0 + (beta_gamma/mass)**2)
    beta       = beta_gamma/(mass*gamma)
    p_norm     = gamma*beta
    light_vert%properties(1:3,property_id,time_id) = p_in%p/velocity
    light_vert%properties(4,property_id,time_id)    = gamma
    light_vert%properties(5,property_id,time_id)    = p_norm
    !> electron rest mass energy [J]
    light_vert%properties(6,property_id,time_id)    = mass*ATOMIC_MASS_UNIT*(SPEED_OF_LIGHT**2)
    !> solid angle of the cone of half-angle theta_c=asin(1/gamma): 2*pi*(1-cos(theta_c))
    !> computed as 2*pi/(gamma*(gamma+p)) to avoid catastrophic cancellation for gamma>>1
    light_vert%properties(7,property_id,time_id)    = TWOPI/(gamma*(gamma+p_norm))
    !> pf prefactor (weight, local density, Zeff and QED constants)
    light_vert%properties(8,property_id,time_id)    = p_in%weight*mhd_fields(1)*&
    light_vert%Zeff*(EL_RAD**2)*ALPHA_FINE_STRUCTURE
  end select
end subroutine compute_bremsstrahlung_light_properties

!> initialise and allocate bremsstrahlung light variables
!> inputs:
!>   light_vert: (bremsstrahlung_light_dist) bremsstrahlung lights class
!> outputs:
!>   light_vert: (bremsstrahlung_light_dist) bremsstrahlung lights class
subroutine setup_bremsstrahlung_light_class(light_vert)
  use mod_particle_types, only: particle_kinetic_relativistic_id
  implicit none
  !> inputs-outputs
  class(bremsstrahlung_light_dist),intent(inout) :: light_vert
  !> set-up the bremsstrahlung light variables
  light_vert%n_property_vertex = 8; light_vert%n_mhd = 1;
  light_vert%n_particle_types = 1;
  if(allocated(light_vert%particle_types)) deallocate(light_vert%particle_types)
  allocate(light_vert%particle_types(light_vert%n_particle_types))
  light_vert%particle_types = [particle_kinetic_relativistic_id]
end subroutine setup_bremsstrahlung_light_class

!> Reconstruct the emitting electron from a bremsstrahlung light. Only
!> the electron's direction, relativistic momentum and macro-particle
!> weight can be recovered from the light properties (charge is fixed
!> to that of an electron by construction of this radiation model).
!> inputs:
!>   light_vert:   (bremsstrahlung_light_dist) bremsstrahlung lights class
!>   fields:       (fields_base) JOREK MHD fields data structure
!>   light_id:     (integer) index to the light to be treated
!>   time_id:      (integer) time index
!>   mass:         (real8) particle mass
!> outputs:
!>   particle_out: (particle_base) reconstructed particle
subroutine compute_particle_from_bremsstrahlung_light(&
light_vert,fields,light_id,time_id,mass,particle_out)
  use constants,                 only: TWOPI,SPEED_OF_LIGHT
  use constants,                 only: EL_RAD,ALPHA_FINE_STRUCTURE
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_base,particle_kinetic_relativistic
  !> used only for unit testing but required for compilation
  use mod_particle_common_test_tools, only: compute_test_n_e
  implicit none
  !> inputs:
  class(bremsstrahlung_light_dist),intent(in) :: light_vert
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: light_id,time_id
  real*8,intent(in)                           :: mass
  !> outputs:
  class(particle_base),intent(out) :: particle_out
  !> variables:
  integer :: ifail
  real*8  :: n_e,T_e
  select type (p_out => particle_out)
  type is (particle_kinetic_relativistic)
    !> initialise unknown variables, assuming an electron
    p_out%i_life = 0; p_out%t_birth = 0.0; p_out%q = int(-1,kind=1);
    !> compute spatial global and local coordinates
    p_out%x = cartesian_to_cylindrical(light_vert%x(:,light_id,time_id))
    if(p_out%x(3).lt.0d0) p_out%x(3) = TWOPI+p_out%x(3)
#ifndef UNIT_TESTS_AFIELDS
    call find_RZ(fields%node_list,fields%element_list,p_out%x(1),p_out%x(2),&
    p_out%x(1),p_out%x(2),p_out%i_elm,p_out%st(1),p_out%st(2),ifail)
    if(p_out%i_elm.le.0) return !< return if the particle is out-of-mesh
    call fields%calc_NeTe(light_vert%times(time_id),p_out%i_elm,p_out%st,&
    p_out%x(3),n_e,T_e)
#else
    call compute_test_n_e(p_out%x,n_e)
#endif
    !> reconstruct the relativistic momentum from direction and gamma*beta
    p_out%p = mass*SPEED_OF_LIGHT*light_vert%properties(5,light_id,time_id)*&
    light_vert%properties(1:3,light_id,time_id)
    !> reconstruct the macro-particle weight by inverting the pf prefactor
    p_out%weight = light_vert%properties(8,light_id,time_id)/&
    (n_e*light_vert%Zeff*(EL_RAD**2)*ALPHA_FINE_STRUCTURE)
  end select
end subroutine compute_particle_from_bremsstrahlung_light

!> Tools ------------------------------------------
!> compute the bremsstrahlung directionality function for a single
!> photon energy: the angle-integrated Born-approximation photon
!> spectrum (Hoppe et al 2018, Eq.9 / Koch & Motz 1959) divided by
!> the emission cone's solid angle and converted from a per-unit-
!> path-length rate to a per-unit-time rate (factor beta*c).
!> inputs:
!>   k_energy:           (real8) photon energy [J]
!>   gamma:              (real8) electron relativistic factor
!>   p_norm:             (real8) electron normalised momentum gamma*beta
!>   rest_mass_energy:   (real8) electron rest mass energy [J]
!>   beta_c_over_omega:  (real8) beta*c/(rest_mass_energy*Omega_cone) [m/(s.J.sr)]
!> outputs:
!>   dir_funct: (real8) bremsstrahlung directionality function [1/(s.sr.J)] per unit pf
subroutine compute_bremsstrahlung_directionality_funct(k_energy,&
gamma,p_norm,rest_mass_energy,beta_c_over_omega,dir_funct)
  implicit none
  !> inputs:
  real*8,intent(in)  :: k_energy,gamma,p_norm,rest_mass_energy,beta_c_over_omega
  !> outputs:
  real*8,intent(out) :: dir_funct
  !> variables:
  real*8 :: k_dim,gamma_p,p_p,eps,eps_p,L
  real*8 :: TT1,TT2,TT3,TT4,TT5,LT1,LT2,LTF,LT3,LT4,LT5,bracket
  !> initialisation
  dir_funct = 0.d0
  !> convert the photon energy to units of the electron rest mass energy
  k_dim = k_energy/rest_mass_energy
  !> kinematically forbidden or non-physical photon energy
  if(k_dim.le.0.d0) return
  if(k_dim.ge.(gamma-1.d0)) return
  !> final electron state after emitting the photon
  gamma_p = gamma-k_dim; p_p = sqrt(gamma_p**2-1.d0);
  eps   = 2.d0*log(gamma  +p_norm)
  eps_p = 2.d0*log(gamma_p+p_p)
  L     = 2.d0*log((gamma*gamma_p-1.d0+p_norm*p_p)/k_dim)
  !> Born-approximation bracket (Hoppe et al 2018, Eq.9)
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
  !> assemble the directionality function
  dir_funct = beta_c_over_omega*(p_p/(k_dim*p_norm))*bracket
end subroutine compute_bremsstrahlung_directionality_funct

!>-------------------------------------------------
end module mod_bremsstrahlung_light_dist_vertices
