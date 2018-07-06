#include"cppdefs.h"
!=======================================================================
MODULE EPBL_GOTM
!-----------------------------------------------------------------------
! Module Name: EPBL_GOTM
! Module Author: Brandon Reichl
! Module Data: May 2, 2018
!-----------------------------------------------------------------------
! Module Description
!-----------------------------------------------------------------------
! Contains the subroutines to call ePBL/JHL from GOTM
!  Options are contained in the namelists MOMturb.nml, epbl.nml, and
!  jhl.nml
!-----------------------------------------------------------------------
! Public subroutines
!   1. epbl_gotm_init
!   2. epbl_gotm_interface
!-----------------------------------------------------------------------
! External subroutines:
!   1. EPBL_aux/cvmix_epbl_interface_init
!   2. EPBL_aux/cvmix_epbl_interface
!     Sets control structures and interfaces to ePBL code.
use EPBL_aux, only: cvmix_epbl_interface_init
use EPBL_aux, only: cvmix_epbl_interface
!   3. EPBL_aux/cvmix_jhl_interface_init
!   4. EPBL_aux/cvmix_jhl_interface
!     Sets control structures and interfaces jhl code
use EPBL_aux, only: cvmix_jhl_interface_init
use EPBL_aux, only: cvmix_jhl_interface
!   5. EPBL_aux/cvmix_interface_init
!      Allows to initialize mean values (gravity and rho0)
use EPBL_aux, only: cvmix_interface_init
!   6. eqstate/eos_alpha
!     Computes alpha (where -alpha is the thermal expansion coefficient)
use eqstate, only: eos_alpha
!   7. eqstate/eos_beta
!     Computes beta ( where beta is the haline contraction coefficient)
use eqstate, only: eos_beta
!   8. eqstate/eqstate1
!     Computes the buoyancy
use eqstate, only: eqstate1
!
!-----------------------------------------------------------------------
! External parameters:
!  1. turbulence/num
!    Turbulent viscosity
use turbulence, only: num
!  2. turbulence/nuh
!    Turbulent diffusivity of heat
use turbulence, only: nuh
!  3. turbulence/nus
!    Turbulent diffusivity of salt
use turbulence, only: nus
!  4. turbulence/gamu/gamv/gamh/gams
!    Nonlocal tendency of u, v, h (heat), and s
use turbulence, only: gamu, gamv, gamh, gams
!  5. turbulence/rig
!    Gradient Richardson number
use turbulence, only: Rig
!  6. turbulence/tke
!    Turbulence Kinetic Energy
use turbulence, only: TKE
!  7. meanflow/gravity
!    Gravitational constant [m/s2]
use meanflow, only: gravity
!  8. meanflow/rho_0
!    Reference density [kg/m3]
use meanflow, only: rho_0
!  9. turbulence/nucl
!  turbulent eddy coefficient for momentum flux down Stokes gradient
!  in second moment closures with Craik-Leibovich vortex force in the
!  algebraic Reynolds stress and flux models.
use turbulence, only: nucl

!
use kpp, only: langmuir_method
use kpp, only: KPP_LT_NOLANGMUIR, KPP_LT_LWF16, KPP_LT_LF17
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
! MODE : Set the MOMturb Mode
!  1 - EPBL+JHL
!  2 - EPBL
!  3 - JHL
INTEGER :: MODE = 0
INTEGER :: MODE_JHL_EPBL = 0, MODE_EPBL = 1, MODE_JHL = 2,&
           MODE_JHL_EPBL_ADD = 3
!--------------------------------------------------------
! This should be a runtime parameter
REALTYPE :: MinLat = 10.0
! These are probably elsewhere in gotm? Should be set for consistency.
REALTYPE, parameter :: Deg2Rad = 3.141592653589793/180.
REALTYPE, parameter :: TwoOmega = 1.45842e-4
!
REALTYPE, public :: EPBL_OSBL
REALTYPE, public :: JHL_OSBL
!
!-----------------------------------------------------------------------
contains
!=======================================================================
subroutine epbl_gotm_init(nlev,namlst)
! EPBL_GOTM_INIT
! This subroutine intializes the turbulence fields and reads in
! namelists
!-----------------------------------------------------------------------
  integer, intent(in) :: nlev, namlst
  integer :: rc
  real :: rho0_4, grav_4 ! Used for converting from real(8) to real(4)
!
  namelist /momturb/ mode, langmuir_method
!
  allocate(num(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (num)'
  num = _ZERO_
  allocate(nuh(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (nuh)'
  nuh = _ZERO_
  allocate(nus(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (nus)'
  nus = _ZERO_
  allocate(nucl(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (nucl)'
  nucl = _ZERO_
  allocate(gamu(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (gamu)'
  gamu = _ZERO_
  allocate(gamv(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (gamv)'
  gamv = _ZERO_
  allocate(gamh(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (gamh)'
  gamh = _ZERO_
  allocate(gams(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (gams)'
  gams = _ZERO_
  allocate(Rig(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (Rig)'
  Rig = _ZERO_
  allocate(tke(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (tke)'
  tke = _ZERO_
!
# ifdef EXTRA_OUTPUT

  allocate(turb1(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (turb1)'
  turb1 = _ZERO_

  allocate(turb2(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (turb2)'
  turb2 = _ZERO_

  allocate(turb3(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (turb3)'
  turb3 = _ZERO_

  allocate(turb4(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (turb4)'
  turb4 = _ZERO_

  allocate(turb5(0:nlev),stat=rc)
  if (rc /= 0) stop 'init_turbulence: Error allocating (turb5)'
  turb5 = _ZERO_

# endif

  open(namlst,file='MOMturb.nml',status='old',action='read')
  read(namlst,nml=momturb)
  close(namlst)
  rho0_4=rho_0
  grav_4=gravity
  call cvmix_interface_init(rho0_=rho0_4,grav_=grav_4)
  call cvmix_epbl_interface_init(namlst)
  call cvmix_jhl_interface_init(namlst)

end subroutine epbl_gotm_init
!=======================================================================
subroutine epbl_gotm_interface( nlev,h,u,v,t,s,u_taus,u_taub,&
     tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f,dt)
! EPBL_GOTM_INTERFACE
! The main interface to GOTM
!-----------------------------------------------------------------------
integer, intent(in) :: nlev ! Number of model levels
!--
! It is assumed that the following arrays are z=0 at index nlev
! as is standard in GOTM.  EPBL is written z=0 at index 1, so these
! arrays must be flipped prior to calling the epbl routines.
! These are also converted to standard "real" type parameters
! from the gotom "REALTYPE" definition.
!--
  REALTYPE, intent(in) :: h(0:nlev),& ! Thicknesses of model levels
                          u(0:nlev),& ! Current (x)
                          v(0:nlev),& ! Current (y)
                          t(0:nlev),& ! Potential Temperature (deg C)
                          s(0:nlev),& ! Salinity (ppt)
                          u_taus,&    ! Surface friction velocity (m/s)
                          u_taub,&    ! Bottom friction velocity (m/s)
                          tFlux,&     ! Surface temperature flux (units?)
                          btFlux,&    ! Surface buoyancy flux due to temp (m2/s3)
                          sFlux,&     ! Surface salinity flux (units?)
                          bsFlux,&    ! Surface buoyancy flux due to salt (m2/s3)
                          tRad(0:nlev),& !Radiative temperature flux (units?)
                          bRad(0:nlev),& !Radiative buoyancy flux (units?)
                          f,&         ! Coriolis acceleration (1/s)
                          dt          ! Timestep (s)
!-----------------------------------------------------------------------
  integer :: ki, kigotm
  real, dimension(nlev) :: H_epbl, U_epbl, V_epbl, S_epbl, T_epbl
  real, dimension(nlev) :: dSV_dS, dSV_dT
  real, dimension(nlev+1) :: dbuoy_dt, dbuoy_ds
  real, dimension(nlev+1) :: NUH_shear, NUH_epbl, NUH_out
  real :: ustar_epbl, bflux_epbl, f_epbl, dt_epbl
  real :: MLD_out
  !-----------------------------------------------------------------------
  ! This loop is flipping the array and converting from REALTYPE to real
  do ki=1,nlev
     kigotm=nlev-ki+1
     H_epbl(ki) = H(kigotm)
     U_epbl(ki) = U(kigotm)
     V_epbl(ki) = V(kigotm)
     T_epbl(ki) = T(kigotm)
     S_epbl(ki) = S(kigotm)
      NUH_shear(ki) = NUH(kigotm)
  end do
  NUH_shear(nlev+1)=0.0
  ! These calls are converting from REALTYPE to real
  ustar_epbl = u_taus
  bflux_epbl = btflux+bsflux+(brad(nlev)-brad(nlev-1))! How to use brad?
  f_epbl = abs(f)
  dt_epbl = dt
  ! Initialize NUH_epbl, NUH_shear, and NUH_out
  NUH_epbl(:) = 0.0; NUH_out(:) = 0.0
  !-
  call compute_state_derivs(T_epbl,S_epbl,H_epbl,dSV_dT,dSV_dS,&
       dbuoy_dt,dbuoy_ds,nlev)
  call cvmix_epbl_interface(H_epbl,U_epbl,V_epbl,T_epbl,S_epbl, &
       dSV_dt, dSV_dS, &
       ustar_epbl, bflux_epbl,f_epbl,nlev,nlev+1, &
       NUH_epbl,MLD_out, dt_epbl)
  EPBL_OSBL = MLD_out

  call cvmix_jhl_interface(H_epbl,U_epbl,V_epbl,T_epbl,S_epbl, &
       dbuoy_dt, dbuoy_ds, ustar_epbl, bflux_epbl,f_epbl,nlev,nlev+1, &
       NUH_shear,MLD_out, dt_epbl)
  JHL_OSBL = MLD_out
  ! Flipping the viscosity/diffusivity back for GOTM and converting
  ! from real to REALTYPE
  do ki=1,nlev+1
     kigotm=nlev-ki+1
     if (MODE==MODE_JHL_EPBL) then
        NUH_out(ki) = max(NUH_shear(ki),NUH_epbl(ki))
     elseif (Mode==MODE_EPBL) then
        NUH_out(ki) = NUH_epbl(ki)
     elseif (Mode==MODE_JHL) then
        NUH_out(ki) = NUH_shear(ki)
     elseif (MODE==MODE_JHL_EPBL_ADD) then
        NUH_out(ki) = NUH_shear(ki)+NUH_epbl(ki)
     endif

     NUH(kigotm) = NUH_out(ki)
     NUS(kigotm) = NUH_out(ki)
     NUM(kigotm) = NUH_out(ki)
  end do
!
  return
!
end subroutine epbl_gotm_interface
!=======================================================================
subroutine compute_state_derivs(T,S,H,dSV_dT, dSV_dS, &
     dbuoy_dt,dbuoy_ds,npts)
! COMPUTE_STATE_DERIVS
! This subroutine is computing derivatives of the state variables
!-----------------------------------------------------------------------
! Potential temperature relative to the surface in C
  REAL,    intent(in),  dimension(:) :: T
  !Salinity in g/kg.
  REAL,    intent(in),  dimension(:) :: S
    !Salinity in g/kg.
  REAL,    intent(in),  dimension(:) :: H
  !The partial derivative of specific volume with salinity, in m3 kg-1 / (g/kg)
  real,    intent(out), dimension(:) :: dSV_dS
  !The partial derivative of specific volume with potential temperature, in m3 kg-1 K-1
  real,    intent(out), dimension(:) :: dSV_dT
    !The partial derivative of specific volume with salinity, in m3 kg-1 / (g/kg)
  real,    intent(out), dimension(:) :: dbuoy_dS
  !The partial derivative of specific volume with potential temperature, in m3 kg-1 K-1
  real,    intent(out), dimension(:) :: dbuoy_dT

  !The number of values to calculate.
  integer, intent(in)                :: npts
!
  REAL :: drho_dt, drho_ds
  REAL :: I_rho, buoy
  REAL :: P
  real :: ZtoP = 0.1 ! Approx conversion of water depth to dbar (equivalent to GOTM)
  integer :: j
  REALTYPE :: Tr8,Sr8,Pr8
!
  P=H(1)/2*ZtoP !Surface pressure approx 0
  do j=1,npts
     if (j>1) then
        P = P + 0.5*(H(j)+H(j-1))*ZtoP
     end if
     Sr8=S(j);Tr8=T(j);Pr8=P
     drho_dt = eos_alpha(Sr8,Tr8,Pr8,gravity,rho_0)
     drho_ds = eos_beta(Sr8,Tr8,Pr8,gravity,rho_0)
     buoy = eqstate1(Sr8,Tr8,Pr8,gravity,rho_0)
     I_rho = 1./(rho_0 * (_ONE_ - buoy/gravity))
     ! +/- signs are correct.
     dSV_dT(j) = dRho_dt * I_rho
     dSV_dS(j) = -dRho_dS * I_rho
  enddo
  P=0
  do j=1,npts+1
     if (j>1) then
        P = P + (H(j-1))*ZtoP
     end if
     Pr8=P
     if (j==1) then
        Sr8=S(j);Tr8=T(j);
     elseif (j>1 .and. j<npts+1) then
        Sr8=.5*(s(j)+s(j-1));Tr8=0.5*(T(j)+T(j-1))
     elseif (j==npts+1) then
        Sr8=s(npts);Tr8=T(npts)
     endif
     drho_dt = eos_alpha(Sr8,Tr8,Pr8,gravity,rho_0)
     drho_ds = eos_beta(Sr8,Tr8,Pr8,gravity,rho_0)
     dbuoy_dt(j) = drho_dt*gravity
     dbuoy_ds(j) = -drho_ds*gravity
  enddo
!
  return
!
end subroutine compute_state_derivs
!=======================================================================
end module EPBL_gotm
