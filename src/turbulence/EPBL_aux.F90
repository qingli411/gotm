module EPBL_aux
!-----------------------------------------------------------------------
! Module Name: EPBL_AUX
! Module Author: Brandon Reichl
! Module Data: May 2, 2018
!-----------------------------------------------------------------------
! Module Description
!-----------------------------------------------------------------------
! This module is meant to serve as the interface to cvmix_energetic_pbl.
!  Because cvmix_energetic_pbl contains control structures (which
!  cannot be dealt with by f2py), an additional interface is required
!  for use in the subroutines called using f2py.
! This module is also used as the interface for the GOTM experiments
!  so that efforts do not need to be duplicated.
!-----------------------------------------------------------------------
! External subroutines/Types
! 1. cvmix_energetic_pbl/cvmix_epbl_column
!    Routine to solve for ePBL mixing coefficients
! 2. cvmix_energetic_pbl/cvmix_energetic_PBL_CS
!    Control structure for cvmix_epbl_column
use cvmix_energetic_pbl, only: cvmix_energetic_PBL_CS, cvmix_epbl_column
! 3. cvmix_kappa_shear/cvmix_kappa_shear_column
!    Routine to solve for JHL08 mixing coefficients.
! 4. cvmix_kappa_shear/cvmix_kappa_shear_CS
!    Control structure for cvmix_kappa_shear_column
use cvmix_kappa_shear, only: cvmix_kappa_shear_CS, cvmix_kappa_shear_column
!
!-----------------------------------------------------------------------
implicit none ; private
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Public subroutines
public cvmix_interface_init
! 1. cvmix_interface_init
!    Optional initialization of some module variables
public cvmix_epbl_interface, cvmix_epbl_interface_init
! 2. cvmix_epbl_interface
!    Interface to be called externally to use cvmix_epbl_column
! 3. cvmix_epbl_interface_init
!    Routine to set control structures stored in this module
public cvmix_jhl_interface, cvmix_jhl_interface_init
! 4. cvmix_jhl_interface
!    Interface to be called externally to use cvmix_kappa_shear_column
! 5. cvmix_jhl_interface_init
!    Routine to set control structures stored in this module
!-----------------------------------------------------------------------
! CSepbl - Control structure target for ePBL
type(cvmix_energetic_PBL_CS), target, save :: CSepbl
! CSjhl - Control structure target for JHL08
type(cvmix_kappa_shear_CS), target, save :: CSjhl
!
! FirstJHL - logical for first call to allocate JHL storage arrays
logical :: FirstJHL = .true.
! TKEjhl, NUHjhl - Arrays to store TKE and NUH for JHL08 to be
!                  used in subsequent iterations
real, allocatable, dimension(:) :: TKEjhl, NUHjhl
! Some parameters (with reasonable initialized values that might
!   be overwritten by the calling models values later)
real :: rho0 = 1000.
real :: idtdr0
real :: grav = 9.80616
real :: m_to_h = 1.0
real :: h_to_m = 1.0
real :: h_neglect = 1.e-20 !??
!
!-----------------------------------------------------------------------
contains
!=======================================================================
subroutine cvmix_interface_init(rho0_,grav_)
! CVMIX_INTERFACE_INIT
! A way to set mean external parameters
  real, intent(in),optional :: rho0_ ! Mean Density from external model
                                     ! [kg/m3]
  real, intent(in),optional :: grav_ ! Gravitational acceleration from
                                     ! external model [m/s2]
  if (present(rho0_)) rho0 = rho0_
  if (present(grav_)) grav = grav_
end subroutine cvmix_interface_init
!=======================================================================
subroutine cvmix_epbl_interface_init(namlst)
! CVMIX_EPBL_INTERFACE_INIT
! This subroutine initializes the control structure based on a
!  namelist.  It should be equivalent to MOM6 read formats
!  except parameter names in the namelist are identical to those
!  in the control structures.
!-----------------------------------------------------------------------
  integer, intent(in) :: namlst ! NAMLST integer (opened and closed
                                ! within this routine, so reusable).
!
! Local parameters: see the cvmix_energetic_pbl_CS for descriptions
!
  integer :: mstar_mode ,lt_enhance_form, vstar_mode
  real :: mstar, nstar, mixlenexponent, tke_decay,mke_to_tke_effic,&
       wstar_ustar_coef,vstar_scale_fac,ekman_scale_coef,&
       translay_scale,mld_tol,min_mix_len,n2_dissipation_scale_neg,&
       n2_dissipation_scale_pos,mstar_cap,mstar_slope,mstar_xint,&
       mstar_at_xint,&
       lt_enhance_coef, lt_enhance_exp, mstar_n, c_ek, mstar_coef,&
       lac_mldoek, lac_mldoob_stab, lac_ekoob_stab,&
       lac_mldoob_un, lac_ekoob_un, ladepthratio,&
       max_enhance_m, cnv_mst_fac,&
       vstar_surf_fac,RH18_CN1,RH18_CN2,RH18_CN3,RH18_CS1,RH18_CS2
  logical :: Use_LA_WindSea
!
  namelist /epbl/ mstar_mode, mstar, nstar, mixlenexponent, &
       tke_decay, mke_to_tke_effic, wstar_ustar_coef, &
       vstar_mode,vstar_surf_fac,vstar_scale_fac, ekman_scale_coef, &
       translay_scale, mld_tol, min_mix_len, &
       n2_dissipation_scale_neg, n2_dissipation_scale_pos, &
       mstar_cap, mstar_slope, mstar_xint, mstar_at_xint,&
       rh18_cn1, rh18_cn2, rh18_cn3, rh18_cs1, rh18_cs2,&
       lt_enhance_coef, &
       lt_enhance_exp, mstar_n, c_ek, mstar_coef, &
       Use_LA_WindSea, lac_mldoek, lac_mldoob_stab, &
       lac_ekoob_stab, lac_mldoob_un, lac_ekoob_un, &
       ladepthratio, max_enhance_m, cnv_mst_fac,lt_enhance_form
!
!-----------------------------------------------------------------------
!
! 1. Read in input file namelist
!
  open(20,file='epbl.nml',status='old',action='read')
  read(20,nml=epbl)
  close(20)
!
! 2. Fill ePBL control structure from namelist input
!
  CSepbl%mstar_mode = mstar_mode
  CSepbl%vstar_mode = vstar_mode
  CSepbl%mstar = mstar
  CSepbl%nstar = nstar
  CSepbl%mixlenexponent=mixlenexponent
  CSepbl%tke_decay=tke_decay
  CSepbl%mke_to_tke_effic=mke_to_tke_effic
  CSepbl%wstar_ustar_coef=wstar_ustar_coef
  CSepbl%vstar_scale_fac=vstar_scale_fac
  CSepbl%vstar_surf_fac=vstar_surf_fac
  CSepbl%ekman_scale_coef=ekman_scale_coef
  CSepbl%translay_scale=translay_scale
  CSepbl%mld_tol=mld_tol
  CSepbl%min_mix_len=min_mix_len
  CSepbl%n2_dissipation_scale_neg=n2_dissipation_scale_neg
  CSepbl%n2_dissipation_scale_pos=n2_dissipation_scale_pos
  CSepbl%mstar_cap=mstar_cap
  CSepbl%mstar_slope=mstar_slope
  CSepbl%mstar_xint=mstar_xint
  CSepbl%mstar_at_xint=mstar_at_xint
  CSepbl%RH18_CN1 = RH18_CN1
  CSepbl%RH18_CN2 = RH18_CN2
  CSepbl%RH18_CN3 = RH18_CN3
  CSepbl%RH18_CS1 = RH18_CS1
  CSepbl%RH18_CS2 = RH18_CS2
  CSepbl%lt_enhance_coef=lt_enhance_coef
  CSepbl%lt_enhance_exp=lt_enhance_exp
  CSepbl%mstar_n=mstar_n
  CSepbl%c_ek=c_ek
  CSepbl%mstar_coef=mstar_coef
  CSepbl%Use_LA_WindSea=Use_LA_WindSea
  CSepbl%lac_mldoek=lac_mldoek
  CSepbl%lac_mldoob_stab=lac_mldoob_stab
  CSepbl%lac_ekoob_stab=lac_ekoob_stab
  CSepbl%lac_mldoob_un=lac_mldoob_un
  CSepbl%lac_ekoob_un=lac_ekoob_un
  CSepbl%ladepthratio=ladepthratio
  CSepbl%max_enhance_m=max_enhance_m
  CSepbl%cnv_mst_fac=cnv_mst_fac
  CSepbl%lt_enhance_form=lt_enhance_form
!
! 3. Fitting some coefficients based on empirical models
!   a. Fitting coefficients to asymptote twoard 0 as MLD -> Ekman depth
!
  CSepbl%MSTAR_A = CSepbl%MSTAR_AT_XINT**(1./CSepbl%MSTAR_N)
  CSepbl%MSTAR_B = CSepbl%MSTAR_SLOPE / (CSepbl%MSTAR_N* &
                       CSepbl%MSTAR_A**(CSepbl%MSTAR_N-1.))
!
!   b. Fitting coefficients to asymptote toward MSTAR_CAP
!      *Fixed to begin asymptote at MSTAR_CAP-0.5 toward MSTAR_CAP
  CSepbl%MSTAR_A2 = 0.5**(1./CSepbl%MSTAR_N)
  CSepbl%MSTAR_B2 = -CSepbl%MSTAR_SLOPE / (CSepbl%MSTAR_N* &
                         CSepbl%MSTAR_A2**(CSepbl%MSTAR_N-1))
!
!   c. Compute value of X (referenced to MSTAR_XINT) where transition
!      to asymptotic regime based on value of X where MSTAR=MSTAR_CAP-0.5
!
  CSepbl%MSTAR_XINT_UP = (CSepbl%MSTAR_CAP-0.5-&
                              CSepbl%MSTAR_AT_XINT)/CSepbl%MSTAR_SLOPE
!
  return
!
endsubroutine cvmix_epbl_interface_init
!=======================================================================
subroutine cvmix_jhl_interface_init(namlst)
! CVMIX_JHL_INTERFACE_INIT
! This subroutine initializes the control structure based on a
!  namelist.  It should be equivalent to MOM6 read formats
!  except parameter names in the namelist are identical to those
!  in the control structures.
!-----------------------------------------------------------------------
  integer, intent(in) :: namlst ! NAMLST integer (opened and closed
                                ! within this routine, so reusable).
!
! Local parameters: see the cvmix_kappa_shear_CS for descriptions
!
  real :: RiNo_Crit, ShearMix_Rate, FRi_Curvature, &
       C_N, C_S, Lambda, Lambda2_N_S, TKE_bg, &
       Kappa_0, Kappa_Tol_Err, Prandtl_Turb
  integer :: Max_RiNo_It, Max_KS_it
!
  namelist /jhl/ RiNo_Crit, ShearMix_Rate, FRi_Curvature, &
       C_N, C_S, Lambda, Lambda2_N_S, TKE_bg, &
       Kappa_0, Kappa_Tol_Err, Prandtl_Turb, &
       Max_RiNo_It, Max_KS_It
!-----------------------------------------------------------------------
!
! 1. Open and read namelist
!
  open(namlst,file='jhl.nml',status='old',action='read')
  read(namlst,nml=jhl)
  close(namlst)
!
! 2. Set control structure parameters
!
  CSjhl%RiNo_Crit = RiNo_Crit
  CSjhl%ShearMix_Rate = ShearMix_Rate
  CSjhl%FRi_Curvature = FRi_Curvature
  CSjhl%C_N = C_N
  CSjhl%C_S = C_S
  CSjhl%Lambda = Lambda
  CSjhl%Lambda2_N_S = Lambda2_N_S
  CSjhl%TKE_bg = TKE_bg
  CSjhl%Kappa_0 = Kappa_0
  CSjhl%Kappa_Tol_Err = Kappa_Tol_Err
  CSjhl%Prandtl_Turb = Prandtl_Turb
  CSjhl%Max_Rino_It = Max_Rino_It
  CSjhl%Max_KS_It = Max_KS_It
!
return
!
endsubroutine cvmix_jhl_interface_init
!=======================================================================
subroutine cvmix_epbl_interface(H, U, V, T, S, dSV_dT, dSV_dS, ustar, BF, cori, NZ, NZZ, NUH, BLD, DT)
! CVMIX_EPBL_INTERFACE
! Calls the EPBL routines to compute mixing coefficient profiles.
!-----------------------------------------------------------------------
  integer, intent(in) :: NZ,& ! Number of grid cells (tracer points)
                         NZZ  ! Numer of grid interfaces (diffusivity points)
  real, intent(in), dimension(NZ) :: H,& ! Grid thickness [m]
                                     U,& ! Velocity (x-dir) [m/s]
                                     V,& ! Velocity (y-dir) [m/s]
                                     T,& ! Temperatue [deg C]
                                     S,& ! Salinity [ppt]
                                     dSV_dT,& ! Change of Specific volume
                                              ! wrt Temp [m3/kg/degC]
                                     dSV_dS ! Change of Specific volume
                                            ! wrt Salt [m3/kg/ppt]
  real, intent(in) :: ustar,& ! Surface friction coefficient [m/s]
                      BF,& ! Surface buoyancy flux [m2/s3]
                      CORI ! Coriolis frequency [1/s]
  real, intent(in) :: DT ! Time step (s)
  real, intent(inout), dimension(NZZ) :: NUH ! Eddy diffusivity on interfaces
                                             ! [m2/s]
  real, intent(out) :: BLD ! Boundary Layer Depth [m]
!
  real, dimension(NZ) :: P ! Pressure (referenced to surface of 0)
  real, dimension(NZ) :: tke_forced ! TKE supply to EPBL (taken as 0)
!
  type(cvmix_energetic_pbl_cs), pointer :: CSp ! Pointer to ePBL CS
!
  real :: h_to_kg_m2 ! if non-Boussinesq (not complete)
  real :: mech_tke ! Surface mech_tke (only used if mstar_mode==0)
  real :: conv_PE_released ! PE released at surface by convection
  real :: ust_use ! A maximum between ustar and a small value
  integer :: zi ! vertical index
  logical :: boussinesq ! Don't make this false
!
! 1. Point to the target control structure
!
  CSp => CSepbl
!
! 2. Set some parameters
!
  boussinesq = .true.
  idtdr0 = 1./(dt*rho0)
  mech_tke = 0.0
  conv_PE_released = 0.0
!
  if (boussinesq) then
     h_to_kg_m2 = rho0 * h_to_m !Needs modified if not Boussinesq
  else
     print*,'Stopping, boussinesq flag false and code not ready'
     stop
  end if
!
  do zi=1,NZ
     TKE_forced(zi) = 0.0
  enddo
!
  call calculate_cTKE(NZ,TKE_forced)
!
  if (CSp%mstar_mode==0) then
     mech_tke = (dt*CSp%mstar*Rho0)*((UStar**3))
  endif
!
  ust_use = max(CSp%ustar_min,ustar)
!
! 3. Actual call to cvmix_epbl_column
!
  call cvmix_epbl_column(NZ,NUH,H,U,V,T,S,&
                         dSV_dT, dSV_dS, TKE_forced,dt, &
                         idtdr0, m_to_h, h_to_m, h_to_kg_m2, grav, rho0, h_neglect, &
                         BF, ust_use, ust_use, cori, mech_tke, conv_PE_released, &
                         CSp, BLD)
!
  return
!
end subroutine cvmix_epbl_interface
!=======================================================================
subroutine cvmix_jhl_interface(H, Ui, Vi, Ti, Si, dbuoy_dt, dbuoy_ds, &
     ustar, BF, cori, NZ, NZZ, &
     NUH, BLD, DT)
! CVMIX_JHL_INTERFACE
! Calls the JHL routine to compute mixing coefficients
  integer, intent(in) :: NZ,& ! Number of grid cells (tracer points)
                         NZZ  ! Numer of grid interfaces (diffusivity points)
  real, intent(in), dimension(NZ) :: H,& ! Grid thickness [m]
                                     Ui,& ! Velocity (x-dir) [m/s]
                                     Vi,& ! Velocity (y-dir) [m/s]
                                     Ti,& ! Temperatue [deg C]
                                     Si   ! Salinity [ppt]
  ! Note the following two are on interfaces because they are used w/ N2
  real, intent(in), dimension(NZZ) :: dbuoy_dT,& ! Change of buoyancy
                                                ! wrt Temp [m/s2/degC]
                                     dbuoy_dS   ! Change of Specific volume
                                                ! wrt Salt [m/s2/ppt]
  real, intent(in) :: ustar,& ! Surface friction coefficient [m/s]
                      BF,& ! Surface buoyancy flux [m2/s3]
                      CORI ! Coriolis frequency [1/s]
  real, intent(in) :: DT ! Time step (s)
  real, intent(inout), dimension(NZZ) :: NUH ! Eddy diffusivity on interfaces
                                             ! [m2/s]
  real, intent(out) :: BLD ! Boundary Layer Depth [m]
!
  real, dimension(NZZ) :: TKE, Kappa_avg, TKE_avg
  real, dimension(NZ)  :: U, V, S, T
!
  type(cvmix_kappa_shear_cs), pointer :: CSp
!
  integer :: zi
!
! 1. Setting some values
!
  U(:) = Ui
  V(:) = Vi
  T(:) = Ti
  S(:) = Si
  CSp => CSjhl
  if (firstjhl) then
     firstjhl=.false.
     allocate(NUHjhl(NZ+1),TKEjhl(NZ+1))
     NUHjhl(:)=1.0;TKEjhl(:)=0.0
  endif
!
! 2. Calling JHL
!
  call cvmix_kappa_shear_column(NZ,NZ,dt,cori**2,.true.,H,  &
                         dbuoy_dt,dbuoy_ds,                 &
                         U,V,T,S,                           &
                         NUHjhl, TKEjhl, Kappa_avg, TKE_avg,&
                         CSp)
!
! 3. Set output values and store values for next timestep
!    to make iteration converge faster
!
  TKE=TKE_avg
  NUH=Kappa_avg
  NUHjhl=NUH
  TKEjhl=TKE
!
!  4. How to set JHL BLD?
!
  BLD=0.0
!
  return
!
end subroutine cvmix_jhl_interface
!=======================================================================
subroutine calculate_cTKE(NZ,cTKE)
  integer, intent(in) :: NZ
  real, intent(out), dimension(NZ) :: cTKE
! Just set cTKE to zero.
cTKE(:) = 0.0
!
return
!
end subroutine calculate_cTKE
!=======================================================================
end module EPBL_aux
