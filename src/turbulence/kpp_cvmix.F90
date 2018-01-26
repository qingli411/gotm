#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: kpp_cvmix: the KPP-turbulence model with CVMix \label{sec:kpp_cvmix}
!
! !INTERFACE:
   module kpp_cvmix
!
! !DESCRIPTION:
!  Routines to do KPP using CVMix routines

! !USES:
   use turbulence,   only: num,nuh,nus
   use turbulence,   only: gamu,gamv,gamh,gams
   use turbulence,   only: Rig
   use turbulence,   only: kappa
!  use variables in module turbulence
!  Qing Li, 20171120
   use turbulence,   only: tke,tkeo,eps,L,kb,epsb,P,B,Pb,gamb,cmue1,cmue2
   use turbulence,   only: gam,an,as,at,r,xRf,uu,vv,ww
!  use variables in module airsea
!  Qing Li, 20171214
   use airsea,       only: u10, v10
   use meanflow,     only: pi, gravity
   use observations, only: nfreq, wav_freq, wav_spec, wav_xcmp, wav_ycmp

#ifdef EXTRA_OUTPUT
   use turbulence,   only: turb1,turb2,turb3,turb4,turb5
#endif

!  use CVMix
   use cvmix_kinds_and_types, only : cvmix_data_type,          &
                                     cvmix_global_params_type
   use cvmix_kpp,             only : cvmix_init_kpp,                           &
                                     cvmix_get_kpp_real,                       &
                                     cvmix_kpp_compute_OBL_depth,              &
                                     cvmix_kpp_compute_kOBL_depth,             &
                                     cvmix_kpp_compute_bulk_Richardson,        &
                                     cvmix_kpp_compute_unresolved_shear,       &
                                     cvmix_kpp_compute_turbulent_scales,       &
                                     cvmix_kpp_compute_shape_function_coeffs,  &
                                     cvmix_kpp_efactor_read,                   &
                                     cvmix_kpp_efactor_model,                  &
                                     cvmix_coeffs_kpp
   use cvmix_put_get,         only : cvmix_put

   IMPLICIT NONE

   private

! !PUBLIC MEMBER FUNCTIONS:
!
   public init_kpp_cvmix, do_kpp_cvmix, clean_kpp_cvmix

! !PUBLIC DATA MEMBERS:
!
!  z-position of surface boundary layer depth
   REALTYPE, public                 :: zsbl

!  z-position of bottom boundary layer depth
   REALTYPE, public                 :: zbbl

! !DEFINED PARAMETERS:
!
!  non-dimensional extent of the surface layer (epsilon=0.1)
   REALTYPE, parameter :: epsilon = 0.1

!  method of Langmuir turbulence parameterization
!  Qing Li, 20171213
   integer, parameter ::  KPP_LT_NOLANGMUIR = 0
   integer, parameter ::  KPP_LT_EFACTOR_MODEL = 1
   integer, parameter ::  KPP_LT_EFACTOR_READ = 2
   integer, parameter ::  KPP_LT_EFACTOR_SPEC = 3
   integer, parameter ::  KPP_LT_ENTRAINMENT = 4

! !LOCAL VARIABLES:
!
!  proportionality coefficient for
!  parameterizing non-local transport
   REALTYPE                              :: Cg

!  coefficient from computation of
!  turbulent shear velocity
   REALTYPE                              :: Vtc

!  acceleration of gravity
   REALTYPE                              :: g

!  reference density
   REALTYPE                              :: rho_0

!  g/rho_0
   REALTYPE                              :: gorho0

!  critical bulk Richardson number
   REALTYPE                              :: Ric

!  compute surface and bottom BBL
   logical                               :: kpp_sbl,kpp_bbl

!  compute internal mixing
   logical                               :: kpp_interior

!  use clipping of MLD at Ekman and Monin-Oboukhov scale
   logical                               :: clip_mld

!  positions of grid faces and centers
   REALTYPE, dimension(:), allocatable   :: z_w,z_r

!  distance between centers
   REALTYPE, dimension(:), allocatable   :: h_r

   integer                               :: ksblOld
   REALTYPE                              :: zsblOld

!  method to parameterize the effects of Langmuir turbulence
!  Qing Li, 20171213
   integer                               :: langmuir_method
   character(len=PATH_MAX)               :: langmuir_file

!  CVMix datatypes
   type(cvmix_data_type)                 :: CVmix_vars
   type(cvmix_global_params_type)        :: CVmix_params
!
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the CVMix KPP module
!
! !INTERFACE:
   subroutine init_kpp_cvmix(namlst,fn,nlev,h0,h,kpp_g,kpp_rho_0)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  namelist reference
   integer,          intent(in)        :: namlst

!  filename containing namelist
   character(len=*), intent(in)        :: fn

!  number of grid cells
   integer,          intent(in)        :: nlev

!  bathymetry (m)
   REALTYPE,         intent(in)        :: h0

!  size of grid cells (m)
   REALTYPE,         intent(in)        :: h(0:nlev)

!  acceleration of gravity (m/s^2)
   REALTYPE,         intent(in)        :: kpp_g

!  reference density (kg/m^3)
   REALTYPE,         intent(in)        :: kpp_rho_0
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: k
   integer                             :: rc

   namelist /kpp_cvmix/                kpp_sbl,kpp_bbl,kpp_interior,    &
                                       clip_mld,Ric,langmuir_method,    &
                                       langmuir_file
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_kpp_cvmix...'

   ! read the variables from the namelist file
   open(namlst,file=fn,status='old',action='read',err=80)

   LEVEL2 'reading kpp_cvmix namelist...'

   read(namlst,nml=kpp_cvmix,err=81)
   close (namlst)

   LEVEL2 'done.'

!  allocate memory for variables defined in other modules
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

   allocate(z_w(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (z_w)'
   z_w = _ZERO_

   allocate(z_r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (z_r)'
   z_r = _ZERO_

   allocate(h_r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (h_r)'
   h_r = _ZERO_

!  allocate memory defined in module turbulence
!  Qing Li, 20171120

   LEVEL2 'allocation memory..'
   allocate(tke(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (tke)'
   tke = _ZERO_

   allocate(tkeo(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (tkeo)'
   tkeo = _ZERO_

   allocate(eps(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (eps)'
   eps = _ZERO_

   allocate(L(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (L)'
   L = _ZERO_

   LEVEL2 'allocation memory..'
   allocate(kb(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (kb)'
   kb = _ZERO_

   LEVEL2 'allocation memory..'
   allocate(epsb(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (epsb)'
   epsb = _ZERO_

   allocate(P(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (P)'
   P = _ZERO_

   allocate(B(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (B)'
   B = _ZERO_

   allocate(Pb(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (Pb)'
   Pb = _ZERO_

   allocate(gamb(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamb)'
   gamb = _ZERO_

   allocate(cmue1(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (cmue1)'
   cmue1 = _ZERO_

   allocate(cmue2(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (cmue2)'
   cmue2 = _ZERO_

   allocate(gam(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gam)'
   gam = _ZERO_

   allocate(an(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (an)'
   an = _ZERO_

   allocate(as(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (as)'
   as = _ZERO_

   allocate(at(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (at)'
   at = _ZERO_

   allocate(r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (r)'
   r = _ZERO_

   allocate(xRf(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (xRf)'
   xRf = _ZERO_

   allocate(uu(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (uu)'
   uu = _ZERO_

   allocate(vv(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (vv)'
   vv = _ZERO_

   allocate(ww(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (ww)'
   ww = _ZERO_

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


!  report model parameters

   LEVEL2 '--------------------------------------------------------'
   LEVEL3 'You are using the KPP turbulence model with CVMix       '
   LEVEL3 'with the following specifications:                      '
   LEVEL3 '                                                        '
   if (kpp_interior) then
      LEVEL3 'Interior mixing algorithm                  - active -   '
# ifdef KPP_SHEAR
      LEVEL4 'KPP shear instability mixing           - active -   '
# else
      LEVEL4 'KPP shear instability mixing       - not active -   '
# endif
# ifdef KPP_INTERNAL_WAVE
      LEVEL4 'KPP internal wave mixing               - active -   '
# else
      LEVEL4 'KPP internal wave mixing           - not active -   '
# endif
# ifdef KPP_CONVEC
      LEVEL4 'KPP convective instability mixing      - active -   '
# else
      LEVEL4 'KPP convective instability mixing  - not active -   '
# endif
# ifdef KPP_DDMIX
      LEVEL4 'KPP double diffusive mixing            - active -   '
# else
      LEVEL4 'KPP double diffusive mixing        - not active -   '
# endif
   else
      LEVEL3 'Interior mixing algorithm              - not active -   '
   endif

   if (kpp_sbl) then
      LEVEL3 'Surface layer mixing algorithm             - active -   '
      if (clip_mld) then
         LEVEL4 'Clipping at Ekman/Oboukhov scale       - active -   '
      else
         LEVEL4 'Clipping at Ekman/Oboukhov scale   - not active -   '
      endif
# ifdef KPP_SALINITY
      LEVEL4 'Compute salinity fluxes                - active -   '
# else
      LEVEL4 'Compute salinity fluxes            - not active -   '
# endif
# ifdef NONLOCAL
      LEVEL4 'Nonlocal fluxes                        - active -   '
# else
      LEVEL4 'Nonlocal fluxes                    - not active -   '
# endif
# ifdef KPP_TWOPOINT_REF
      LEVEL4 'Ri_b from 2-point interpolation        - active -   '
# else
      LEVEL4 'Ri_b from 2-point interpolation    - not active -   '
# endif
# ifdef KPP_IP_FC
      LEVEL4 'F_c =0 criterion for SL-depth          - active -   '
# else
      LEVEL4 'Ri_b - Ri_c =0 criterion for SL-depth  - active -   '
# endif
# ifdef KPP_CLIP_GS
      LEVEL4 'Clipping G''(sigma) for matching        - active -   '
# else
      LEVEL4 'Clipping G''(sigma) for matching    - not active -   '
# endif

      LEVEL4 'Ri_c  = ', Ric

   else
      LEVEL3 'Surface layer mixing algorithm         - not active -   '
   endif

   if (kpp_bbl) then
      LEVEL3 'Bottom layer mixing algorithm              - active -   '
      LEVEL4 '(Same parameters as surface layer mixing)'

   else
      LEVEL3 'Bottom layer mixing algorithm          - not active -   '
   endif

!  message of Langmuir turbulence parameterization
!  Qing Li, 20171213
   if (langmuir_method .eq. KPP_LT_NOLANGMUIR) then
      LEVEL3 'Langmuir turbulence parameterization   - not active -   '
   else
      LEVEL3 'Langmuir turbulence parameterization       - active -   '
      if (langmuir_method .eq. KPP_LT_EFACTOR_MODEL) then
         LEVEL4 'Approximate enhancement factor from simple model'
      else if (langmuir_method .eq. KPP_LT_EFACTOR_READ) then
         LEVEL4 'Read enhancement factor from file'
      else if (langmuir_method .eq. KPP_LT_EFACTOR_SPEC) then
         LEVEL4 'Calculate enhancement factor from wave spectrum'
      else if (langmuir_method .eq. KPP_LT_ENTRAINMENT) then
         LEVEL4 'Langmuir mixing + Langmuir enhanced entrainment'
      else
         LEVEL4 'Method not supported'
      end if
   end if

   LEVEL2 '--------------------------------------------------------'

   ! Initialize parameter datatype and set up column
   call cvmix_init_kpp(ri_crit=Ric,                                &
                       ! interp_type=interp_type,                    &
                       ! lEkman=lcheckekmo,                          &
                       ! lMonOb=lcheckekmo,                          &
                       ! llangmuirEF=llangmuir_efactor,              &
                       MatchTechnique='MatchGradient',             &
                       ! lnoDGat1=lnoDGat1,                          &
                       surf_layer_ext = epsilon)
   call cvmix_put(CVmix_vars, 'nlev', nlev)
   call cvmix_put(CVmix_vars, 'max_nlev', nlev)
   call cvmix_put(CVmix_vars, 'ocn_depth', h0)

   LEVEL1 'done.'

   return

80 FATAL 'I could not open "kpp_cvmix.nml"'
   stop 'init_kpp_cvmix'
81 FATAL 'I could not read "kpp_cvmix" namelist'
   stop 'init_kpp_cvmix'
!
 end subroutine init_kpp_cvmix
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do KPP with CVMix
!
! !INTERFACE:
   subroutine do_kpp_cvmix(nlev,h0,h,rho,u,v,NN,NNT,NNS,SS,  &
                u_taus,u_taub,tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  number of grid cells
   integer                                       :: nlev

!  bathymetry (m)
   REALTYPE                                      :: h0

!  thickness of grid cells (m)
   REALTYPE                                      :: h(0:nlev)

!  potential density at grid centers (kg/m^3)
   REALTYPE                                      :: rho(0:nlev)

!  velocity components at grid centers (m/s)
   REALTYPE                                      :: u(0:nlev),v(0:nlev)

!  square of buoyancy frequency (1/s^2)
   REALTYPE                                      :: NN(0:nlev)

!  square of buoyancy frequency caused by
!  temperature and salinity stratification
   REALTYPE                                      :: NNT(0:nlev),NNS(0:nlev)

!  square of shear frequency (1/s^2)
   REALTYPE                                      :: SS(0:nlev)

!  surface and bottom friction velocities (m/s)
   REALTYPE                                      :: u_taus,u_taub

!  surface temperature flux (K m/s) and
!  salinity flux (psu m/s) (negative for loss)
   REALTYPE                                      :: tFlux,sFlux

!  surface buoyancy fluxes (m^2/s^3) due to
!  heat and salinity fluxes
   REALTYPE                                      :: btFlux,bsFlux

!  radiative flux [ I(z)/(rho Cp) ] (K m/s)
!  and associated buoyancy flux (m^2/s^3)
   REALTYPE                                      :: tRad(0:nlev),bRad(0:nlev)

!  Coriolis parameter (rad/s)
   REALTYPE                                      :: f
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!BOC
   call cvmix_put(CVmix_vars, 'Coriolis', f)
   continue
!
   end subroutine do_kpp_cvmix
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Clean up the CVMix KPP module
!
! !INTERFACE:
   subroutine clean_kpp_cvmix()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!BOC
!
   end subroutine clean_kpp_cvmix
!EOC

 end module kpp_cvmix

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

