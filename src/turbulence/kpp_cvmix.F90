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
   call cvmix_put(CVmix_params, 'Gravity', kpp_g)

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
!
!-----------------------------------------------------------------------
! Update model grid
!-----------------------------------------------------------------------

!  Compute distance between centers (between rho-points)
!  Note that h is the distance between faces (between w-points)
   do k=1,nlev-1
      h_r(k) = 0.5*(h(k)+ h(k+1))
   enddo

!  Compute position of interfaces (w-points)
   z_w(0) = - h0
   do k=1,nlev
      z_w(k) = z_w(k-1) + h(k)
   enddo

!  Compute position of centers (rho-points)
   z_r(1) = - h0 + 0.5*h(1)
   do k=2,nlev
      z_r(k) = z_r(k-1) + h_r(k-1)
   enddo

   CVmix_vars%zw_iface => z_w(0:nlev)
   CVmix_vars%zt_cntr  => z_r(1:nlev)

!-----------------------------------------------------------------------
! Update Cvmix variables
!-----------------------------------------------------------------------

   call cvmix_put(CVmix_vars, 'Coriolis', f)


!-----------------------------------------------------------------------
! compute interior mixing
!-----------------------------------------------------------------------

   if (kpp_interior) then
      call interior(nlev,NN,NNT,NNS,SS)
   else
      num = _ZERO_
      nuh = _ZERO_
      nus = _ZERO_
   endif

!-----------------------------------------------------------------------
! compute surface boundary layer mixing
!-----------------------------------------------------------------------

   if (kpp_sbl) then
      call surface_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,            &
                         tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
   endif

!-----------------------------------------------------------------------
! compute bottom boundary layer mixing
!-----------------------------------------------------------------------

   if (kpp_bbl) then
      call bottom_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,             &
                        _ZERO_,_ZERO_,_ZERO_,_ZERO_,tRad,bRad,f)
   endif

!
 end subroutine do_kpp_cvmix
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute turbulence in the surface layer \label{sec:kppSurface}
!
! !INTERFACE:
   subroutine surface_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,       &
                            tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
!
! !DESCRIPTION:
! In this routine all computations related to turbulence in the surface layer
! are performed. The algorithms are described in \sect{sec:kpp}. Note that these
! algorithms are affected by some pre-processor macros defined in {\tt cppdefs.inp},
! and by the parameters set in {\tt kpp.nml}, see \sect{sec:kpp}.
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

!  surface and bottom friction velocities (m/s)
   REALTYPE                                      :: u_taus,u_taub

!  surface temperature flux (K m/s) and
!  salinity flux (sal m/s) (negative for loss)
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
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE, parameter          :: eps      = 1.0E-10

   integer                      :: k,ksbl
   integer                      :: kk,kref

   REALTYPE                     :: Bo
   REALTYPE                     :: Bfsfc
   REALTYPE                     :: tRadSrf,tRadSbl
   REALTYPE                     :: bRadSrf,bRadSbl
   REALTYPE                     :: Gm1
   REALTYPE                     :: Gt1
   REALTYPE                     :: Gs1
   REALTYPE                     :: dGm1dS
   REALTYPE                     :: dGt1dS
   REALTYPE                     :: dGs1dS
   REALTYPE                     :: f1
   REALTYPE                     :: sl_dpth,sl_z
   REALTYPE                     :: swdk
   REALTYPE                     :: wm
   REALTYPE                     :: ws
   REALTYPE                     :: zgrid,depth

   REALTYPE                     :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk, Rref
   REALTYPE                     :: Uk, Uref, Ustarb, Vk, Vref
   REALTYPE                     :: dR,dU,dV
   REALTYPE                     :: a1, a2, a3
   REALTYPE                     :: cff,cff_up,cff_dn
   REALTYPE                     :: c1,c2,c3
   REALTYPE                     :: dK_bl, hekman, hmonob, sigma, zbl

   REALTYPE, dimension (0:nlev) :: Bflux
   REALTYPE, dimension (0:nlev) :: FC

!  Langmuir enhancement factor
!  Qing Li, 20171213
   REALTYPE                     :: efactor
!  Thickness of surface layer
!  Qing Li, 20171213
   REALTYPE                     :: surfthick
!  10-meter wind
   REALTYPE                     :: wind10m

!-----------------------------------------------------------------------
!BOC
!
!
!-----------------------------------------------------------------------
!  Get approximation of surface layer depth using "epsilon" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
   sl_dpth = epsilon*(z_w(nlev)-zsbl)
   sl_z    = epsilon*zsbl
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux at W-points.
!  Bo is negative, if heat is lost or salinity is gained (convection).
!  It does not include the short wave radiative flux at the surface.
!-----------------------------------------------------------------------
!
    tRadSrf   =   tRad(nlev)
    bRadSrf   =   bRad(nlev)

!  surface buoyancy flux (negative for buoyancy loss)
   Bo         = btFlux + bsFlux

!  include effect of short wave radiation
!  prepare non-local fluxes
   do k = 0,nlev
      Bflux(k)  = Bo  + ( bRadSrf - bRad(k) )
# ifdef NONLOCAL
      cff       = _ONE_-(0.5+sign(0.5d0,Bflux(k)))
      gamh(k)   =  -cff*( tFlux + tRadSrf - tRad(k) )
#  ifdef KPP_SALINITY
      gams(k)   =  -cff*sFlux
#  endif
# endif

   enddo

!-----------------------------------------------------------------------
!  Get Langmuir enhancement factor
!-----------------------------------------------------------------------
! Qing Li, 20171213
   ! 10-meter wind
   wind10m = sqrt(u10**2+v10**2)
   if (langmuir_method .eq. KPP_LT_EFACTOR_MODEL) then
      efactor = kpp_efactor_model(wind10m, u_taus, z_w(nlev)-zsbl)
   else if (langmuir_method .eq. KPP_LT_EFACTOR_READ) then
       ! TODO: read from file <13-12-17, Qing Li> !
      efactor = _ONE_
   else if (langmuir_method .eq. KPP_LT_EFACTOR_SPEC) then
      efactor = kpp_efactor_spec(wav_freq, wav_spec, wav_xcmp, wav_ycmp, &
          u10, v10, u_taus, z_w(nlev)-zsbl)
   else
      ! no enhancement factor applied to ws in bulk Richardson number
      ! when Langmuir enhanced entrainment is on
      efactor = _ONE_
   end if

!-----------------------------------------------------------------------
!  Compute potential density and velocity components surface reference
!  values.
!-----------------------------------------------------------------------
!  update the reference potential density and velocity below if the
!  surface layer is thicker than the uppermost grid and if the tag
!  KPP_AVGSLAYER_REF is defined in cppdefs.h
!  Qing Li, 20171213
!  simply use uppermost value
   Rref = rho(nlev)
   Uref =   u(nlev)
   Vref =   v(nlev)

!-----------------------------------------------------------------------
!  Compute critical function, FC, for bulk Richardson number
!-----------------------------------------------------------------------

   FC(0   ) = _ZERO_
   FC(nlev) = _ZERO_

   do k=nlev-1,2,-1

      depth=z_w(nlev)-z_w(k)

      if (Bflux(k).lt._ZERO_) then
         sigma=min(sl_dpth,depth)
      else
         sigma=depth
      endif

      call wscale (Bflux(k),u_taus,sigma,efactor,wm,ws)

#ifdef KPP_AVGSLAYER_REF
!-----------------------------------------------------------------------
!  Update potential density and velocity components surface reference
!  values.
!-----------------------------------------------------------------------
! Qing Li, 20171213
      ! determine which layer contains surface layer
      surfthick = epsilon*depth
      do kk = nlev,k,-1
         if (z_w(nlev)-z_w(kk-1) .ge. surfthick) then
            kref = kk
            exit
         end if
      end do
      ! Update Rref, Uref and Vref
      if (kref < nlev) then
         Rref = rho(kref)*(surfthick+z_w(kref))
         Uref =   u(kref)*(surfthick+z_w(kref))
         Vref =   v(kref)*(surfthick+z_w(kref))
         do kk = nlev,kref+1,-1
            Rref = Rref + rho(kk)*h(kk)
            Uref = Uref +   u(kk)*h(kk)
            Vref = Vref +   v(kk)*h(kk)
         end do
         Rref = Rref/surfthick
         Uref = Uref/surfthick
         Vref = Vref/surfthick
      end if
#endif

#ifdef KPP_TWOPOINT_REF

!     interpolate reference value at grid face "k"
!     from values at grid centers

      cff = _ONE_/h_r(k)

      dR  = cff*( rho(k+1) - rho(k) )
      dU  = cff*( u  (k+1) - u  (k) )
      dV  = cff*( v  (k+1) - v  (k) )

      cff = _ONE_/2.0

      Rk  = rho(k) + h(k)*cff*dR
      Uk  =   u(k) + h(k)*cff*dU
      Vk  =   v(k) + h(k)*cff*dV

#else
!     identify reference value at grid face "k"
!     with value at center above

      Rk = rho(k+1)
      Uk =   u(k+1)
      Vk =   v(k+1)

!!$      c1 =  0.6
!!$      c2 =  0.2
!!$      c3 =  0.2
!!$
!!$      Rk = c1*rho(k+1) + c2*rho(k) + c3*rho(k-1)
!!$      Uk = c1*  u(k+1) + c2*  u(k) + c3*  u(k-1)
!!$      Vk = c1*  v(k+1) + c2*  v(k) + c3*  v(k-1)
#endif

!     compute numerator and denominator of Ri_b
      Ritop = -gorho0*(Rref-Rk)*depth

      Ribot = (Uref-Uk)**2+(Vref-Vk)**2 +                               &
              Vtc*depth*ws*sqrt(abs(NN(k)))

# ifdef KPP_IP_FC
      FC(k) = Ritop-Ric*Ribot
# else
      FC(k) = Ritop/(Ribot+eps)
# endif

   enddo  ! inner grid faces



!-----------------------------------------------------------------------
! Linearly interpolate to find "zsbl" where Rib/Ric=1.
!-----------------------------------------------------------------------

   ksbl = 1         ! ksbl is index of cell containing zsbl
   zsbl = z_w(0)

# ifdef KPP_IP_FC
!  look for position of vanishing F_crit
   do k=nlev,2,-1
      if ((ksbl.eq.1).and.(FC(k-1).gt._ZERO_)) then
         zsbl = (z_w(k)*FC(k-1)-z_w(k-1)*FC(k))/(FC(k-1)-FC(k))
         ksbl = k
      endif
   enddo
# else
!  look for position of vanishing Ri_b
   do k=nlev,2,-1
      if ((ksbl.eq.1).and.((FC(k).lt.Ric).and.(FC(k-1).ge.Ric))) then
         zsbl = ((FC(k-1)-Ric)*z_w(k) -                                 &
                         (FC(k)-Ric)*z_w(k-1))/(FC(k-1)-FC(k))
         ksbl = k
      endif
   enddo
# endif



!-----------------------------------------------------------------------
!  Compute total buoyancy flux at surface boundary layer depth
!-----------------------------------------------------------------------

!  interpolate from interface values to zsbl
   bRadSbl = ( bRad(ksbl-1)*(z_w(ksbl)-zsbl) +                                &
               bRad(ksbl  )*(zsbl-z_w(ksbl-1) ) )/ h(ksbl)

   Bfsfc   = Bo + (bRadSrf - bRadsbl)


!-----------------------------------------------------------------------
!  Limit boundary layer depth by Ekman and Monin-Obukhov depths
!  (under neutral and stable conditions)
!-----------------------------------------------------------------------

   if (clip_mld) then
      if ((u_taus.gt._ZERO_).and.(Bfsfc.gt._ZERO_)) then
         hekman = cekman*u_taus/max(abs(f),eps)
         hmonob = cmonob*u_taus*u_taus*u_taus/max(kappa*Bfsfc,eps)
         zsbl   = (z_w(nlev)-min(hekman,hmonob,z_w(nlev)-zsbl))
      endif
   endif


   zsbl = min(zsbl,z_w(nlev))
   zsbl = max(zsbl,z_w(0   ))

!  find new boundary layer index "ksbl".
   ksbl=1
   do k=nlev,2,-1
      if ((ksbl.eq.1).and.(z_w(k-1).lt.zsbl)) then
         ksbl = k
      endif
   enddo


!-----------------------------------------------------------------------
!  Compute total buoyancy flux at surface boundary layer depth
!-----------------------------------------------------------------------


   bRadSbl = ( bRad(ksbl-1)*(z_w(ksbl)-zsbl) +                          &
               bRad(ksbl  )*(zsbl-z_w(ksbl-1) ) )/ h(ksbl)


   Bfsfc   = Bo + (bRadSrf - bRadSbl)

!-----------------------------------------------------------------------
!  Update Langmuir enhancement factor
!-----------------------------------------------------------------------
! Qing Li, 20171213
   if (langmuir_method .eq. KPP_LT_EFACTOR_MODEL) then
      efactor = kpp_efactor_model(wind10m, u_taus, z_w(nlev)-zsbl)
   else if (langmuir_method .eq. KPP_LT_EFACTOR_READ) then
       ! TODO: read from file <13-12-17, Qing Li> !
      efactor = _ONE_
   else if (langmuir_method .eq. KPP_LT_EFACTOR_SPEC) then
      efactor = kpp_efactor_spec(wav_freq, wav_spec, wav_xcmp, wav_ycmp, &
          u10, v10, u_taus, z_w(nlev)-zsbl)
   else if (langmuir_method .eq. KPP_LT_ENTRAINMENT) then
       ! TODO: efactor for entrainment <13-12-17, Qing Li> !
      efactor = _ONE_
   else
      efactor = _ONE_
   end if
   ! DEBUG QL
   ! LEVEL2 'efactor = ', efactor

!-----------------------------------------------------------------------
!  Compute tubulent velocity scales (wm,ws) at "zsbl".
!-----------------------------------------------------------------------

   zbl     = z_w(nlev)-zsbl
   sl_dpth = epsilon*zbl

   if (Bfsfc.gt._ZERO_) then
      cff  = _ONE_
   else
      cff  = epsilon
   endif

   sigma=cff*(z_w(nlev)-zsbl)

   call wscale (Bfsfc,u_taus,sigma,efactor,wm,ws)


!-----------------------------------------------------------------------
!  Compute nondimensional shape function Gx(sigma) in terms of the
!  interior diffusivities at sigma=1 (Gm1, Gs1, Gt1) and its vertical
!  derivative evaluated "zsbl" via interpolation.
!-----------------------------------------------------------------------

!   original code with kappa-bug
!   f1 = 5.0*max(_ZERO_,Bfsfc)*kappa/(u_taus*u_taus*u_taus*u_taus+eps)

!  new code code without kappa-bug
   f1 = 5.0*max(_ZERO_,Bfsfc)/(u_taus*u_taus*u_taus*u_taus+eps)

   if (zsbl.gt.z_w(1)) then
      ! if boundary layer does not touch lowest grid box

!     abbreviation for "ksbl"
      k      = ksbl
      cff    = _ONE_/h(k)
      cff_dn = cff*(zsbl-z_w(k-1))
      cff_up = cff*(z_w(k)-zsbl  )

!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "zsbl".

      K_bl   = cff_dn*num(k)+cff_up*num(k-1)
      dK_bl  = cff*(num(k)-num(k-1))
      Gm1    = K_bl/(zbl*wm+eps)

# ifdef KPP_CLIP_GS
      dGm1dS = min(_ZERO_,-dK_bl/(wm+eps)-K_bl*f1)
# else
      dGm1dS = -dK_bl/(wm+eps) + K_bl*f1
# endif

!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "zsbl".
!
      K_bl   = cff_dn*nuh(k)+cff_up*nuh(k-1)
      dK_bl  = cff*(nuh(k)-nuh(k-1))
      Gt1    = K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
      dGt1dS = min(_ZERO_,-dK_bl/(ws+eps)-K_bl*f1)
# else
      dGt1dS = -dK_bl/(ws+eps) + K_bl*f1
# endif

# ifdef KPP_SALINITY
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "zsbl".
!
      K_bl   = cff_dn*nus(k)+cff_up*nus(k-1)
      dK_bl  = cff*(nus(k)-nus(k-1))
      Gs1    = K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
      dGs1dS = min(_ZERO_,-dK_bl/(ws+eps)-K_bl*f1)
# else
      dGs1dS = -dK_bl/(ws+eps) + K_bl*f1
# endif

# endif


   else
!     If the surface boundary layer extends to the bottom, assume that
!     the neutral boundary layer similarity theory holds at the bottom.
!
      dK_bl  = kappa*u_taub
      K_bl   = dK_bl*(zsbl-z_w(0))

!     Compute nondimensional bottom shape function for diffusion of
!     momentum

      Gm1    = K_bl/(zbl*wm+eps)

# ifdef KPP_CLIP_GS
      dGm1dS = min(_ZERO_,-dK_bl/(wm+eps)-K_bl*f1)
# else
      dGm1dS = -dK_bl/(wm+eps) + K_bl*f1
# endif


!     Compute nondimensional bottom shape function for diffusion of
!     temperature.

      Gt1    = K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
      dGs1dS = min(_ZERO_,-dK_bl/(ws+eps)-K_bl*f1)
# else
      dGs1dS = -dK_bl/(ws+eps) + K_bl*f1
# endif

!     Compute nondimensional bottom shape function for diffusion of
!     salinity.
# ifdef KPP_SALINITY
      Gs1    = Gt1
      dGs1dS = dGt1dS
# endif

   endif
!
!-----------------------------------------------------------------------
!  Compute surface boundary layer mixing coefficients
!-----------------------------------------------------------------------
!
!  loop over the inner interfaces
!  (the outermost are not needed)
   do k=1,nlev-1

      if (k.ge.ksbl) then    ! interface above cell containing zsbl

!        Compute turbulent velocity scales at vertical W-points.
         depth=z_w(nlev)-z_w(k)
         if (Bflux(k).lt._ZERO_) then
            sigma=min(sl_dpth,depth)
         else
            sigma=depth
         endif

         call wscale(Bflux(k),u_taus,sigma,efactor,wm, ws)
!
!        Set polynomial coefficients for shape function.
         sigma = depth/(zbl+eps)

         a1    = sigma-2.0
         a2    = 3.0-2.0*sigma
         a3    = sigma-1.0
!
!        Compute nondimesional shape functions.
         Gm = a1+a2*Gm1+a3*dGm1dS
         Gt = a1+a2*Gt1+a3*dGt1dS
# ifdef KPP_SALINITY
         Gs = a1+a2*Gs1+a3*dGs1dS
# endif

!        Compute boundary layer mixing coefficients, combine them
!        with interior mixing coefficients.
         num(k) = depth*wm*(_ONE_+sigma*Gm)
         nuh(k) = depth*ws*(_ONE_+sigma*Gt)
# ifdef KPP_SALINITY
         nus(k) = depth*ws*(_ONE_+sigma*Gs)
# endif
# ifdef NONLOCAL
!        Compute boundary layer nonlocal transport (m/s2).
         cff      = Cg*(_ONE_-(0.5+sign(0.5d0,Bflux(k))))/(zbl*ws+eps)
         gamh(k) = cff*nuh(k)*gamh(k)
#  ifdef KPP_SALINITY
         gams(k) = cff*nus(k)*gams(k)
#  endif
# endif


!     Set non-local transport terms to zero outside boundary  layer.
      else
# ifdef NONLOCAL
         gamh(k) = _ZERO_
#  ifdef KPP_SALINITY
         gams(k) = _ZERO_
#  endif
# endif
      endif


   enddo

!     not non-local fluxes at top and bottom
      gamh(0   ) = _ZERO_
      gams(0   ) = _ZERO_
      gamh(nlev) = _ZERO_
      gams(nlev) = _ZERO_


 end subroutine surface_layer
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

