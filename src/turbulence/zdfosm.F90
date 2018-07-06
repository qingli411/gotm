!$Id: zdfosm.F90,v 1.34 2007-03-15 10:52:07 kbk Exp $
#include"cppdefs.h"

! BOP
!
! ! MODULE: osmosis: the OSMOSIS-turbulence model
!
! ! INTERFACE
MODULE zdfosm
!
! !DESCRIPTION
!
!
! Version 1.1
!
! Model for Langmuir turbulence in the surface boundary layer. The model is a KPP model in that it
! parametrizes the eddy diffusion coefficients in terms of velocity and length scales. It includes
! non-gradient terms which include the contribution to the transports by the Stokes shear.
!
! The full model uses interior mixing and bottom boundary routines from the Large et al (1994) KPP model.
! This is a deep water version in that it does not attempt to represent overlaps between the surface and
! bottom boundary layers.
!
! HISTORY
!
! 13/11/17 - New version, Vn 1.1
! 14/11/17 - Change to stable boundary layer depth equation. Removes need for initial boundary
!            layer depth.
! 15/11/17 - Remove limit on Langmuir number (needed because flux-gradient relationships are
!             for Langmuir turbulence, and can fail as Langmuir number gets large).
!              Alter uw and vw scales in Stokes component of flux-gradient relationship. Modify
!              so (1-6.5 La**(8/3)) limited to minimum of 0.2.
!              Remove La**(8/3) in other scales.
!              Limit maximum value of La to 4. Large enough, avoids problems with very large
!              values.
!
   use meanflow,         only: u,v,w,t,s             ! mean flow variables.
   use meanflow,         only: buoy                  ! buoyancy
   use meanflow,         only: ff => cori            ! Coriolis parameter
   use meanflow,         only: NN,NNT,NNS,SS         ! stratification and shear.
   use meanflow,         only: cp                    ! specific heat of sea water
   use turbulence,       only: num,nuh,nus,nucl      ! eddy viscosity and diffusivities
   use turbulence,       only: gamu,gamv,gamh,gams   ! non-gradient terms in flux-gradient
! initialize tke in module turbulence
! Qing Li, 20180404
   use turbulence,       only: tke

   ! use turbulence,       only: gamh_f,gams_f
   use turbulence,       only: Rig
   use turbulence,       only: kappa
   use observations,     only: rn_abs => a
   use observations,     only: rn_si0 => g1
   use observations,     only: rn_si1 => g2          ! double exponential radiation profile.
! surface Stokes drift and penetration depth are now defined in observation.F90
! Qing Li, 20180405
   use observations,     only: us_x,us_y             ! components of Stokes drift (ms-1)
   use observations,     only: delta                 ! Stokes Penetration Depth (m)
                                                     ! surface fluxes
   use airsea,           only: tx,ty                 ! kinematic stress components (m2s-2)
   use airsea,           only: I_0                   ! Shortwave radiations
   use airsea,           only: heat                  ! heat flux
                                                     ! (Wm-2)
   ! use airsea,           only: p_e                   ! precip-evap (ms-1)
   ! Qing Li, 20180403
   use airsea,           only: precip, evap

   IMPLICIT NONE

   private                                          ! sets default to private

! !PUBLIC MEMBER FUNCTIONS:
!
   public init_osm, do_osm, clean_osm

! !PUBLIC DATA MEMBERS:                            ! these are retained for original KPP routines.
!
! prognostic depth of the boundary layer
   REALTYPE, public                 :: hbl
! depth of mixed layer
   REALTYPE, public                 :: hml
! thickness of stable pycnocline
   REALTYPE, public                 :: dh
! depth of mixed layer on last unstable timestep
   REALTYPE, public                 :: hbli
! depth of boundary layer at change to stable conditions
   REALTYPE, public                 :: h_i
! z-position of surface boundary layer depth - grid value
   REALTYPE, public                 :: zsbl
   REALTYPE, public                 :: zsml
   REALTYPE, PUBLIC                 :: wth_0
   REALTYPE, PUBLIC                 :: ws_0
   REALTYPE, PUBLIC                 :: wb_0
   REALTYPE, public                 :: uw_0,vw_0
   REALTYPE, PUBLIC                 :: ustar,wstrl

   REALTYPE, public, dimension(:), allocatable  :: dtdz_pyc,dsdz_pyc

   LOGICAL, public                  ::  lconv ! unstable/stable bl

! !REVISION HISTORY:
!  Original author(s): Alan Grant
!
!  $Log: osmosis.F90,v $
!
! !LOCAL VARIABLES:
!
!
   INTEGER :: kk

!  reference density
   REALTYPE                              ::    rho_0

!  g/rho_0
   REALTYPE                              ::    gorho0

!  compute surface and bottom BBL
   logical                               ::    osm_sbl,osm_bbl

!  compute internal mixing
   logical                               ::    osm_interior

   REAL                                  ::    emp,rnf,qns,qsr,sfx,utau,vtau   ! flux variables for NEMO code
!  positions of grid faces and centers
   REALTYPE, dimension(:), allocatable   ::    z_w,z_r

!  distance between centers
   REALTYPE, dimension(:), allocatable   ::    h,h_r

! NEMO code arrays

   REALTYPE, dimension(:), allocatable   :: fs_depw,fs_dept
   REALTYPE, dimension(:), allocatable   :: e3w_n,e3t_n
   REALTYPE, dimension(:), allocatable   :: ub,vb,wb
   REALTYPE, dimension(:,:), allocatable :: tsn
   REALTYPE, dimension(:), allocatable   :: rhd  ! (rho-rau0)/rau0=0
   REALTYPE, dimension(:,:), allocatable :: rab_n ! thermal/haline coeffs
   REALTYPE, dimension(:), allocatable :: rhop   ! density
   REALTYPE, dimension(:), allocatable :: zghamt,zghams,zgamu,zgamv
   REALTYPE, dimension(:), allocatable :: zavt,zavmu,zavmv
!
! replaces USE oce in NEMO version, should be in SR zdf_osm
!
   REALTYPE, dimension(:), allocatable   :: zviscos   ! temp. array for viscosities use ua as workspace
   REALTYPE, dimension(:), allocatable   :: zdiffut ! temp. array for diffusivities use sa as workspace

   REALTYPE, DIMENSION(:), allocatable :: zdtdz_pyc ! parametrized gradient of temperature in pycnocline
   REALTYPE, DIMENSION(:), allocatable :: zdsdz_pyc ! parametrised gradient of salinity in pycnocline
   REALTYPE, DIMENSION(:), allocatable :: zdbdz_pyc ! parametrised gradient of buoyancy in the pycnocline
   REALTYPE, DIMENSION(:), allocatable :: zdudz_pyc ! u-shear across the pycnocline
   REALTYPE, DIMENSION(:), allocatable :: zdvdz_pyc ! v-shear across the pycnocline
   REALTYPE, DIMENSIOn(:), allocatable :: zavtb, zavmb ! minimum diffusivity and viscosity
   ! REALTYPE, DIMENSION(:), allocatable :: zghamt_f,zghams_f
   ! parameters for NEMO code

   INTEGER, parameter :: wp = KIND( _ONE_ )  ! single precision - defined in cppdefs.h

   REAL, parameter :: grav=9.81
   REAL, parameter :: rau0=1026.0,r1_rau0=1.0/rau0
   REAL, parameter :: rcp=4000.0,r1_rcp=1.0/rcp
   REAL, parameter :: rau0_rcp=rau0*rcp,r1_rau0_rcp=1.0/rau0_rcp
   REAL, parameter :: rcs=1.0
   REAL, parameter :: epsln=1.0e-20
!
! Parameters use by interior mixing roufine
!
!!  critical gradient Richardson number below which turbulent
!  mixing occurs (Ri0=0.7)
   REALTYPE, parameter :: Ri0     = 0.7

!  buoancy frequency (1/s2) limit for convection (bvfcon=-2.0E-5)
   REALTYPE, parameter :: bvfcon  = -2.0E-5

!  value of double-diffusive density ratio where mixing goes
!  to zero in salt fingering (Rrho0=1.9)
   REALTYPE, parameter :: Rrho0    = 1.9
!  scaling factor for double diffusion of temperature in salt
!  fingering case (fdd=0.7)
   REALTYPE, parameter :: fdd     = 0.7

!  maximum interior convective viscosity and diffusivity
!  due to shear instability (nu0c=0.01)
   REALTYPE, parameter :: nu0c    = 0.01

!  maximum interior viscosity (m2/s) due to shear
!  instability (nu0m=10.0E-4)
   REALTYPE, parameter :: nu0m    = 10.0E-4

!  maximum interior diffusivity (m2/s) due to shear
!  instability (nu0s=10.0E-4)
   REALTYPE, parameter :: nu0s    = 10.0E-4

!  scaling factor for double diffusion in salt fingering (nu=1.5E-6)
   REALTYPE, parameter :: nu      = 1.5E-6

!  scaling factor for double diffusion in salt fingering (nuf=10.0E-4)
   REALTYPE, parameter :: nuf     = 10.0E-4

!  interior viscosity (m2/s) due to wave breaking (nuwm=1.0E-5)
   REALTYPE, parameter :: nuwm    = 1.0E-5

!  interior diffusivity (m2/s) due to wave breaking (nuwm=1.0E-6)
   REALTYPE, parameter :: nuws    = 1.0E-6

!  double diffusion constant for salinity in diffusive
!  convection case (sdd1=0.15)
   REALTYPE, parameter :: sdd1    = 0.15

!  double diffusion constant for salinity in diffusive convection
!  (sdd2=1.85)
   REALTYPE, parameter :: sdd2    = 1.85

!  double diffusion constant for salinity in diffusive convection
!  (sdd3=0.85)
   REALTYPE, parameter :: sdd3    = 0.85

!  double diffusion constant for temperature in diffusive convection
!  (tdd1=0.909)
   REALTYPE, parameter :: tdd1    = 0.909

!  double diffusion constant for temperature in diffusive convection
!  (tdd2=4.6)
   REALTYPE, parameter :: tdd2    = 4.6

!  double diffusion constant for temperature in diffusive convection case
!  (tdd3=0.54).
   REALTYPE, parameter :: tdd3    = 0.54

!EOP
!-------------------------------------------------------------------------

   contains

  !
! !IROUTINE: Loop over the KPP-algorithm
!
! !INTERFACE:
   subroutine do_osm(nlev,h0,h,dt)
!
! !DESCRIPTION:
! Here, the time step for the KPP model is managed. If {\tt kpp\_interior=.true.}
! in {\tt kpp.nml}, the mixing algorithm for the computation of the interior
! diffusivities is called first. This algorithm is described in \sect{sec:kppInterior}.
! Then, if {\tt kpp\_sbl=.true.} and {\tt kpp\_bbl=.true.}, the algorithms
! for the surface and bottom boundary layer are called. They are described in
! \sect{sec:kppSurface} and \sect{sec:kppBottom}, respectively.
!
! If this routine is called from a three-dimensional code, it is essential to
! pass the correct arguments. The first 3 parameters relate to the numerical
! grid, discussed in \sect{SectionNumericsMean}. Note that {\tt h0} denotes
! the local bathymetry, i.e.\ the positive distance between the reference level
! $z=0$ and the bottom.

! The next three parameters denote the potential density, $\rho$,
! and the two mean velocity components, $U$ and $V$. The buoyancy frequency, $N**2$,
! and the different contributions to it, $N_\Theta**2$ and $N_S**2$, have to be computed
! from the potential density as discussed in \sect{sec:stratification}. The shear frequency,
! $M**2$, is defined in \eq{MSquared}. The vertical discretisation does not necessarly
! have to follow \eq{shearsquared}, since in the KPP model no TKE equation is solved
! and thus energy conservation is not an issue. All three-dimensional fields have to
! be interpolated "in a smart way" to the water column defined by GOTM. The corresponding
! interpolation schemes may be quite different for the different staggered grids,
! finite volume, and finite element approaches used in the horizontal. Therefore,
! we cannot offer a general recipe here.

! The bottom friction velocity is computed as described in \sect{sec:friction}. If this
! parameter is passed from a three-dimensional code, it has to be insured that the parameter
! $r$ in \eq{uStar} is computed consistently, see \eq{rParam}.
!
! All fluxes without exception are counted positive, if they enter the water body.
! Note that for consistency, the equations of state in GOTM cannot
! be used if the KPP routines are called from a 3-D model. Therefore,
! it is necessary to pass the temperature and salinity fluxes, as well as the
! corresponding buoyancy fluxes. The same applies
! to the radiative fluxes. The user is responsible for
! performing the flux conversions in the correct way. To get an idea
! have a look at \sect{sec:convertFluxes}.
!
! The last argument is the Coriolis parameter, $f$. It is only used for clippling the mixing
! depth at the Ekman depth.
!
!
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

! model timestep
   REALTYPE                                      :: dt

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE                            :: cff
   integer                             :: k

!
!-----------------------------------------------------------------------
!BOC
!
!-----------------------------------------------------------------------
! Update model grid
!-----------------------------------------------------------------------

!  Compute distance between centers (between rho-points)
!  Note that h is the distance between faces (between w-points)
   do k=1,nlev
      h_r(k) = 0.5*(h(k)+ h(k+1))
   enddo

!  Compute position of interfaces (w-points)
   z_w(0) =  -h0
   do k=1,nlev
      z_w(k) = z_w(k-1) + h(k)
   enddo

!  Compute position of centers (rho-points)
   z_r(1) = - h0 + 0.5*h(1)
   do k=2,nlev
      z_r(k) = z_r(k-1) + h_r(k-1)
   enddo

!-----------------------------------------------------------------------
! compute interior mixing
!-----------------------------------------------------------------------

   if (osm_interior) then
      call zdfosm_interior(nlev)
   else
      num = 0._wp
      nuh = 0._wp
      nus = 0._wp
   endif


!-----------------------------------------------------------------------
! compute surface boundary layer mixing
!-----------------------------------------------------------------------

   if (osm_sbl) then

! NEMO code interface depth array, name, reverse order and sign.

      kk=1
      do k=nlev,1,-1
         fs_depw(kk) = -z_w(k)
         fs_dept(kk) = -z_r(k)
         e3w_n(kk) = h(k)
         e3t_n(kk) = h_r(k)
         kk=kk+1
      enddo

! reverse indices and names of model variables, put temperature/salinity into 2D array.

      kk=1
      do k=nlev,1,-1
         ub(kk) = u(k)
         vb(kk) = v(k)
         wb(kk) = w(k)
         tsn(kk,1) = t(k)
         tsn(kk,2) = s(k)
         kk=kk+1
      enddo

! properties of sea water

      do k=1,nlev
         rhop(k) = rau0
         rhd(k) = 0._wp
         rab_n(k,1) = 0.230    ! value for comparison with Georges NEMO simulation
         rab_n(k,2) = 0.76554
      enddo

      kk=1
      do k=nlev,1,-1
         zavt(kk) = nuh(k)
         zavmu(kk) = num(k)
         zavmv(kk) = num(k)
         kk=kk+1
      enddo
      do k=1,nlev
         zghamt(k)=0._wp
         zghams(k)=0._wp
         zgamu(k)=0._wp
         zgamv(k)=0._wp
         ! zghamt_f(k)=0._wp
         ! zghams_f(k)=0._wp
      enddo

      ! emp = -p_e
      ! Qing Li, 20180403
      emp = -(precip+evap)
      rnf=0.0
      qns = heat
      qsr = I_0
      utau = tx * rau0   ! undoes division by rau0 in gotm
      vtau = ty * rau0
      sfx=0.0
      call zdf_osm( dt, nlev )
!
! put NEMO flux-gradient arrays in GOTM arrays.
!
      kk=nlev
      do k=1,nlev
         nuh(kk) = zavt(k)
         nus(kk) = zavt(k)
         num(kk) = zavmu(k)
!       num(kk) = zavmv(k) ! shouldn't matter which viscosity is used - test NEMO by trying out both, results should be the same
         gamh(kk) = zghamt(k)
         gams(kk) = zghams(k)
         gamu(kk) = zgamu(k)
         gamv(kk) = zgamv(k)
         dtdz_pyc(kk) = zdtdz_pyc(k)
         dsdz_pyc(kk) = zdbdz_pyc(k)
         ! gamh_f(kk) = zghamt_f(k)
         ! gams_f(kk) = zghams_f(k)
         kk=kk-1
      enddo
      nuh(0)=0._wp
      nus(0)=0._wp
      num(0)=0._wp
      gamh(0)=0._wp
      gams(0)=0._wp
      gamu(0)=0._wp
      gamv(0)=0._wp
   endif
!-----------------------------------------------------------------------
! compute bottom boundary layer mixing
!-----------------------------------------------------------------------

!   if (osm_bbl) then
!      call bottom_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,             &
!                        0._wp,0._wp,0._wp,0._wp,tRad,bRad,f)
!   endif

 end subroutine do_osm
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the KPP module
!
! !INTERFACE:
   subroutine init_osm(namlst,fn,nlev,h0,h)
!
! !DESCRIPTION:
! This routine first reads the namelist {\tt kpp}, which has to be contained
! in a file with filename specified by the string {\tt fn} (typically called
! {\tt kpp.nml}). Since the {\tt kpp} module uses fields defined in the
! {\tt turbulence} module, it has to allocate dynamic memory for them.
! Apart from this, this routine reports the model settings and initialises a
! number of parameters needed later in the time loop.
!
! If you call the GOTM KPP routines from a three-dimensional model, make sure
! that this function is called \emph{after} the call to {\tt init\_turbulence()}.
! Also make sure that you pass the correct arguments.
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

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: k
   integer                             :: rc

   namelist /osm/                      osm_sbl,osm_bbl,osm_interior,h_i
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_osm...'

   ! read the variables from the namelist file
write(*,*) fn
   open(namlst,file=fn,status='old',action='read',err=80)

   LEVEL2 'reading osm namelist...'

   read(namlst,nml=osm,err=81)
   close (namlst)

   LEVEL2 'done.'

   LEVEL2 'initial value for hbl ',h_i
   hbl=h_i
   hbli=h_i
   dh=0.2*hbl
   lconv = .FALSE.

   allocate(num(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (num)'
   num = 1.0D-6

   allocate(nuh(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (nuh)'
   nuh = 1.0D-6

   allocate(nus(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (nus)'
   nus = 1.0D-6

   allocate(nucl(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (nucl)'
   nucl = _ZERO_

   allocate(gamu(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamu)'
   gamu = 0._wp

   allocate(gamv(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamv)'
   gamv = 0._wp

   allocate(gamh(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamh)'
   gamh = 0._wp

   ! allocate(gamh_f(0:nlev),stat=rc)
   ! if (rc /= 0) stop 'init_turbulence: Error allocating (gamh)'
   ! gamh_f = 0._wp

   allocate(gams(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gams)'
   gams = 0._wp

   ! allocate(gams_f(0:nlev),stat=rc)
   ! if (rc /= 0) stop 'init_turbulence: Error allocating (gams_f)'
   ! gams_f = 0._wp

   allocate(Rig(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (Rig)'
   Rig = _ZERO_

!  allocate memory for tke
!  Qing Li, 20180404

   LEVEL2 'allocation memory..'
   allocate(tke(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (tke)'
   tke = _ZERO_

   allocate(dtdz_pyc(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (dtdz_pyc)'
   dtdz_pyc = 0._wp

   allocate(dsdz_pyc(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (dsdz_pyc)'
   dsdz_pyc = 0._wp

! local variables for NEMO flux-gradient relationship

   allocate(fs_depw(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (fs_depw) '
   fs_depw=0._wp

   allocate(fs_dept(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (fs_dept) '
   fs_dept=0._wp

   allocate(e3w_n(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (e3w_n) '
   e3w_n=0._wp

   allocate(ub(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (ub) '
   ub=0._wp

   allocate(vb(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (vb) '
   vb=0._wp

   allocate(wb(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (vb) '
   wb=0._wp

   allocate(tsn(1:nlev,2),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (tsn) '
   tsn=0._wp

   allocate(rhd(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (rhd) '
   rhd=0._wp

   allocate(rhop(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (rhop) '
   rhop=0._wp

   allocate(rab_n(1:nlev,2),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (rab_n) '
   rab_n=0._wp

   allocate(e3t_n(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (e3t_n) '
   e3t_n=0._wp

   allocate(zavt(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zavt)'
   zavt = 0._wp

   allocate(zviscos(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zviscos)'
   zviscos = 0._wp

   allocate(zdiffut(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zdiffut)'
   zdiffut = 0._wp

   allocate(zavmu(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zavmu)'
   zavmu = 0._wp

   allocate(zavmv(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zavmv)'
   zavmv = 0._wp

   allocate(zgamu(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zgamu)'
   zgamu = 0._wp

   allocate(zgamv(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zgamv)'
   zgamv = 0._wp

   allocate(zghamt(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zghamt)'
   zghamt = 0._wp

   ! allocate(zghamt_f(1:nlev),stat=rc)
   ! if (rc /= 0) stop 'init_turbulence: Error allocating (zghamt_f)'
   ! zghamt_f = 0._wp

   allocate(zghams(1:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (zghams)'
   zghams = 0._wp

   ! allocate(zghams_f(1:nlev),stat=rc)
   ! if (rc /= 0) stop 'init_turbulence: Error allocating (zghams_f)'
   ! zghams_f = 0._wp

   allocate(zdtdz_pyc(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zdtdz_pyc)'
   zdtdz_pyc=0._wp

   allocate(zdbdz_pyc(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zdbdz_pyc)'
   zdbdz_pyc=0._wp

   allocate(zdsdz_pyc(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zdsdz_pyc)'
   zdsdz_pyc=0._wp

   allocate(zdudz_pyc(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zdudz_pyc)'
   zdudz_pyc=0._wp

   allocate(zdvdz_pyc(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zdvdz_pyc)'
   zdvdz_pyc=0._wp

   allocate(zavtb(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zavtb)'
   zavtb = 1.e-15

   allocate(zavmb(1:nlev),stat=rc)
   if (rc /= 0 ) stop 'init_turbulence: Error allocating (zavmb)'
   zavmb = 1.e-15

!
!  allocate memory for variables defined in other modules
!

   allocate(z_w(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (z_w)'
   z_w = 0._wp

   allocate(z_r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (z_r)'
   z_r = 0._wp

   allocate(h_r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (h_r)'
   h_r = 0._wp

!  report model parameters

   LEVEL2 '--------------------------------------------------------'
   LEVEL3 'You are using the KPP turbulence model          '
   LEVEL3 'with the following specifications:                      '
   LEVEL3 '                                                        '
   if (osm_interior) then
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

   if (osm_sbl) then
      LEVEL3 'Surface layer mixing algorithm             - active -   '
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


   else
      LEVEL3 'Surface layer mixing algorithm         - not active -   '
   endif

   if (osm_bbl) then
      LEVEL3 'Bottom layer mixing algorithm              - active -   '
      LEVEL4 '(Same parameters as surface layer mixing)'

   else
      LEVEL3 'Bottom layer mixing algorithm          - not active -   '
   endif

   LEVEL2 '--------------------------------------------------------'

   LEVEL1 'done.'

   return

80 FATAL 'I could not open "kpp.nml"'
   stop 'init_kpp'
81 FATAL 'I could not read "kpp" namelist'
   stop 'init_kpp'

 end subroutine init_osm
!EOC

!
   SUBROUTINE zdf_osm( rn_rdt, nlev )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_osm  ***
      !!
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!      coefficients and non local mixing using the OSMOSIS scheme
      !!
      !! ** Method :   The boundary layer depth hkpp is diagnosed at tracer points
      !!      from profiles of buoyancy, and shear, and the surface forcing.
      !!      Above hbl (sigma=-z/hbl <1) the mixing coefficients are computed from
      !!
      !!                      Kx =  hkpp  Wx(sigma) G(sigma)
      !!
      !!             and the non local term ghamt = Cs / Ws(sigma) / hkpp
      !!      Below hkpp  the coefficients are the sum of mixing due to internal waves
      !!      shear instability and double diffusion.
      !!
      !!      -1- Compute the now interior vertical mixing coefficients at all depths.
      !!      -2- Diagnose the boundary layer depth.
      !!      -3- Compute the now boundary layer vertical mixing coefficients.
      !!      -4- Compute the now vertical eddy vicosity and diffusivity.
      !!      -5- Smoothing
      !!
      !!        N.B. The computation is done from jk=2 to jpkm1
      !!             Surface value of avt avmu avmv are set once a time to zero
      !!             in routine zdf_kpp_init.
      !!
      !! ** Action  :   update the non-local terms ghamts
      !!                update avt, avmu, avmv (before vertical eddy coef.)
      !!
      !! References : Large W.G., Mc Williams J.C. and Doney S.C.
      !!         Reviews of Geophysics, 32, 4, November 1994
      !!         Comments in the code refer to this paper, particularly
      !!         the equation number. (LMD94, here after)
      !!----------------------------------------------------------------------
      !!
      IMPLICIT NONE
      INTEGER, INTENT( in  ) ::   nlev   ! number of levels
      REAL(wp), INTENT( in ) :: rn_rdt   ! ocean timestep
      !!
      INTEGER ::   ji, jj, jk              ! dummy loop indices
      INTEGER ::   ikbot, jkmax, jpkm1, jkp2   !

      REAL(wp) ::   ztx, zty, zstabl, zbuofdep,zucube     !
      REAL(wp) ::   zrhos, zalbet, zbeta, zthermal, zatt1   !
      REAL(wp) ::   zref, zt, zs, zh, zu, zv, zrh                   ! Bulk richardson number
      REAL(wp) ::   zrib, zrinum, zdVsq, zVtsq                      !
      REAL(wp) ::   zehat, zeta, zhrib, zsig, zscale, zwst, zws, zwm   ! Velocity scales
      REAL(wp) ::   zwsun, zwmun, zcons, zconm, zwcons, zwconm      !
      REAL(wp) ::   zsr, zbw, ze, zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zcomp , zrhd,zrhdr,zbvzed   ! In situ density
      ! REAL(wp) ::   zt,zs,zu,zv,zrh ! variables used in constructing averages
      REAL(wp) ::   zdb,zvel_max
! Scales
      REAL(wp) :: zBo
      REAL(wp) :: zBosol
      REAL(wp) :: zustar ! friction velocity
      REAL(wp) :: zwstrl ! Langmuir velcotiy scale
      REAL(wp) :: zvstr ! Velocity scale tends to ustr for large Langmuir numbers.
      REAL(wp) :: zwstrc  ! Convective velocity scale
      REAL(wp) :: zuw0 ! Surface u-momentum flux
      REAL(wp) :: zvw0 ! Surface v-momentum flux
      REAL(wp) :: zwth0 ! Surface heat flux (Kinematic)
      REAL(wp) :: zws0 ! Surface freshwater flux
      REAL(wp) :: zwb0  ! Surface buoyancy flux
! 26/05/17 radiative buoyancy flux removed
! change to NEMO code ???
      REAL(wp) :: zwthav ! Heat flux - bl average
      REAL(wp) :: zwsav ! freshwater flux - bl average
      REAL(wp) :: zwbav ! Buoyancy flux - bl average
      REAL(wp) :: zustke ! Surface Stokes drift
      REAL(wp) :: zla    ! Trubulent Langmuir number
      REAL(wp) :: zcos_wind ! Cos angle of surface stress
      REAL(wp) :: zsin_wind ! Sin angle of surface stress
      REAL(wp) :: zhol ! Stability parameter for boundary layer

     ! mixed-layer variables

      INTEGER :: ibld ! level of boundary layer base
      INTEGER :: imld ! level of mixed-layer depth (pycnocline top)

      REAL(wp) :: ztgrad,zsgrad,zbgrad ! Temporary variables used to calculate pycnocline gradients
      REAL(wp) :: zugrad,zvgrad  ! temporary variables for calculating pycnocline shear

      REAL(wp) :: zhbl ! bl depth - grid
      REAL(wp) :: zhml ! ml depth - grid
      REAL(wp) :: zdh ! pycnocline depth - grid
      REAL(wp) :: zhbl_s,zhbl_t   ! temporary varisbles for BL depth
      REAL(wp) :: zthick
      REAL(wp) :: zt_bl,zs_bl,zu_bl,zv_bl,zrh_bl,zb_bl  ! averages over the depth of the blayer
      REAL(wp) :: zdt_bl,zds_bl,zdu_bl,zdv_bl,zdrh_bl,zdb_bl ! difference between blayer average and parameter at base of blayer
      REAL(wp) :: zt_ml,zs_ml,zu_ml,zv_ml,zrh_ml,zb_ml  ! averages over the depth of the mixed layer
      REAL(wp) :: zdt_ml,zds_ml,zdu_ml,zdv_ml,zdrh_ml,zdb_ml ! difference between mixed layer average and parameter at base of blayer
      REAL(wp) :: zdtdz_ext,zdsdz_ext,zdbdz_ext      ! External Gradients
      REAL(wp) :: zzeta_m,zgamma_b_nd,zgamma_s_nd,zgamma_t_nd,zgamma_b_f,zgamma_b
      REAL(wp) :: zwth_ent,zws_ent ! heat and salinity fluxes at the top of the pycnocline
      REAL(wp) :: zwb_ent, zwb_min
      REAL(wp) :: zuw_bse,zvw_bse ! momentum fluxes at the top of the pycnocline
      REAL(wp) :: zz0,zz1,zrad0,zradh,zradav
      REAL(wp) :: zbeta_d_sc,zbeta_v_sc

     ! Flux-gradient relationship variables

      REAL(wp) :: zl_c,zl_l,zl_eps  ! Used to calculate turbulence length scale.

      REAL(wp) :: zdifml_sc,zvisml_sc,zdifpyc_sc,zvispyc_sc,zbeta_d,zbeta_v ! Scales for eddy diffusivity/viscosity
      REAL(wp) :: zsc_wth_1,zsc_ws_1 ! Temporary scales used to calculate scalar non-gradient terms.
      REAL(wp) :: zsc_uw_1,zsc_uw_2,zsc_vw_1,zsc_vw_2 ! Temporary scales for non-gradient momentum flux terms.
! Tenporary veriables
      INTEGER :: inhml,jm
      REAL(wp) :: zari, zdhdt, zdhdt_2, zddhdt, ztau ! rates of change of OSBL depth and thickness of the pycnocline
      REAL(wp) :: znd,zznd,znd_d,zznd_ml,zznd_pyc,zznd_d, zdhoh  ! temporary non-dimensional depths used in various routines
      REAL(wp) :: ztemp, zpert, ztemp_2, ztemp_3, ztemp_4, ztemp_5 ! temporary variables

      INTEGER, parameter :: jp_tem=1, jp_sal=2, ibld_ext=0
      REAL, parameter :: pthird=1.0/3.0, p2third=2.0/3.0, zflageos=0.0
      ! set maximum boundary layer depth
      ! Qing Li, 20180409
      REAL(wp), parameter :: hbl_max=200.0
     !!--------------------------------------------------------------------

      ztemp=0.
      ztemp_2=0.
      ztemp_3=0.
      ztemp_4=0.
      ztemp_5=0.
      jpkm1=nlev-1
      ibld = 0
      imld = 0
      zBo  = 0._wp
      zBosol = 0._wp
      zustar = 0._wp
      zwstrl = 0._wp
      zvstr = 0._wp
      zwstrc = 0._wp
      zuw0 = 0._wp
      zvw0 = 0._wp
      zwth0 = 0._wp
      zws0  = 0._wp
      zwb0 = 0._wp
      zwthav = 0._wp
      zwsav = 0._wp
      zwbav = 0._wp
      zustke = 0._wp
      zla = 0._wp
      zcos_wind = 0._wp
      zsin_wind = 0._wp
      zhol = 0._wp
     ! mixed layer
      zt_bl = 0._wp
      zs_bl = 0._wp
      zu_bl = 0._wp
      zv_bl = 0._wp
      zrh_bl = 0._wp
      zt_ml = 0._wp
      zs_ml = 0._wp
      zu_ml = 0._wp
      zv_ml = 0._wp
      zrh_ml = 0._wp
      zdt_bl = 0._wp
      zds_bl = 0._wp
      zdu_bl = 0._wp
      zdv_bl = 0._wp
      zdrh_bl = 0._wp
      zdt_ml = 0._wp
      zds_ml = 0._wp
      zdu_ml = 0._wp
      zdv_ml = 0._wp
      zdrh_ml = 0._wp
      zwth_ent = 0._wp
      zws_ent = 0._wp
      zuw_bse = 0._wp
      zvw_bse = 0._wp
      zdtdz_ext = 0._wp
      zdsdz_ext = 0._wp
      zdbdz_ext = 0._wp
      zzeta_m=0._wp
      zgamma_b_nd=0._wp
      zgamma_s_nd=0._wp
      zgamma_t_nd=0._wp
      zgamma_b_f=0._wp
      zgamma_b=0._wp
      zdtdz_pyc(:) = 0._wp
      zdsdz_pyc(:) = 0._wp
      zdbdz_pyc(:) = 0._wp
      zdudz_pyc(:) = 0._wp
      zdvdz_pyc(:) = 0._wp
     ! Flux-Gradient arrays.
      zdifml_sc = 0._wp
      zvisml_sc = 0._wp
      zdifpyc_sc = 0._wp
      zvispyc_sc = 0._wp
      zbeta_d = 0._wp
      zbeta_v = 0._wp
      zsc_wth_1 =0._wp
      zsc_ws_1 = 0._wp
      zsc_uw_1 = 0._wp
      zsc_uw_2 = 0._wp
      zsc_vw_1 = 0._wp
      zsc_vw_2 = 0._wp

      zdiffut(:) = 0._wp
      zviscos(:) = 0._wp
      zghamt(:) = 0._wp
      zghams(:) = 0._wp
      zgamu(:) = 0._wp
      zgamv(:) = 0._wp
      ! zghamt_f(:) = 0._wp
      ! zghams_f(:) = 0._wp
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! Calculate boundary layer scales
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     ! Assume two-band radiation model for depth of OSBL
      zz0 =        rn_abs   * r1_rau0_rcp      ! surface equi-partition in 2-bands
      zz1 = ( 1. - rn_abs ) * r1_rau0_rcp
     ! Surface irradiance
      zrad0 = qsr * r1_rau0_rcp
     ! Irradiance at base of boundary layer
      zradh = qsr * ( zz0 * EXP( -hbl/rn_si0 ) + zz1 * EXP( -hbl/rn_si1) )
     ! Irradiance averaged over depth of the OSBL
      zradav = qsr * ( zz0 * ( 1.0 - EXP( -hbl/rn_si0 ) )*rn_si0 &
                          + zz1 * ( 1.0 - EXP( -hbl/rn_si1 ) )*rn_si1 ) / hbl
     ! Turbulent surface fluxes and fluxes averaged over depth of the OSBL
      zrhos    = rau0 * ( 1._wp + rhd(1) )
      zthermal = rab_n(1,jp_tem)
      zbeta    = rab_n(1,jp_sal)
     ! Surface Temperature flux for non-local term
      zwth0 =  -qns * r1_rau0_rcp
     ! Surface salinity flux for non-local term
      zws0 =  ( -( emp-rnf ) * tsn(1,jp_sal)                          &
              &             + sfx                                     ) * rcs
! 26/05/17 radiative buoyancy flux removes
! change to NEMO code 08/06/17
     ! Non radiative surface buoyancy force
      zwb0 = grav * r1_rau0 * zthermal * zwth0 - &
             grav * r1_rau0 * zbeta * zws0
     ! turbulent heat flux averaged over depth of OSBL
      zwthav = 0.5 * zwth0 - ( 0.5*( zrad0 + zradh ) - zradav )
     ! turbulent salinity flux averaged over depth of the OBSL
      zwsav = 0.5 * zws0
     ! turbulent buoyancy flux averaged over the depth of the OBSBL
      zwbav = grav * r1_rau0 * zthermal * zwthav - grav * r1_rau0 * zbeta * zwsav
      zuw0 = -utau * r1_rau0
      zvw0 = -vtau * r1_rau0
     !  Reference surface density = density at first T point level
      zrhos = rhop(1) + zflageos * rau0
     ! Friction velocity (zustar), at T-point : LMD94 eq. 2
!
! 25/5/17 Add minimum values to ustar and ustke to avoid floating exceptions
! Change to NEMO 08/06/17
      zustar = MAX( SQRT( SQRT( zuw0 * zuw0 + zvw0 * zvw0 ) ), 1.0e-8_wp)
     ! Stokes drift (zustke), at T-point
      zustke = MAX( SQRT( us_x**2 + us_y**2 ), 1.e-8_wp)    !   zustar/zla**2
     ! Langmuir velocity scale (zwstrl), at T-point
      zwstrl = ( zustar * zustar * zustke )**pthird
      zla = MIN( SQRT(zustar / zwstrl)**3, 4._wp )
     ! direction of shear stress
!
! 10/05/17 limit added to ensure that the Stokes part of gamu doesn't become infinite, or
! negative if the Langmuir number becomes to large.
! Changed NEMO code ????
!
      zcos_wind = -zuw0 / ( zustar * zustar )
      zsin_wind = -zvw0 / ( zustar * zustar )
     ! Langmuir velocity scale (zwstrl), at T-point
!
! 11/05/17 velocity scale that tends  to ustar as Langmuir number become large
! Changed NEMO code 08/06/17
!
      zvstr = ( zustar * zustar * zustke + ( 1.0 - EXP( -0.5 * zla**2 ) ) * zustar* zustar* zustar )**pthird
!
! 11.05/17 Use velocity scale zvstr in stability parameter zhol
! Changed NEMO code 08/06/17
!
      IF ( zwbav > 0.0) THEN
         zwstrc = ( 2.0 * zwbav * 0.9 * hbl )**pthird       ! approximate mixed-layer depth as 0.9 * hbl
         zhol = -0.9 * hbl / ( zvstr**3 / ( 2.0 * ( zwbav + epsln ) ) )
         lconv = .TRUE.
      ELSE
         zhol = -hbl / ( zvstr**3 / ( 2.0 * ( zwbav - epsln ) ) )
!         IF ( lconv ) hbl=hml
         lconv = .FALSE.
      ENDIF
!
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! Mixed-layer model - calculate averages over the boundary layer, and the change in the boundary layer depth
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      hbl = MAX(hbl, fs_depw(4))
      IF (hbl > hbl_max) hbl=hbl_max
      ibld=4
      DO jk = 5, jpkm1
         IF ( hbl >= fs_depw(jk) ) THEN
            ibld = jk
         ENDIF
      END DO
      zhbl=fs_depw(ibld)
      inhml = MAX( INT( dh / e3t_n(ibld) ) , 1  )
      imld = ibld - inhml
      zhml = fs_depw(imld)
!
! Calculate vertical averages over boundary layer and mixed layer
!
      CALL vert_av( ibld, zt_bl, zs_bl, zb_bl, zu_bl, zv_bl, zdt_bl, zds_bl, zdb_bl,&
             &      zdu_bl, zdv_bl )
      CALL vert_av( imld-1, zt_ml, zs_ml, zb_ml, zu_ml, zv_ml, zdt_ml, zds_ml, zdb_ml,&
             &      zdu_ml, zdv_ml )
! Calculate external gradients
      CALL grad_ext( zdtdz_ext, zdsdz_ext, zdbdz_ext )
      CALL calc_dhdt( zdhdt, zdhdt_2 )
!
      imld = ibld   ! use imld to hold previous boundary layer index
      ibld = 4
      zhbl_t = hbl + zdhdt * rn_rdt
!
! Not part of scheme, but used on OSMOSIS tests to limit OSBL depth in the winter.
!
      IF (zhbl_t > hbl_max) THEN
         zhbl_t=hbl_max
         zdhdt=(zhbl_t-hbl)/rn_rdt
      ENDIF

      DO jk = 5, jpkm1
        IF ( zhbl_t >= fs_depw(jk) ) THEN
           ibld = jk
        ENDIF
      END DO
!
      CALL depth_tstep( zdhdt, zdhdt_2 )
      zhbl = fs_depw(ibld)
! Recalculate averages over boundary layer after depth updated
      CALL vert_av( ibld, zt_bl, zs_bl, zb_bl, zu_bl, zv_bl, zdt_bl, zds_bl, zdb_bl, &
             & zdu_bl, zdv_bl )
!
      CALL grad_ext( zdtdz_ext, zdsdz_ext, zdbdz_ext )
!
      CALL pyc_depth( dh, zdh )
!
! Check that pycnocline exists if not encroachment
!
      IF ( zdb_bl > 0.0 ) THEN
         zwth_ent = zwb_ent * zdt_bl / zdb_bl
         zws_ent = zwb_ent * zds_bl / zdb_bl
      ELSE
         zwth_ent=0.0
         zws_ent=0.0
      ENDIF
     ! Average over the depth of the mixed layer in the convective boundary layer
      CALL vert_av( imld-1, zt_ml, zs_ml, zb_ml,zu_ml, zv_ml, zdt_ml, zds_ml, zdb_ml, &
             & zdu_ml, zdv_ml )
     ! rotate mean currents and changes onto wind align co-ordinates
      CALL vel_rot( zcos_wind, zsin_wind, zu_ml, zv_ml, zdu_ml, zdv_ml )
      CALL vel_rot( zcos_wind, zsin_wind, zu_bl, zv_bl, zdu_bl, zdv_bl )
     ! Calculate profiles in the pycnocline.
      CALL pyc_profiles( zdtdz_pyc, zdsdz_pyc, zdbdz_pyc, &
                   &     zdudz_pyc, zdvdz_pyc )


       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ! Eddy viscosity/diffusivity and non-gradient terms in the flux-gradient relationship
       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF ( lconv ) THEN
         zuw_bse = -0.0075*((zvstr**3+0.5*zwstrc**3)**pthird*zdu_ml + 1.5*zustar**2*(zhbl-zhml)/(zhml*MIN(zla**(8./3.),1._wp)))
         zvw_bse = 0.01*(-(zvstr**3+0.5*zwstrc**3)**pthird*zdv_ml+2.0*ff*zustke*delta*zla)
      ELSE
         zwth_ent = -2.0*zwthav * ( ( 1.0 - 0.8 ) - ( 1.0 - 0.8)**(3.0/2.0) )
         zws_ent = -2.0*zwsav * ( (1.0-0.8 ) - ( 1.0 - 0.8 )**(3.0/2.0) )
         zuw_bse = 0.0
         zvw_bse = 0.0
      ENDIF
!
      IF ( lconv ) THEN                                                ! Replaces NEMO version WHERE-ELSEWHERE-ENDWHERE
         zdifml_sc = zhml * ( zvstr**3 + 0.5 * zwstrc**3 )**pthird
         zvisml_sc = zdifml_sc
         zdifpyc_sc = 0.165 * ( zvstr**3 + 0.5 * zwstrc**3 )**pthird * ( zhbl - zhml )
         zvispyc_sc = 0.142 * ( zvstr**3 + 0.5 * zwstrc**3 )**pthird * ( zhbl - zhml )
         zbeta_d_sc = 1.0 - (zdifpyc_sc/(0.8 * zdifml_sc))**p2third
         zbeta_v_sc = 1.0 - ( 2.0 * 0.142 /0.375 * (zhbl - zhml ) / zhml )
      ELSE
         zdifml_sc = zvstr * zhbl * EXP ( -( zhol / 0.6 )**2 )
         zvisml_sc = zvstr * zhbl * EXP ( -( zhol / 0.6 )**2 )
      ENDIF
!
      IF ( lconv ) THEN
         DO jk = 2, imld   ! mixed layer diffusivity
            zznd_ml = fs_depw(jk) / zhml
                    !
            zdiffut(jk) = 0.8 * zdifml_sc * zznd_ml * &
                                   ( 1.0 - zbeta_d_sc * zznd_ml )**1.5
                       !
            zviscos(jk) = 0.375 * zvisml_sc * zznd_ml * &
                              ( 1.0 - zbeta_v_sc * zznd_ml ) * ( 1.0 - 0.5 * zznd_ml**2 )
         END DO
             ! pycnocline - linear profile
! ??/05/17
! Change NEMO code 08/06/17
         IF ( zdh > 0.0 ) THEN
            zgamma_b =6.0
            DO jk = imld+1 , ibld
               zznd_pyc = -( fs_depw(jk) - zhml ) / ( zhbl - zhml )
                       !
               zdiffut(jk) = zdifpyc_sc * EXP( zgamma_b * zznd_pyc )
                     !
               zviscos(jk) = zvispyc_sc * EXP( zgamma_b * zznd_pyc )

            END DO
         ENDIF
         IF ( ibld_ext .eq. 0 ) THEN
! 23/10/17 : Set eddy diffusivity/viscosity to zero at base of boundary layer when no interface layer
!            included. Prevents unwanted effects on prediction of boundary layer depth
            zdiffut(ibld) = 0._wp
            zviscos(ibld) = 0._wp
         ELSE
            zdiffut(ibld) = zdhdt * (hbl-fs_dept(ibld-1))
            zviscos(ibld) = zdhdt * (hbl-fs_dept(ibld-1))
         ENDIF
      ELSE
        ! stable conditions
         DO jk = 2, ibld
            zznd_ml = fs_depw(jk) / zhbl
            zdiffut(jk) = 0.75 * zdifml_sc * zznd_ml * ( 1.0 - zznd_ml )**1.5
            zviscos(jk) = 0.375 * zvisml_sc * zznd_ml * (1.0 - zznd_ml) * ( 1.0 - zznd_ml**2 )
         END DO
         IF ( ibld_ext .eq. 0 ) THEN
! 23/10/17 : Set eddy diffusivity/viscosity to zero at base of boundary layer when no interface layer
!            included. Prevents unwanted effects on prediction of boundary layer depth
            zdiffut(ibld) = 0._wp
            zviscos(ibld) = 0._wp
         ELSE
            zdiffut(ibld) = zdhdt * e3w_n(ibld)
            zviscos(ibld) = zdhdt * e3w_n(ibld)
         ENDIF
      ENDIF   ! end if ( lconv )
!
      !
      ! calculate non-gradient components of the flux-gradient relationships
       !
! Stokes term in scalar flux, flux-gradient relationship
      IF ( lconv ) THEN
         zsc_wth_1 = zwstrl**3 * zwth0 / ( zvstr**3 + 0.5 * zwstrc**3 )
        !
         zsc_ws_1 = zwstrl**3 * zws0 / ( zvstr**3 + 0.5 * zwstrc**3 )
      ELSE
         zsc_wth_1 = 2.0 * zwthav * SQRT ( SQRT( delta / hbl ) ) ! ?
        !
         zsc_ws_1 = 2.0 * zwsav * SQRT ( SQRT( delta / hbl ) ) ! ? The delta/hbl factor ?
      ENDIF
      IF ( lconv ) THEN
         DO jk = 2, imld
            zznd_d = fs_depw(jk) / delta
            zghamt(jk) = zghamt(jk) + 1.35 * EXP ( -zznd_d ) * &
                            ( 1.0 - EXP ( -2.0 * zznd_d ) ) * zsc_wth_1
            !
            zghams(jk) = zghams(jk) + &
                            1.35 * EXP ( -zznd_d ) * ( 1.0 - EXP ( -2.0 * zznd_d ) ) &

                                 * zsc_ws_1
         END DO ! end jk loop
      ELSE     ! else for if (lconv)
 ! Stable conditions
         DO jk = 2, ibld
            zznd_d=fs_depw(jk) / delta
            zghamt(jk) = zghamt(jk) + 1.5 * EXP ( -0.9 * zznd_d ) * ( 1.0 - EXP ( -4.0 * zznd_d ) ) * zsc_wth_1
                      !
            zghams(jk) = zghams(jk) + 1.5 * EXP ( -0.9 * zznd_d ) * ( 1.0 - EXP ( -4.0 * zznd_d ) ) * zsc_ws_1
         END DO
      ENDIF               ! endif for check on lconv
! Stokes term in flux-gradient relationship
      IF ( lconv ) THEN
!
! problem here 10/05/17. The La factor can become large and/or negative for large Langmuir
! number. Causes problems with currents. Testing limiting Langmuir number.
! Changed coefficient 1.5 to 1.0 10/05/17 added NEMO 08/06/17
!
         zsc_uw_1 = ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * zustke / &
                    & MAX( ( 1.0 - 1.0 * 6.5 * zla**(8.0/3.0) ), 0.2_wp )
         zsc_uw_2 = ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * zustke / &
                         MIN( zla**(8.0/3.0), 0.12_wp )
         zsc_vw_1 = ff * zhml * zustke**3 * MIN ( zla**(8.0/3.0), 0.12_wp ) / ( zvstr**3 + 0.5 * zwstrc**3 )**(2.0/3.0)
!         zsc_vw_1 = ff * zhml * zustke**(5.0/3.0) * zustar**(4.0/3.0)/ ( zvstr**3 + 0.5 * zwstrc**3 )**(2.0/3.0)
      ELSE
         zsc_uw_1 = zustar**2
         zsc_vw_1 = ff * zhbl * zustke**3.0 * MIN ( zla**(8./3.), 0.12_wp ) /zvstr**2
!              zsc_vw_1 = ff * zhbl * zustke**(5.0/3.0) * zustar**(4.0/3.0) / zvstr**2
      ENDIF
      IF ( lconv ) THEN
         DO jk = 2, imld
            zznd_d = fs_depw(jk) / delta
            zgamu(jk) = zgamu(jk) + ( -0.05 * EXP ( -0.4 * zznd_d ) * zsc_uw_1 &
                                          + 0.00125 * EXP ( -zznd_d ) * zsc_uw_2 ) * &
                                         ( 1.0 - EXP ( -2.0 * zznd_d ) )
!
            zgamv(jk) = zgamv(jk) - 0.65 * 0.15 * EXP ( -zznd_d ) * zsc_vw_1  * &
                                              ( 1.0 - EXP ( -2.0 * zznd_d ) )
         END DO   ! end jk loop
      ELSE
! Stable conditions
         DO jk = 2, imld
            zznd_d = fs_depw(jk) / delta
            zgamu(jk) = zgamu(jk) - 0.75 * 1.3 * EXP ( -0.5 * zznd_d ) * ( 1.0 - EXP ( -4.0 * zznd_d ) ) * zsc_uw_1
            IF ( zznd_d <= 1.0 ) THEN
               zgamv(jk) = zgamv(jk) + 0.0*2.15 * zznd_d * (1.0 - zznd_d )**1.5 * EXP ( -2.0 * zznd_d ) * zsc_vw_1
            ENDIF
         END DO   ! end jk loop
      ENDIF

! Buoyancy term in flux-gradient relationship [note : includes ROI ratio (X0.3) and pressure (X0.5)]

      IF ( lconv ) THEN
         zsc_wth_1 =  zwbav * zwth0 * ( 1.0 + EXP ( 0.2 * zhol ) ) / ( zvstr**3 + 0.5 * zwstrc**3 )
         zsc_ws_1 = zwbav * zws0 * ( 1.0 + EXP ( 0.2 * zhol ) ) / ( zvstr**3 + 0.5 * zwstrc**3 )
      ELSE

      ENDIF
      IF ( lconv ) THEN
         DO jk = 2, imld
            znd=fs_depw(jk) / zhml
         ! calculate turbulent length scale
            zl_c = 0.9 * ( 1.0 - EXP ( -7.0 * ( znd - znd**3 / 3.0 ) ) ) * &
                         ( 1.0 - EXP ( -15.0 * ( 1.1 - znd ) ) )

            zl_l = 2.0 * ( 1.0 - EXP ( -2.0 * ( znd - znd**3 / 3.0 ) ) ) * &
                         (1.0 -  EXP ( -5.0 * ( 1.0 -znd ) ) ) &
                       * ( 1.0 + delta / zhml )
            zl_eps = zl_l + ( zl_c - zl_l ) / &
                            ( 1.0 + EXP ( 3.0 * LOG10 ( - zhol ) ) ) ** (3.0/2.0)
         ! non-gradient buoyancy terms
            zghamt(jk) = zghamt(jk) + 0.3 * 0.5 * zsc_wth_1 * zl_eps * zhml / ( 0.15 + znd )
            zghams(jk) = zghams(jk) + 0.3 * 0.5 * zsc_ws_1 * zl_eps * zhml / ( 0.15 + znd )
         END DO
      ELSE
      ! ************************** add stable buoyancy term
      ENDIF

      IF ( lconv ) THEN
         zsc_uw_1 = -zwb0 * zustar**2 * zhml / ( zvstr**3 + 0.5 * zwstrc**3 )
         zsc_uw_2 = zwb0 * zustke * zhml / ( zvstr*3 + 0.5 * zwstrc**3 )**(2.0/3.0)
      ELSE
      !********************** add stable buoyancy term
      ENDIF
      IF ( lconv ) THEN
         DO jk = 2 , imld
            zznd_d = fs_depw(jk) / delta
            znd = -( fs_depw(jk) - zhml ) / zhml
            zgamu(jk) = zgamu(jk) + 0.3 * 0.5 * ( zsc_uw_1 + 0.125 * EXP( -0.5 * zznd_d ) &
                                     * ( 1.0 - EXP ( -0.5 * zznd_d) ) * zsc_uw_2 )
         END DO  ! jk loop
      ELSE
          ! stable conditions
      ENDIF

! Transport term in flux-gradient relationship [note : includes ROI ratio (X0.3) ]

      IF ( lconv ) THEN
         zsc_wth_1 = zwth0
         zsc_ws_1 = zws0
      ELSE
         zsc_wth_1 = 2.0 * zwthav
         zsc_ws_1 = 2.0 * zws0
      ENDIF
      IF ( lconv ) THEN
         DO jk = 2, imld
            zznd_ml=fs_depw(jk) / zhml
            zghamt(jk) = zghamt(jk) + 0.3 * zsc_wth_1 &
                          * ( -2.0 + 2.75 * ( ( 1.0 + 0.6 * zznd_ml**4 ) &
                             - EXP ( -6.0 *zznd_ml ) ) )  &
                          * ( 1.0 -EXP ( -15.0 * ( 1.0 - zznd_ml )   ) )
                  !
            zghams(jk) = zghams(jk) + 0.3 * zsc_ws_1 * &
                          ( -2.0 + 2.75 * ( ( 1.0 + 0.6 * zznd_ml**4 ) &
                            - EXP ( -6.0 *zznd_ml ) ) )  &
                        * ( 1.0 -EXP ( -15.0 * ( 1.0 - zznd_ml )   ) )
         END DO
      ELSE
         IF( zdhdt .gt. 0.) THEN
! stable transport term only when boundary layer depth increasing
            DO jk = 2, ibld
               zznd_d = -fs_depw(jk) / delta
               zznd = -fs_depw(jk) / zhbl
               zghamt(jk) = zghamt(jk) + zsc_wth_1 * 0.3 * ( -4.06 * EXP( 2.0 * zznd_d ) * &
                                 ( 1.0 - EXP( 4.0 * zznd_d ) ) + &
                                 7.5 * EXP( -10.0 * ( 0.95+zznd )**2 ) * ( 1.0 + zznd ) )
               zghams(jk) = zghams(jk) + zsc_ws_1 * 0.3 *( -4.06 * EXP( 2.0 * zznd_d ) * &
                                     ( 1.0 - EXP( 4.0 * zznd_d ) ) + &
                                    7.5 * EXP( -10.0 * ( 0.95+zznd )**2 ) * ( 1.0 + zznd ) )
            ENDDO
         ENDIF
      ENDIF

      IF ( lconv ) THEN
         zsc_uw_1 = zustar**2
         zsc_vw_1 = ff * zustke * zhml
      ELSE
         zsc_uw_1 = zustar**2
         zsc_uw_2 = (2.25 - 3.0 * ( 1.0 - EXP( -1.25 * 2.0 ) ) ) * ( 1.0 - EXP( -4.0 * 2.0 ) ) * zsc_uw_1
         zsc_vw_1 = ff * zustke * zhbl
         zsc_vw_2 = -0.11 * SIN( 3.14159 * ( 2.0 + 0.4 ) ) * EXP(-( 1.5 + 2.0 )**2 ) * zsc_vw_1
      ENDIF
      IF ( lconv ) THEN
         DO jk = 2, imld
            zznd_ml = fs_depw(jk) / zhml
            zznd_d = fs_depw(jk) / delta
            znd = -( fs_depw(jk) - zhml ) / zhml
            zgamu(jk) = zgamu(jk) - 0.3 * ( -2.0 + 2.5 * ( 1.0 + 0.1 * zznd_ml**4 ) - &
                                         EXP ( -8.0 * zznd_ml ) ) * zsc_uw_1
              !
            zgamv(jk) = zgamv(jk) + 0.3 * 0.1 * ( EXP( -zznd_d ) + &
                                           EXP( -5.0 * ( 1.0 - zznd_ml ) ) ) * zsc_vw_1
         END DO
      ELSE
         DO jk = 2, ibld
            zznd = fs_depw(jk) / zhbl
            zznd_d = fs_depw(jk) / delta
            zgamv(jk) = zgamv(jk) +  0.3 * 0.15 * SIN( 3.14159 * ( 0.65 * zznd_d ) ) * &
                               EXP( -0.25 * zznd_d**2 ) * zsc_vw_1
                !
            zgamv(jk) = zgamv(jk) + 0.3 * 0.15 * EXP( -5.0 * ( 1.0 - zznd ) ) * ( 1.0 - EXP( -20.0 * ( 1.0 - zznd ) ) ) * zsc_vw_2
            IF ( zznd_d <= 2.0 ) THEN
               zgamu(jk) = zgamu(jk) + 0.5 * 0.3 * ( 2.25 - 3.0 * ( 1.0 - &
                                                    EXP( - 1.25 * zznd_d ) ) * &
                                              ( 1.0 - EXP( -2.0 * zznd_d ) ) ) * zsc_uw_1
            ELSE
               zgamu(jk) = zgamu(jk) + 0.5 * 0.3 * ( 1.0 - EXP( -5.0 * ( 1.0- zznd ) ) ) * zsc_uw_2
            ENDIF
         END DO
      ENDIF

      IF ( lconv ) THEN
!
! make velocity non-gradient terms go to  zero at the base of the nixed layer
!
         DO jk = 2, ibld
            znd = -( fs_depw(jk) - zhml ) / zhml
            IF ( znd >= 0.0 ) THEN
               zgamu(jk) = zgamu(jk) * ( 1.0 - EXP ( -30.0 * znd**2 ) )
               zgamv(jk) = zgamv(jk) * ( 1.0 - EXP ( -30.0 * znd**2 ) )
            ELSE
               zgamu(jk) = 0._wp
               zgamv(jk) = 0._wp
            ENDIF
         END DO
      ELSE
!
! make velocity non-gradient terms go to  zero at the base of the boundary layer
!
         DO jk = 2, ibld
            znd = -( fs_depw(jk) - zhbl ) / zhbl
!               IF ( znd >= 0.0 ) THEN
            zgamu(jk) = zgamu(jk) * ( 1.0 - EXP ( -1.0 * znd ) )
            zgamv(jk) = zgamv(jk) * ( 1.0 - EXP ( -0.5 * znd ) )
!               ELSE
!                zgamu(jk) = 0._wp
!                zgamv(jk) = 0._wp
!               ENDIF
         END DO
      ENDIF

! pynocline contributions

      DO jk= 2, ibld
         znd = fs_depw(jk)/zhbl
         zghamt(jk) = zghamt(jk) + zdiffut(jk) * zdtdz_pyc(jk)
         zghams(jk) = zghams(jk) + zdiffut(jk) * zdsdz_pyc(jk)
         zgamu(jk) = zgamu(jk) + zviscos(jk) * zdudz_pyc(jk)
! taken out because it can lead to large spike in shear at base of BL, presumably when
! zdb_bl is small. Need to check so left in as a comment.
!IF ( znd .lt. 1.0 ) THEN
!        zgamu(jk) = zgamu(jk) + 50.0 * zla**(8./3.) * zustar**2 * ( 1.0 - znd )**(7./4.) * &
 !                                                ( zhbl / zdb_bl ) * zdbdz_pyc(jk)
!ENDIF
         zgamv(jk) = zgamv(jk) + zviscos(jk) * zdvdz_pyc(jk)
      END DO
! Entrainment flux contribution.
      IF (lconv ) THEN
         DO jk = 1, imld-1
            znd=fs_depw(jk) / zhml
            zghamt(jk) = zghamt(jk) + zwth_ent * znd
            zghams(jk) = zghams(jk) + zws_ent * znd
            zgamu(jk) = zgamu(jk) + zuw_bse * znd
            zgamv(jk) = zgamv(jk) + zvw_bse * znd
         END DO
         DO jk = imld, ibld-1
            znd = -( fs_depw(jk) - zhml ) / ( zhbl - zhml )
            zghamt(jk) = zghamt(jk) + zwth_ent * ( 1.0 + znd )
            zghams(jk) = zghams(jk) + zws_ent * ( 1.0 + znd )
            zgamu(jk) = zgamu(jk) + zuw_bse * ( 1.0 + znd )
            zgamv(jk) = zgamv(jk) + zvw_bse * ( 1.0 + znd )
         END DO
      ENDIF
      zghamt(ibld+ibld_ext)=0._wp
      zghams(ibld+ibld_ext)=0._wp
      zgamu(ibld+ibld_ext)=0._wp
      zgamv(ibld+ibld_ext)=0._wp

     ! Interpolate eddy viscosity onto u and v points
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Code removed from GOTM  u,v grid interpolation
!                         lateral boundary layers
! transfer from local zviscos, zdiffut to zavmu, zavmv and zavt
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      DO jk = 2, jpkm1
         zavt(jk) = MAX( zdiffut(jk), zavt(jk) )
         zavmu(jk) = MAX( zviscos(jk), zavmu(jk) )
         zavmv(jk) = MAX( zviscos(jk), zavmv(jk) )
      END DO

      DO jk = 2, jpkm1
        !
        !  Minimum value on the eddy diffusivity
        ! ----------------------------------------
         zavt(jk) = MAX( zavt(jk), zavtb(jk) )

        !
        ! Minimum value on the eddy viscosity
        ! ----------------------------------------
         zavmu(jk) = MAX( zavmu(jk), zavmb(jk) )
         zavmv(jk) = MAX( zavmv(jk), zavmb(jk) )
      END DO
!
      DO jk = ibld, jpkm1
         zgamu(jk) = zgamu(jk) + zavmu(jk) * zdudz_pyc(jk)
      END DO
!
! rotate non-gradient velocity terms back to model reference frame

      DO jk = 2, ibld
         ztemp = zgamu(jk)
         zgamu(jk) = zgamu(jk) * zcos_wind - zgamv(jk) * zsin_wind
         zgamv(jk) = zgamv(jk) * zcos_wind + ztemp * zsin_wind
      END DO

CONTAINS

!BOP
      SUBROUTINE vert_av( jnlev_av, zt, zs, zb, zu, zv, zdt, zds, zdb, zdu, zdv )
!
!  Calculates vertical averages.
      IMPLICIT NONE
!
      INTEGER :: jnlev_av
      REAL(wp) :: zt, zs, zu, zv, zb, zdt, zds, zdb, zdu, zdv

      INTEGER :: jm
      REAL(wp) :: zthick,zthermal,zbeta

      zt = 0.
      zs = 0.
      zu = 0.
      zv = 0.
      zthick=0.0
!
      zthermal = rab_n(1,jp_tem)
      zbeta    = rab_n(1,jp_sal)
!
      DO jm = 2, jnlev_av
         zthick=zthick+e3t_n(jm)
         zt   = zt  + e3w_n(jm) * tsn(jm,jp_tem)
         zs   = zs  + e3w_n(jm) * tsn(jm,jp_sal)
         zu   = zu  + e3w_n(jm) * ub(jm)         ! interpolation from NEMO u-grid removed
         zv   = zv  + e3w_n(jm) * vb(jm)         ! interpolation from NEMO v-grid removed
      END DO
      zt = zt / zthick
      zs = zs / zthick
      zu = zu / zthick
      zv = zv / zthick
      zb = grav * r1_rau0 * zthermal * zt - grav * r1_rau0 * zbeta * zs
      zdt = zt - tsn(ibld+ibld_ext,jp_tem)
      zds = zs - tsn(ibld+ibld_ext,jp_sal)
      zdu = zu - ub(ibld+ibld_ext)
      zdv = zv - vb(ibld+ibld_ext)
      zdb = grav * r1_rau0 * zthermal * zdt - grav * r1_rau0 * zbeta * zds

      END SUBROUTINE vert_av
!EOC
!
!BOP
      SUBROUTINE grad_ext( zdtdz, zdsdz, zdbdz )
!
! Subroutine to calculate external gradients
!
      IMPLICIT NONE
      REAL(wp) :: zdtdz, zdsdz, zdbdz
      REAL(wp) :: zthermal,zbeta

      zthermal = rab_n(1,jp_tem)
      zbeta    = rab_n(1,jp_sal)
!
      zdtdz = -( tsn(ibld+ibld_ext+1,jp_tem) - tsn(ibld+ibld_ext,jp_tem) ) / &
         &           ( fs_dept(ibld+ibld_ext+1) - fs_dept(ibld+ibld_ext) )
      zdsdz = -( tsn(ibld+ibld_ext+1,jp_sal) - tsn(ibld+ibld_ext,jp_sal) ) / &
         &           ( fs_dept(ibld+ibld_ext+1) - fs_dept(ibld+ibld_ext) )
      zdbdz = grav * r1_rau0 * zthermal * zdtdz - grav * r1_rau0 * zbeta * zdsdz

      END SUBROUTINE grad_ext
!EOC
!BOP
      SUBROUTINE vel_rot( zcos_w, zsin_w, zu, zv, zdu, zdv )
!
! Rotate velocity components to frame aligned with the surface stress
!
      IMPLICIT NONE
      REAL(wp) :: zcos_w, zsin_w                ! cos and sin of rotation angle
      REAL(wp) :: zu, zv                        ! Components of current
      REAL(wp) :: zdu, zdv                      ! Change in current across pycnocline

      REAL(wp) :: ztemp
    !
    ! rotate mean currents and changes onto wind align co-ordinates
    !

      ztemp = zu
      zu = zu * zcos_wind + zv * zsin_wind
      zv = zv * zcos_wind - ztemp * zsin_wind
      ztemp = zdu
      zdu = zdu * zcos_wind + zdv * zsin_wind
      zdv = zdv * zcos_wind - ztemp * zsin_wind

      END SUBROUTINE vel_rot
!EOC
!BOP
      SUBROUTINE calc_dhdt( zdhdt, zdhdt_2 )
!
! Subroutine calculates dhdt
!
      IMPLICIT NONE
      REAL(wp) :: zdhdt, zdhdt_2

      INTEGER :: jj, ji
      REAL(wp) :: zgamma_b_nd, zgamma_dh_nd, zpert
      REAL(wp) :: zvel_max
      REAL(wp) :: zzeta_m = 0.3
      REAL(wp) :: zgamma_b = 2.0
      REAL(wp) :: zdhod = 0.1
      REAL(wp) :: alpha_bc = 0.5
      REAL(wp) :: cori=1.09e-4
      zdhdt = 0._wp
      zdhdt_2 = 0._wp
!
      IF ( lconv ) THEN    ! Convective
!
         zwb_ent = -( 2.0 * 0.2 * zwbav + &
                (0.15 * ( 1.0 - EXP( -1.5 * zla ) ) * zustar**3  + &
                 0.03 * zwstrl**3 ) &
                 / hbl )
         zwb_min =  dh * zwb0 / zhml + zwb_ent

         zvel_max = -0.1 * ( 1.0 + 1.0 * ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * rn_rdt / hbl )
         zvel_max = zvel_max * zwb_ent / ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird
         IF ( zdb_bl > 0.0 .and. zdbdz_ext > 0.0) THEN
            zdh = zhbl - zhml
            zdhdt = -wb(ibld) - zwb_ent / ( zvel_max + MAX(zdb_bl,1.0e-15_wp) )
            IF ( zdb_ml > 0._wp ) THEN
               zzeta_m = 0.3
               zgamma_b_nd = zdbdz_ext * zhml / zdb_ml
               zgamma_dh_nd =  zdbdz_ext * zdh / zdb_ml
               zdhdt_2 = MAX( -alpha_bc / (4.0 * zgamma_b) * (alpha_bc + zgamma_dh_nd ) * ( 1.0  - &
                  & SQRT( 3.14159 / ( 4.0 * zgamma_b ) * zdhoh ) ) * &
                  & ( zwb0 - ( 1.0 - 1.0 / alpha_bc * zgamma_b_nd ) * zwb_min ) * zdh / zhml, 0.0_wp )
!
! Effect of finite depth of pycnocline must increase growth of OSBL. Equation for the effect
! derived and tested for quasi-steady conditions, hence restriction.
!
! since zdhdt_2 not well contrained need to limit magnitude
                  IF ( zdhdt_2/( zvel_max + MAX(zdb_bl,1.0e-15_wp) ) <= 0.2*zdhdt ) &
                    &         zdhdt = zdhdt + zdhdt_2/ ( zvel_max + MAX( zdb_bl, 1.0e-15_wp) )

            ENDIF
         ELSE
            zdhdt = -wb(ibld) - zwb_ent / zvel_max
         ENDIF
      ELSE                        ! Stable
         zari = MAX( zdb_bl, -zwbav / zvstr ) * hbl / zvstr**2
         zdhdt = ( 0.06 + 0.52 * zhol / 2.0 ) * zwstrl**3 / hbl
         zdhdt = zdhdt + zwbav !* ( 1.0 - 0.25 * ( 1.0 - exp(-0.5*( zari - 1.0 ) ))* &
!& (1.0 - exp(-4.0 * (hbli / hbl - 1.0 ) ) ) )
         IF ( zdhdt < 0 ) THEN
            zpert = 2.0 * ( 1.0 + 0.5*2.0 * ( zvstr * rn_rdt / hbl ) ) * zvstr**2 / hbl
         ELSE
            IF ( zdb_bl <= 0.0 ) THEN
               zpert = 2.0 * ( 1.0 + 2.0 * ( zvstr * rn_rdt / hbl ) ) * zvstr**2 / hbl
            ELSE
               zpert = 2.0 * ( 1.0 + 2.0 * ( zvstr * rn_rdt / hbl ) ) * zvstr**2 / hbl + zdb_bl
            ENDIF
         ENDIF
         zdhdt = -wb(ibld) + 2.0 * zdhdt / zpert
      ENDIF    ! IF ( lconv )

      END SUBROUTINE calc_dhdt
!EOC
!BOP
      SUBROUTINE depth_tstep( zdhdt, zdhdt_2 )

      IMPLICIT NONE
      REAL(wp) :: zdhdt     ! rate of change of OSBL thickness
      REAL(wp) :: zdhdt_2   ! correction for pycnocline thickness
!
      INTEGER :: jm,jk
      REAL(wp) :: zhbl_s
      REAL(wp) :: zdb
      REAL(wp) :: zvel_max
!
! Timestep through model levels, taking account of buoyancy.
!
      IF ( ibld-imld > 1 ) THEN
!
! boundary layer depth changes by more than one level
!
         zhbl_s = hbl
         jm = imld
         IF ( lconv ) THEN
!
! Unstable
!
            zvel_max = 0.1 * ( 1.0 + 1.0 * ( zvstr**3 + 0.5 * zwstrc**3 )**pthird * rn_rdt / hbl )
            DO jk = imld, ibld
               zdb = max( grav * r1_rau0 * ( zthermal * ( zt_bl - tsn(jm,1) ) - &
                                      zbeta * ( zs_bl - tsn(jm,2) ) ) , 1.0e-10_wp )
!
! change 09/05/17 one level is maximum step per go around loop. prevents growing beyond
! stable layer.
! Added to NEMO code - 08/06/17
!
               zdb = zdb - &
                  zvel_max * zwb_ent / ( zvstr**3 + 0.5 * zwstrc**3 )**(1.0/3.0)
               zhbl_s=zhbl_s + min((-zwb_ent + zdhdt_2) / zdb * rn_rdt / &
         &                   FLOAT(ibld - imld ), e3w_n(jk))
               IF ( zhbl_s >= fs_depw(jm+1) ) jm = jm + 1
            END DO
            hbl = zhbl_s
            ibld = jm
!
            hbli=hbl
         ELSE
!
! Stable
!
            DO jk = imld, ibld
               zdb = max( grav * r1_rau0 * ( zthermal * ( zdt_bl - tsn(jm,1) ) - &
                                      zbeta * ( zds_bl - tsn(jm,2) ) ) , 0.0_wp ) + &
                     2.0 * ( 1.0 + 2.0 * ( zvstr * rn_rdt / hbl ) ) * zvstr**2 / zhbl_s
               zari = MAX( zdb_bl, -zwbav / zvstr ) * hbl / zvstr**2
               zdhdt = ( 0.06 + 0.52 * zhol / 2.0 ) * zwstrl**3 / hbl
               zdhdt = zdhdt + zwbav! * ( 1.0 - 0.25 * ( 1.0 - exp(-0.5*( zari - 1.0 ) )) * &
!& (1.0 - exp(-4.0 * (hbli / hbl - 1.0 ) ) ) )
               zhbl_s=zhbl_s + min( zdhdt / zdb * rn_rdt / &
         &                   FLOAT(ibld - imld ), e3w_n(jk))
               IF ( zhbl_s >= fs_depw(jm+1) ) jm = jm + 1
            END DO
            hbl = zhbl_s
            ibld = jm
            IF ( hbl < hbli ) hbli = hbl
         ENDIF
      ELSE                               ! ibld-imld <= 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 03/02/17
         hbl = MAX(zhbl_t, fs_depw(4))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 03/02/17
         IF ( lconv ) THEN
            hbli = hbl
         ELSE
            IF ( hbl > hbli ) hbli = hbl  ! when stable blayer grows above hbli update hbli
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 03/02/17
      ENDIF
      END SUBROUTINE depth_tstep
!EOC
!BOP
      SUBROUTINE pyc_depth( dh, zdh )
!
      IMPLICIT NONE
      REAL(wp) :: dh, zdh     ! pycnocline thickness.
!
      INTEGER :: inhml
      REAL(wp) :: ztau        ! relaxation timescale for pycnocline.
      REAL(wp) :: zddhdt      ! rate of change of pycnocline thickness.

      IF ( lconv ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!07/02/17
!
! 25/05/17 test on zdb_bl added
! change to NEMO code 08/06/17
!
         IF ( zdb_bl > 0.0 ) THEN
            IF ( ( zwstrc / zvstr )**3 <= 0.5 ) THEN
               zari = 4.5 /( ( zdb_bl * hbl / &
                  ( zvstr**2 )) + 0.01 )
            ELSE
               zari = 4.5 /( ( zdb_bl * hbl / &
                  ( zwstrc**2 )) + 0.01 )
            ENDIF

            ztau = -dh * ( zb_ml - zb_bl )/ zwb_min
            zddhdt = MIN (zari, 0.2_wp) * hbl
            IF ( ztau <= 0.0 ) ztau = hbl / ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird
         ELSE
            ztau = hbl / ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird
            zddhdt = 0.2 * hbl
         ENDIF
         dh = dh * exp(-rn_rdt / ztau ) + zddhdt * ( 1.0 - exp(-rn_rdt / ztau ) )
         hml = hbl - dh
         inhml = MAX( INT( dh / e3t_n(ibld) ) , 1  )
         imld = ibld - inhml
         zhml = fs_depw(imld)
         zdh = zhbl - zhml
      ELSE
      ! pycnocline depth fraction of boundary layer depth
!
         IF ( zdhdt >= 0.0 ) THEN
!
! stable but growing OSBL
!
            IF ( zdb_bl > 0 ) THEN
               zari = MIN( 4.5 /( ( zdb_bl * hbl / &
                    ( zwstrl**2 )) + 0.01 ), 0.2_wp)
               ztau = hbl / zwstrl
               zddhdt = MIN(zari, 0.2_wp) * hbl
            ELSE
!
! set pycnocline depth to 0.2*hbl
!
               ztau = hbl / zwstrl
               zddhdt = 0.2 * hbl
            ENDIF
! Use exponential relaxation to avoid problems withlong timesteps.
            dh = dh * exp(-rn_rdt / ztau ) + zddhdt * ( 1.0 - exp(-rn_rdt / ztau ) )
            hml=hbl-dh
            inhml = MAX( INT( dh / e3t_n(ibld) ) , 1 )
            imld = ibld - inhml
            zhml = fs_depw(imld)
            zdh = zhbl-zhml
         ELSE
!
! stable but collapsing.
!
            ztau = hbl / zwstrl
            zddhdt = 0.2 * hbl
            dh = dh * exp(-rn_rdt / ztau ) + zddhdt * ( 1.0 - exp(-rn_rdt / ztau ) )
            IF( dh >= hbl ) dh = zddhdt ! problem with rapid collapse
            inhml = MAX( INT( dh / e3t_n(ibld) ) , 1 )
            imld = ibld - inhml
            zhml = fs_depw(imld)
            zdh = zhbl-zhml
         ENDIF
      ENDIF

      END SUBROUTINE pyc_depth
!EOC
!BOP
      SUBROUTINE pyc_profiles( zdtdz_pyc, zdsdz_pyc, zdbdz_pyc, &
                            &  zdudz_pyc, zdvdz_pyc )
!
      IMPLICIT NONE
      REAL(wp), DIMENSION(:) :: zdtdz_pyc
      REAL(wp), DIMENSION(:) :: zdsdz_pyc
      REAL(wp), DIMENSION(:) :: zdbdz_pyc
      REAL(wp), DIMENSION(:) :: zdudz_pyc
      REAL(wp), DIMENSION(:) :: zdvdz_pyc
!
      INTEGER :: jk
      REAL(wp) :: ztgrad, zsgrad, zbgrad, zugrad, zvgrad
      REAL(wp) :: zgamma_b_nd, zgamma_b, znd
      REAL(wp) :: zzeta_m = 0.3



     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     !  Pycnocline gradients for scalars and velocity
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF ( lconv ) THEN
          ! Unstable conditions
!
! 25/05/17  Check that stable pycnocline is present
! change to NEMO code 08/06/17
         IF ( zdb_bl > 0._wp .and. zdbdz_ext >0._wp) THEN
            ztgrad = 0.5 * ( zdt_ml / ( zhbl - zhml ) ) + zdtdz_ext
            zsgrad = 0.5 * ( zds_ml / ( zhbl - zhml ) ) + zdsdz_ext
            zbgrad = 0.5 * ( zdb_ml / ( zhbl - zhml ) ) + zdbdz_ext
            zgamma_b_nd = zdbdz_ext * ( zhbl - zhml) / zdb_ml
            zgamma_b = ( 3.14159 / 4.0 ) * ( ( 0.5 + zgamma_b_nd ) / &
                         ( 1.0 - 0.25 * SQRT(3.14159 / 6.0 ) - 2.0 * zgamma_b_nd * zzeta_m ) )**2
            DO jk = 2 , ibld+ibld_ext
               znd = -( fs_depw(jk) - zhbl ) / ( zhbl - zhml )
               IF ( znd <= zzeta_m ) THEN
                  zdtdz_pyc(jk) =  zdtdz_ext + 0.5 * zdt_ml / ( zhbl - zhml ) *  EXP( -6.0 * ( znd - zzeta_m )**2 )
                  zdbdz_pyc(jk) =  zdbdz_ext + 0.5 * zdb_ml / ( zhbl - zhml ) * EXP( -6.0 * ( znd - zzeta_m )**2 )
                  zdsdz_pyc(jk) =  zdsdz_ext + 0.5 * zds_ml / ( zhbl - zhml ) * EXP( -6.0 * ( znd - zzeta_m )**2 )
               ELSE
                  zdtdz_pyc(jk) =  ztgrad * EXP( -zgamma_b * ( znd - zzeta_m )**2 )
                  zdbdz_pyc(jk) =  zbgrad * EXP( -zgamma_b * ( znd - zzeta_m )**2 )
                  zdsdz_pyc(jk) =  zsgrad * EXP( -zgamma_b * ( znd - zzeta_m )**2 )
               ENDIF
            END DO
         ELSE

         ENDIF
      ELSE
             ! stable conditions
         IF (zdhdt < 0.0) THEN
         ELSE
            IF( zdb_bl > 0 ) THEN
               IF (zhol >= 0.5) then
                  ztgrad = zdt_bl / zhbl
                  zsgrad = zds_bl / zhbl
                  zbgrad = zdb_bl / zhbl
                  DO jk = 2, ibld
                     znd = fs_depw(jk) / zhbl
                     zdtdz_pyc(jk) =  1.0 * ztgrad * EXP( -15.0 * ( znd - 0.9 )**2 )
                     zdbdz_pyc(jk) =  1.0 * zbgrad * EXP( -15.0 * ( znd - 0.9 )**2 )
                     zdsdz_pyc(jk) =  1.0 * zsgrad * EXP( -15.0 * ( znd - 0.9 )**2 )
                  END DO
               ELSE
                  ztgrad = ( zdt_ml / ( zhbl - zhml ) )
                  zsgrad = ( zds_ml / ( zhbl - zhml ) )
                  zbgrad = ( zdb_ml / ( zhbl - zhml ) )
                  DO jk = 2 , ibld
                     znd = -( fs_depw(jk) - zhml ) / ( zhbl - zhml )
                     zdtdz_pyc(jk) =  ztgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                     zdbdz_pyc(jk) =  zbgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                     zdsdz_pyc(jk) =  zsgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                  END DO
               ENDIF
            ELSE
            ENDIF
         ENDIF
      ENDIF
!
        !
      IF ( lconv ) THEN
          ! Unstable conditions
         zugrad = 0.7 * zdu_ml / ( zhbl - zhml ) + 0.3 * zustar**2 / &
                       ( ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * MIN( zla**(8./3.), 0.12_wp)  * zhml )
         zvgrad = ( 0.7 * zdv_ml + 2.0 * ff * zustke * delta / &
                        ( ( zwstrl**3 + 0.5 * zwstrc**3 )**(1./3.) ) )/ ( zhbl - zhml )
         zzeta_m=0.45
         DO jk = 2 ,  ibld + ibld_ext
            znd = -( fs_depw(jk) - zhbl ) / ( zhbl - zhml ) - zzeta_m
            IF ( znd <= 0.0 ) THEN
               zdudz_pyc(jk) = 1.25 * zugrad * EXP( 3.0 * znd )
               zdvdz_pyc(jk) = 1.25 * zvgrad * EXP( 3.0 * znd )
            ELSE
               zdudz_pyc(jk) = 1.25 * zugrad * EXP( -2.0 * znd )
               zdvdz_pyc(jk) = 1.25 * zvgrad * EXP( -2.0 * znd )
            ENDIF
         END DO
      ELSE
             ! stable conditions
         zugrad = 3.25 * zdu_bl / zhbl
         zvgrad = 2.75 * zdv_bl / zhbl
         DO jk = 2, ibld
            znd = fs_depw(jk) / zhbl
            IF ( znd <= 1.0 ) THEN
               zdudz_pyc(jk) = zugrad * EXP( -40.0 * ( znd - 1.0 )**2 )
            ELSE
               zdudz_pyc(jk) = zugrad * EXP( -20.0 * ( znd - 1.0 )**2 )
            ENDIF
            IF ( znd <= 0.85 ) THEN
               zdvdz_pyc(jk) = zvgrad * EXP( -20.0 * ( znd - 0.85 )**2 )
            ELSE
               zdvdz_pyc(jk) = zvgrad * EXP( -20.0 * ( znd - 0.85 )**2 )
            ENDIF
         END DO
      ENDIF

      END SUBROUTINE pyc_profiles
 END SUBROUTINE zdf_osm

!EOC
! !IROUTINE: Clean up the kpp module
!
! !INTERFACE:
   subroutine clean_osm()
!
! !DESCRIPTION:
!  De-allocate all memory allocated in init\_kpp().
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_osm'

   LEVEL2 'de-allocating OSMOSIS memory ...'
   if (allocated(nuh)) deallocate(nuh)
   if (allocated(nus)) deallocate(nus)
   if (allocated(num)) deallocate(num)
   if (allocated(nucl)) deallocate(nucl)
   if (allocated(gamh)) deallocate(gamh)
   if (allocated(gams)) deallocate(gams)
   if (allocated(gamu)) deallocate(gamu)
   if (allocated(gamv)) deallocate(gamv)
   if (allocated(Rig)) deallocate(Rig)
   if (allocated(tke)) deallocate(tke)
   if (allocated(zavt)) deallocate(zavt)
   if (allocated(zgamu)) deallocate(zgamu)
   if (allocated(zgamv)) deallocate(zgamv)
   if (allocated(zghamt)) deallocate(zghamt)
   if (allocated(zghams)) deallocate(zghams)
   ! if (allocated(zghamt_f)) deallocate(zghamt_f)
   ! if (allocated(zghams_f)) deallocate(zghams_f)
   if (allocated(z_w)) deallocate(z_w)
   if (allocated(z_r)) deallocate(z_r)
   if (allocated(h_r)) deallocate(h_r)
   if (allocated(fs_depw)) deallocate(fs_depw)
   if (allocated(fs_dept)) deallocate(fs_dept)
   if (allocated(e3w_n)) deallocate(e3w_n)
   if (allocated(e3t_n)) deallocate(e3t_n)
   if (allocated(ub)) deallocate(ub)
   if (allocated(vb)) deallocate(vb)
   if (allocated(wb)) deallocate(wb)
   if (allocated(tsn)) deallocate(tsn)
   if (allocated(rhd)) deallocate(rhd)
   if (allocated(rab_n)) deallocate(rab_n)
   if (allocated(rhop)) deallocate(rhop)
   if (allocated(zavtb)) deallocate(zavtb)
   if (allocated(zavmb)) deallocate(zavmb)

# ifdef EXTRA_OUTPUT

   if (allocated(turb1)) deallocate(turb1)
   if (allocated(turb2)) deallocate(turb2)
   if (allocated(turb3)) deallocate(turb3)
   if (allocated(turb4)) deallocate(turb4)
   if (allocated(turb5)) deallocate(turb5)

# endif

   return
 end subroutine clean_osm

!EOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute interior fluxes \label{sec:kppInterior}
!
! !INTERFACE:
   subroutine zdfosm_interior(nlev)
!
    !!  This is a modified version of the interior routine from the KPP model. Only the
    !!  diffusivity

! !DESCRIPTION:
! Here, the interior diffusivities (defined as the diffusivities outside the
! surface and bottom boundary layers) are computed. The algorithms are identical
! to those suggested by \cite{Largeetal94}. For numerical efficiency, the
! algorithms for different physical processes are active only if certain
! pre-processor macros are defined in {\tt cppdefs.h}.
! \begin{itemize}
!  \item The shear instability algorithm is active if the macro
!        {\tt KPP\_SHEAR} is defined.
!  \item The internal wave algorithm is active if the macro
!        {\tt KPP\_INTERNAL\_WAVE} is defined.
!  \item The convective instability algorithm is active if the macro
!        {\tt KPP\_CONVEC} is defined.
!  \item The double-diffusion algorithm is active if the macro
!        {\tt KPP\_DDMIX} is defined. Note that in this case, the
!        macro {\tt SALINITY} has to be defined as well.
! \end{itemize}

!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of grid cells
   integer                                       :: nlev


! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE , parameter       :: eps=1.0E-14

   integer                    :: i
   REALTYPE                   :: cff,shear2
   REALTYPE                   :: nu_sx,nu_sxc
   REALTYPE                   :: iwm,iws
   REALTYPE                   :: drhoT,drhoS,Rrho,nu_dds,nu_ddt

!
!-----------------------------------------------------------------------
!BOC
!
!-----------------------------------------------------------------------
! Compute gradient Richardson number
!-----------------------------------------------------------------------
!
   do i=1,nlev-1
      Rig(i) = NN(i) / (SS(i) + eps)
   enddo
   Rig(0)    = 0._wp
   Rig(nlev) = 0._wp
!
!-----------------------------------------------------------------------
!  Compute "interior" viscosities and diffusivities everywhere as
!  the superposition of three processes: local Richardson number
!  instability due to resolved vertical shear, internal wave
!  breaking, and double diffusion.
!-----------------------------------------------------------------------
!
   do i=1,nlev-1
!
!     Smooth gradient Richardson number
      Rig(i)=0.25*Rig(i-1) + 0.50*Rig(i) + 0.25*Rig(i+1)
!
!     Compute interior diffusivity due to shear instability mixing.
# ifdef KPP_SHEAR
      cff=min(1._wp,max(0._wp,Rig(i))/Ri0)
      nu_sx  = 1._wp-cff*cff
      nu_sx  = nu_sx*nu_sx*nu_sx
!
!     The shear mixing should be also a function of the actual magnitude
!     of the shear, see Polzin (1996, JPO, 1409-1425).
      shear2 = SS(i)
      cff    = shear2*shear2/(shear2*shear2+16.0E-10)
      nu_sx  = cff*nu_sx
# else
      nu_sx=0._wp
! KPP_SHEAR
# endif

#ifdef KPP_INTERNAL_WAVE
!
!      Compute interior diffusivity due to wave breaking
!
!      Version A, see Gargett and Holloway (1984)
!      cff  =  1._wp/sqrt(max(NN(i),1.0d-7))
!      iwm  =  1.0E-6*cff
!      iws  =  1.0E-7*cff

!     Version B, see Large et al. (1994)
      iwm  =  nuwm
      iws  =  nuws
#else
      iwm  =  0._wp
      iws  =  0._wp
! KPP_INTERNAL_WAVE
#endif


# ifdef KPP_CONVEC
!     Compute interior convective diffusivity due to static instability
!     mixing
      cff    =  max(NN(i),bvfcon)
      cff    =  min(1._wp,(bvfcon-cff)/bvfcon)
      nu_sxc =  1._wp-cff*cff
      nu_sxc =  nu_sxc*nu_sxc*nu_sxc
# else
      nu_sxc =  0._wp
! KPP_CONVEC
# endif
!
!     Sum contributions due to internal wave breaking, shear instability
!     and convective diffusivity due to shear instability.
      num(i)=iwm+nu0m*nu_sx+nu0c*nu_sxc
      nuh(i)=iws+nu0s*nu_sx+nu0c*nu_sxc
      nus(i)=nuh(i)
!
# ifdef KPP_DDMIX
!
!-----------------------------------------------------------------------
!  Compute double-diffusive mixing.  It can occur when vertical
!  gradient of density is stable but the vertical gradient of
!  salinity (salt figering) or temperature (diffusive convection)
!  is unstable.
!-----------------------------------------------------------------------
!
!     Compute double-diffusive density ratio, Rrho.
      drhoT = NNT(i)
      drhoS = NNS(i)
      Rrho= - drhoT/drhoS
      LEVEL2 'double diffusion'
!
!
!     Salt fingering case.
      if ((Rrho.gt.1._wp).and.(drhoS.gt.0._wp)) then
!
!        Compute interior diffusivity for double diffusive mixing of
!        salinity.  Upper bound "Rrho" by "Rrho0"; (Rrho0=1.9, nuf=0.001).
         Rrho=min(Rrho,Rrho0)
         nu_dds=1._wp-((Rrho-1._wp)/(Rrho0-1._wp))**2.0
         nu_dds=nuf*nu_dds*nu_dds*nu_dds
!
!        Compute interior diffusivity for double diffusive mixing
!        of temperature (fdd=0.7).
         nu_ddt=fdd*nu_dds
!
!
!     Diffusive convection case.
      elseif ((0._wp.lt.Rrho).and.(Rrho.lt.1._wp).and.(drhoS.lt.0._wp)) then
!
!        Compute interior diffusivity for double diffusive mixing of
!        temperature (Marmorino and Caldwell, 1976); (nu=1.5e-6,
!        tdd1=0.909, tdd2=4.6, tdd3=0.54).
         nu_ddt=nu*tdd1*exp(tdd2*exp(-tdd3*((1._wp/Rrho)-1._wp)))

!        Compute interior diffusivity for double diffusive mixing
!        of salinity (sdd1=0.15, sdd2=1.85, sdd3=0.85).
!
         if (Rrho.lt.0.5) then
            nu_dds=nu_ddt*sdd1*Rrho
         else
            nu_dds=nu_ddt*(sdd2*Rrho-sdd3)
         endif
      else
         nu_ddt=0._wp
         nu_dds=0._wp
      endif
!
!     Add double diffusion contribution to temperature and salinity
!     mixing coefficients.
      nuh(i)=nuh(i)  + nu_ddt
      nus(i)=nuh(i)  + nu_dds

! KPP_DDMIX
# endif

   enddo ! loop over interior points

 end subroutine zdfosm_interior
!EOC

END MODULE zdfosm

