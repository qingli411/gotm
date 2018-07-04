#include"cppdefs.h"
!/ ------------------------------------------------------------------- /
   module langmuir
!/  ------------------------------------------------------------------- /
!/ Contains a routine for computing Langmuir number moved from Stokes
!   drift module in ../observations/stokes.F90
!   Moved here from that routine to fix a compilation bug in June 2018
!   - This could be placed inside kpp, or possibly we can move other
!     Langmuir number bits of code from kpp here.
!   - The reliance on all time/space Stokes drift data from Stokes for
!     this routine is why I have given it its own module for now,
!     perhaps prior to a code release this should be improved upon to
!     make the hurricane cases faster.
!/
!/   List of subroutines:
!/    1. Get_LaNum
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
! Need access to all the Stokes drift data.  Much of this lacks an
!  approach w/ efficiency in mind, but it works for now.
   use Stokes,       only: &
        JDATESPC,&! JDATESPC - Julian Date for spectra
        SPC,&! SPC - Spectra (f'n of wavenum, dir, and time)
        STOKESx,STOKESy,&! STOKESx(/y) - Stokes drift in x(/y)
        !direction (f'n of depth, time)
        WD,&!WD - Wave directions
        WN,&!WN - Wavenumbers
        Z_CEN,&! Z_CEN - Depth of center
        NSPC, NDIR, NK, NZ
   ! NSPC - number spectra (total times)
   ! NDIR - number of spectral directions
   ! NK   - numer of spectral wavenumbers
   ! NZ   - number of levels

   use airsea,       only: u10, v10
   use meanflow,     only: pi, gravity, z, zi
   use turbulence,   only: kappa
   use observations, only: ustokes, vstokes, us_x, us_y, Hs
   use observations, only: nfreq, wav_freq, wav_spec, wav_xcmp, wav_ycmp

#ifdef KPP_CVMIX
!  use CVMix
   use cvmix_kinds_and_types, only: cvmix_global_params_type
   use cvmix_kpp,             only: cvmix_kpp_ustokes_SL_model
   use cvmix_put_get,         only: cvmix_put
#endif

   private

!  turbulent Langmuir number
   REALTYPE, public                      :: La_Turb
!  surface layer averaged Langmuir number
   REALTYPE, public                      :: La_SL
!  surface layer averaged, projected Langmuir number (Van Roekel et al., 2012)
   REALTYPE, public                      :: La_SLP1
!  surface layer averaged, projected Langmuir number (Reichl et al., 2016)
   REALTYPE, public                      :: La_SLP2
!  angles between wind and waves (radian)
   REALTYPE, public                      :: theta_WW
!  angles between wind and Langmuir cells (radian)
   REALTYPE, public                      :: theta_WL
!  method to compute Langmuir number
   integer, public                       :: langmuir_number_method
!  name of the file which contains the Langmuir number
   character(len=PATH_MAX), public       :: langmuir_number_file

!  method of enhancement factor
   integer, parameter, public ::  LANGMUIR_NUMBER_OPT_NONE    = 0
   integer, parameter, public ::  LANGMUIR_NUMBER_OPT_MODEL   = 1
   integer, parameter, public ::  LANGMUIR_NUMBER_OPT_READ    = 2
   integer, parameter, public ::  LANGMUIR_NUMBER_OPT_SPEC    = 3
   integer, parameter, public ::  LANGMUIR_NUMBER_OPT_USTOKES = 4
   integer, parameter, public ::  LANGMUIR_NUMBER_OPT_HURR    = 5

   public init_langmuir, langmuir_number

!-----------------------------------------------------------------------
   contains
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize Langmuir number
!
! !INTERFACE:
   subroutine init_langmuir(namlst,fn)
! !DESCRIPTION:
!  This routine reads the namelist {\tt langmuir} and initialize
!  the Langmuir number.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  namelist reference
   integer,          intent(in)        :: namlst

!  filename containing namelist
   character(len=*), intent(in)        :: fn

!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   namelist /langmuir/               langmuir_number_method,            &
                                     langmuir_number_file

!
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------

   LEVEL1 'init_langmuir'

   ! default values
   langmuir_number_method = LANGMUIR_NUMBER_OPT_NONE
   langmuir_number_file = ''
   La_Turb  = _ONE_/SMALL
   La_SL    = _ONE_/SMALL
   La_SLP1  = _ONE_/SMALL
   La_SLP2  = _ONE_/SMALL
   theta_WW = _ZERO_
   theta_WL = _ZERO_

   ! read the variables from the namelist file
   open(namlst,file=fn,status='old',action='read',err=80)

   LEVEL2 'reading langmuir namelist...'

   read(namlst,nml=langmuir,err=81)
   close (namlst)

   LEVEL2 'done.'

   select case(langmuir_number_method)
   case(LANGMUIR_NUMBER_OPT_NONE)
      LEVEL3 'No Langmuir number is computed'
   case(LANGMUIR_NUMBER_OPT_MODEL)
      LEVEL3 'Approximate Langmuir number from simple model (Li et al., 2017)'
   case(LANGMUIR_NUMBER_OPT_READ)
      LEVEL3 'Read Langmuir number from file: ', trim(langmuir_number_file)
   case(LANGMUIR_NUMBER_OPT_SPEC)
      LEVEL3 'Compute Langmuir number from wave spectrum'
   case(LANGMUIR_NUMBER_OPT_USTOKES)
      LEVEL3 'Compute Langmuir number from Stokes drift'
   case (LANGMUIR_NUMBER_OPT_HURR)
      LEVEL3 'Compute Langmuir number from hurricane spectrum'
   case default
      stop 'init_langmuir: unsupported langmuir_number_method'
   end select

   return

80 FATAL 'I could not open "langmuir.nml"'
   stop 'init_langmuir'
81 FATAL 'I could not read "langmuir" namelist'
   stop 'init_langmuir'

   end subroutine init_langmuir
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute Langmuir number
!
! !INTERFACE:
   subroutine langmuir_number(nlev, u_taus, hbl)
! !DESCRIPTION:
!  This routine computes Langmuir number in various definition and
!  returns the corresponding one according to the Langmuir turbulence
!  parameterization option.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: u_taus, hbl
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   REALTYPE                            :: wind10m, ussl_model

#ifdef KPP_CVMIX
!  CVMix datatypes
   type(cvmix_global_params_type)      :: CVmix_params
#endif
!
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------

   select case (langmuir_number_method)
   case (LANGMUIR_NUMBER_OPT_NONE)
      La_Turb  = _ONE_/SMALL
      La_SL    = _ONE_/SMALL
      La_SLP1  = _ONE_/SMALL
      La_SLP2  = _ONE_/SMALL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
   case (LANGMUIR_NUMBER_OPT_MODEL)
#ifdef KPP_CVMIX
      ! 10-meter wind speed
      wind10m = sqrt(u10**2+v10**2)
      ! assume surface Stokes drift is proportional to U_10
      La_Turb  = sqrt(u_taus/wind10m/0.0162)
      ! this returns the Li et al. 2017 SL averaged Langmuir number
      call cvmix_put(CVmix_params, 'Gravity', gravity)
      ussl_model = cvmix_kpp_ustokes_SL_model(wind10m, hbl, CVmix_params)
      La_SL = sqrt(u_taus/ussl_model)
      La_SLP1  = La_SL
      La_SLP2  = La_SL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
#else
      stop "langmuir_number: compile with CVMix to use the efactor model"
#endif
   case (LANGMUIR_NUMBER_OPT_READ)
      ! TODO: read from file <13-12-17, Qing Li> !
      La_Turb  = _ONE_/SMALL
      La_SL    = _ONE_/SMALL
      La_SLP1  = _ONE_/SMALL
      La_SLP2  = _ONE_/SMALL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
   case (LANGMUIR_NUMBER_OPT_SPEC)
      call langmuir_number_spec(wav_freq, wav_spec, wav_xcmp, wav_ycmp, &
                            u10, v10, u_taus, hbl)
   case (LANGMUIR_NUMBER_OPT_USTOKES)
      call langmuir_number_ustokes(nlev, zi, ustokes, vstokes, &
                               us_x, us_y, u10, v10, Hs, u_taus, hbl)
   case (LANGMUIR_NUMBER_OPT_HURR)
      call langmuir_number_hurr(hbl*0.2,u_taus)
   case default
      stop 'langmuir_number: unsupported efactor_method'
   end select

   end subroutine langmuir_number
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute the Langmuir number from Stokes drift profile
!
! !INTERFACE:
   subroutine langmuir_number_ustokes(nlev,z_w, &
                        uus,vus,uus0,vus0,uwnd,vwnd,hsw,ustar,hbl)
!
! !DESCRIPTION:
!  This routine computes different definitions of Langmuir number
!  from Stokes drift profile.
!
! TODO: More detailed description <09-04-18, Qing Li> !
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   ! number of grid
   integer, intent(in)                 :: nlev
   ! depth at the grid interface
   REALTYPE, intent(in)                :: z_w(0:nlev)
   ! Stokes drift profile (m/s)
   REALTYPE, intent(in)                :: uus(0:nlev), vus(0:nlev)
   ! surface Stokes drift (m/s)
   REALTYPE, intent(in)                :: uus0, vus0
   ! surface wind at 10 meter (m/s)
   REALTYPE, intent(in)                :: uwnd, vwnd
   ! significant wave height (m)
   REALTYPE, intent(in)                :: hsw
   ! friction velocity (m/s)
   REALTYPE, intent(in)                :: ustar
   ! boundary layer depth (m)
   REALTYPE, intent(in)                :: hbl
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   integer                             :: k, kk, kref
   REALTYPE                            :: ussl, vssl, us0
   REALTYPE                            :: hsl, dz
!
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------

!  magnitude of surface Stokes drift
   us0 = sqrt(uus0**2.+vus0**2.)

!  surface layer
   hsl = 0.2*abs(hbl)

!  determine which layer contains surface layer
   do kk = nlev,k,-1
      if (z_w(nlev)-z_w(kk-1) .ge. hsl) then
         kref = kk
         exit
      end if
   end do

!  calculate the surface layer averaged Stokes drift
   if (kref < nlev) then
      ussl =   uus(kref)*(hsl+z_w(kref))
      vssl =   vus(kref)*(hsl+z_w(kref))
      do kk = nlev,kref+1,-1
         dz = z_w(kk)-z_w(kk-1)
         ussl = ussl + uus(kk)*dz
         vssl = vssl + vus(kk)*dz
      end do
      ussl = ussl/hsl
      vssl = vssl/hsl
   else
      ussl = uus(nlev)
      vssl = vus(nlev)
   end if

   if (us0 .gt. _ZERO_) then
      La_Turb = sqrt(ustar/us0)
      La_SL = sqrt(ustar/sqrt(ussl**2.+vssl**2.))
      ! angles between wind and waves
      theta_WW = atan2(vssl,ussl)-atan2(vwnd,uwnd)
      ! angles between wind and LCs
      ! (approximate from law of the wall, Van Roekel et al., 2012)
      theta_WL = atan(sin(theta_WW) &
            /(ustar/us0/kappa*log(max(abs(hbl/4./hsw),_ONE_))+cos(theta_WW)))
      La_SLP1 = La_SL*sqrt(abs(cos(theta_WL))/abs(cos(theta_WW-theta_WL)))
      La_SLP2 = La_SL*sqrt(_ONE_/max(abs(cos(theta_WW-theta_WL)),SMALL))
   else
      La_Turb = _ONE_/SMALL
      La_SL = _ONE_/SMALL
      La_SLP1 = _ONE_/SMALL
      La_SLP2 = _ONE_/SMALL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
   end if

   end subroutine langmuir_number_ustokes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute the Langmuir number from wave spectrum
!
! !INTERFACE:
   subroutine langmuir_number_spec(freq,spec,xcmp,ycmp,uwnd,vwnd,ustar,hbl)
!
! !DESCRIPTION:
!  This routine computes different definitions of Langmuir number
!  from wave spectrum.
!
! TODO: More detailed description  <22-12-17, Qing Li> !
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: spec(:), xcmp(:), ycmp(:)
   REALTYPE, intent(in)                :: freq(:)
   REALTYPE, intent(in)                :: uwnd, vwnd, ustar, hbl
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   integer                             :: i
   REALTYPE                            :: ussl, vssl, us0
   REALTYPE                            :: hsw, hsl
   REALTYPE                            :: tmp, factor, factor2
   REALTYPE                            :: dfreqc, freqc
!
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------
!  initialization
   ussl = _ZERO_
   vssl = _ZERO_
   us0  = _ZERO_
   hsw  = _ZERO_
   tmp  = _ZERO_
!  surface layer
   hsl = 0.2*abs(hbl)
!  calculate the significant wave height, surface Stokes drift and
!  surface layer averaged Stokes drift
   factor = 2.*pi/hsl
   do i=1,nfreq
      factor2 = 8.*pi**2.*freq(i)**2./gravity
      tmp  = tmp + spec(i)
      us0 = us0 + 2.*pi*freq(i)*factor2*spec(i)
      ussl = ussl + factor*freq(i)*spec(i)*xcmp(i)*(_ONE_-exp(-factor2*hsl))
      vssl = vssl + factor*freq(i)*spec(i)*ycmp(i)*(_ONE_-exp(-factor2*hsl))
   end do
   hsw = 4.*sqrt(tmp)
!  cutoff frequency
   freqc = 1.5*freq(nfreq)-0.5*freq(nfreq-1)
   dfreqc = freq(nfreq)-freq(nfreq-1)
!  add contribution from a f^-5 tails
   factor = 16.*pi**3.*freqc**4./gravity
   us0 = us0 + factor*spec(nfreq)/dfreqc
   factor = 4.*pi**2.*freqc**2./3./hsl &
       *(_ONE_-(_ONE_-16.*pi**2.*freqc**2.*hsl/gravity) &
       *exp(-8.*pi**2.*freqc**2.*hsl/gravity))
   ussl = ussl + factor*xcmp(nfreq)*spec(nfreq)/dfreqc
   vssl = vssl + factor*ycmp(nfreq)*spec(nfreq)/dfreqc

   if (us0 .gt. _ZERO_) then
      La_Turb = sqrt(ustar/us0)
      La_SL = sqrt(ustar/sqrt(ussl**2.+vssl**2.))
      ! angles between wind and waves
      theta_WW = atan2(vssl,ussl)-atan2(vwnd,uwnd)
      ! angles between wind and LCs
      ! (approximate from law of the wall, Van Roekel et al., 2012)
      theta_WL = atan(sin(theta_WW) &
            /(ustar/us0/kappa*log(max(abs(hbl/4./hsw),_ONE_))+cos(theta_WW)))
      La_SLP1 = La_SL*sqrt(abs(cos(theta_WL))/abs(cos(theta_WW-theta_WL)))
      La_SLP2 = La_SL*sqrt(_ONE_/max(abs(cos(theta_WW-theta_WL)),SMALL))
   else
      La_Turb = _ONE_/SMALL
      La_SL = _ONE_/SMALL
      La_SLP1 = _ONE_/SMALL
      La_SLP2 = _ONE_/SMALL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
   end if

   end subroutine langmuir_number_spec
!EOC

   subroutine langmuir_number_hurr(AVGDPTH,USTin)
!/ ------------------------------------------------------------------- /
!/ Name    : Get_LaNum
!/ Purpose : Recieves u* and calculates Langmuir number.  Langmuir
!/           number uses averaged Stokes drift over SL and accounts
!/           for misalignment between shear and Stokes drift.  SL
!/           integrated Stokes drift computed from spectra using SL
!/           depth.
!/ Author  : Brandon Reichl (URI-GSO)
!/ Intent(in) : AVGDPTH - SL depth for Stokes drift averaging.
!/ Intent(in) : USTin - u* magnitude
!/ Intent(out): LAout - Langmuir number output
!/ ------------------------------------------------------------------- /
   use time,     only: julianday,secondsofday
   use meanflow, only: u, v
!/
   IMPLICIT NONE
   REALTYPE, INTENT(IN)  :: AVGDPTH,USTin
   REAL, PARAMETER :: GRAV=9.806
   REALTYPE              :: JDAYIn
   REALTYPE              :: DT,DZ
   INTEGER               :: K,T,I, ZBL,Z
   REALTYPE              :: USX, USY
   REALTYPE :: DKDTH(NK), USU(NZ), VSU(NZ)
   REALTYPE :: SPC2(NK,NDIR)
   REALTYPE :: la_shear_x, la_shear_y, la_shear_d
   REALTYPE :: STK_d, LA_LOWER, LAdir_L
   !/
   JDAYIn=real(secondsofday/86400.)
   JDAYIN=JDAYIN+REAL(JULIANDAY)
   IF (JDayIn.GT.JDATESPC(1).AND.JDayIn.LT.JDATESPC(Nspc)) THEN
      I=1
      DO WHILE(JDayIN.GT.JDATESPC(I+1))
         I=I+1
      ENDDO
      DT=JDATESPC(I+1)-JDATESPC(I)
      DO K=1,NK
         DO T=1,NDIR
            SPC2(K,T) =  SPC(I+1,K,T)*(JDAYIN-JDATESPC(I))/DT &
                 + SPC(I,K,T)*  (JDATESPC(I+1)-JDAYIN)/DT
         ENDDO
      ENDDO
      USX=_ZERO_
      USY=_ZERO_
      DO K=2,NK-1
         DKDTH(K) = 0.5*(WN(K+1)-WN(k-1))*(WD(2)-WD(1))
      ENDDO
      DKDTH(1)= (WN(2)-WN(1)) *( WD(2) - WD(1))
      DKDTH(NK)=(WN(nk)-WN(nk-1))*(WD(2) - WD(1))
      DO K=1,NK
         IF (2.*3.1415/WN(k).gt.1.0) then
            DO T=1,NDIR
               USX = USX + SPC2(K,T) * SQRT(GRAV*WN(K))*COS(WD(T)) &
                    * (1.0 - EXP(2.*WN(K)*(0.-AVGDPTH))) * DKDTH(K)
               USY = USY + SPC2(K,T) * SQRT(GRAV*WN(K))*SIN(WD(T)) &
                    * (1.0 - EXP(2.*WN(K)*(0.-AVGDPTH))) * DKDTH(K)
            ENDDO
         ENDIF
      ENDDO
      USX = USX / AVGDPTH
      USY = USY / AVGDPTH
      ! find depth of 20% of ML
      ZBL = NZ-1
      IF (AvgDpth.lt.Z_CEN(ZBL)) THEN
         DO WHILE(AvgDpth.lt.Z_CEN(ZBL-1))
            ZBL = ZBL - 1
         ENDDO
      ENDIF
      DO Z = NZ,1, -1
         USU(Z)=Stokesx(I+1,Z)*(JDAYIN-JDATESPC(I))/DT+Stokesx(I,Z)* &
              (JDATESPC(I+1)-JDAYIN)/DT
         VSU(Z)=Stokesy(I+1,Z)*(JDAYIN-JDATESPC(I))/DT+Stokesy(I,Z)* &
              (JDATESPC(I+1)-JDAYIN)/DT
      ENDDO

      LA_Shear_x = (USU(NZ)+U(NZ)-USU(NZ-1)-U(NZ-1)) / (Z_CEN(NZ)-Z_CEN(ZBL))
      LA_Shear_y = (VSU(NZ)+V(NZ)-VSU(NZ-1)-V(NZ-1)) / (Z_CEN(NZ)-Z_CEN(ZBL))
      DO Z = NZ, max(2,(min(ZBL+1,NZ-1))),-1
         LA_shear_x = LA_shear_x + &
              (USU(Z)+U(Z)-USU(Z-1)-U(Z-1)) / (Z_CEN(NZ)-Z_CEN(ZBL))
         LA_shear_y = LA_shear_y + &
              (VSU(Z)+V(Z)-VSU(Z-1)-V(Z-1)) / (Z_CEN(NZ)-Z_CEN(ZBL))
      enddo
      LA_SHEAR_D = ATAN2(LA_SHEAR_Y,LA_SHEAR_X)
      STK_d=atan2(USY,USX)
      LA_LOWER=sqrt(USX**2+USY**2)
      LAdir_L = cos(STK_d-LA_shear_d )
      LAdir_L = max(LAdir_L,1.E-8)
      La_SLP2 = sqrt (USTin / LA_LOWER) * min(10.0,(sqrt(_ONE_/LAdir_L)))
      La_Turb = sqrt(USTin/sqrt(USU(NZ)**2.+VSU(NZ)**2.))
      La_SL = sqrt(USTin/LA_LOWER)
      La_SLP1 = La_SLP2
      theta_WW =  STK_d-LA_shear_d
      theta_WL = _ZERO_
   ELSE
      La_Turb = _ONE_/SMALL
      La_SL = _ONE_/SMALL
      La_SLP1 = _ONE_/SMALL
      La_SLP2 = _ONE_/SMALL
      theta_WW = _ZERO_
      theta_WL = _ZERO_
   ENDIF

   end subroutine langmuir_number_hurr

 end module langmuir
