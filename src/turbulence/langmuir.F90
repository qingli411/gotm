#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: langmuir --- Langmuir turbulence \label{sec:langmuir}
!
! !INTERFACE:
   module langmuir
!
! !DESCRIPTION:
!
! !USES:
   use airsea,       only: u10, v10
   use meanflow,     only: zi
   use turbulence,   only: kappa
   use observations, only: ustokes, vstokes, us_x, us_y, Hs

   IMPLICIT NONE

!  default: all is private.
   private
!
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
!
! !PUBLIC MEMBER FUNCTIONS:
   public langmuir_number

!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
!
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------

   contains

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
!
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------

   call langmuir_number_ustokes(nlev, zi, ustokes, vstokes, &
                                us_x, us_y, u10, v10, Hs, u_taus, hbl)

   end subroutine langmuir_number
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute the Langmuir number from Stokes drift profile
!
! !INTERFACE:
   subroutine langmuir_number_ustokes(nlev,z_w,uus,vus,uus0,vus0,uwnd,vwnd,hsw,ustar,hbl)
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

 end module langmuir
