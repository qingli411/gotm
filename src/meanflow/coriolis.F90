#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The Coriolis rotation \label{sec:coriolis}
!
! !INTERFACE:
   subroutine coriolis(nlev,dt)
!
! !DESCRIPTION:
!  This subroutine carries out the Coriolis rotation by applying a
!  $2\times 2$ rotation matrix with the angle $f\Delta t$ on the
!  horizontal velocity vector $(U,V)$.
!
! !USES:
   USE meanflow, only: u,v,cori
   USE meanflow, only: stokes_coriolis
   ! USE observations, only: ustokes,vstokes
   USE observations,   only: us_x,us_y,delta
   USE meanflow, only: z
!
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  Add Stokes-Coriolis term
!  Qing Li, 20180509
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: ua,omega,cosomega,sinomega
   REALTYPE                  :: depth
!
!-----------------------------------------------------------------------
!BOC

   omega=cori*dt
   cosomega=cos(omega)
   sinomega=sin(omega)

   if (stokes_coriolis) then
      do i=1,nlev
         ! ua=u(i)
         ! u(i)= u(i)*cosomega + (v(i)+vstokes(i))*sinomega
         ! v(i)=-(ua  +ustokes(i))*sinomega + v(i)*cosomega
         depth=(z(nlev)-z(i))
         ua=u(i)
         u(i)= u(i) *cosomega+(v(i)+us_y*exp(-depth/delta))*sinomega
         v(i)=-(ua+us_x*exp(-depth/delta))   *sinomega+v(i)*cosomega
      end do
   else
      do i=1,nlev
         ua=u(i)
         u(i)= u(i) *cosomega+v(i)*sinomega
         v(i)=-ua   *sinomega+v(i)*cosomega
      end do
   end if

   return
   end subroutine coriolis
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
