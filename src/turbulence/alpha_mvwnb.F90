#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update dimensionless alpha's\label{sec:alpha}
!
! !INTERFACE:
   subroutine alpha_mvwnb(nlev,NN,SS,CSSTK,SSSTK)
!
! !DESCRIPTION:
! This subroutine updates the dimensionless numbers $\alpha_M$, $\alpha_N$,
! and $\alpha_b$ according to \eq{alphaMN}. Note that according to \eq{Nbar}
! and \eq{NbarVertical} the following identities are valid
! \begin{equation}
!  \label{alphaIdentities}
!    \alpha_M = \overline{S}^2 \comma
!    \alpha_V = dudz*dusdz + dvdz*dvsdz \comma
!    \alpha_W = dusdz*dusdz + dvsdz*dvsdz \comma
!    \alpha_N = \overline{N}^2 \comma
!    \alpha_b = \overline{T}   \point
! \end{equation}
!
!
! !USES:
  use turbulence,  only:     tke,eps,kb
  use turbulence,  only:     as,an,at
  use turbulence,  only:     av,aw
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer,  intent(in)      :: nlev
  REALTYPE, intent(in)      :: NN(0:nlev),SS(0:nlev)
  REALTYPE, intent(in)      :: CSSTK(0:nlev),SSSTK(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!     Qing Li, 20180419, code from Ramsey Harcout
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  integer              :: i
  REALTYPE             :: tau2

!-----------------------------------------------------------------------
!BOC

  do i=0,nlev
     tau2   = tke(i)*tke(i) / ( eps(i)*eps(i) )
     as(i)  = tau2 * SS(i)
     av(i)  = tau2 * CSSTK(i)
     aw(i)  = tau2 * SSSTK(i)
     an(i)  = tau2 * NN(i)
     at(i)  = tke(i)/eps(i) * kb(i)/eps(i)

!    clip negative values
     as(i) = max(as(i),1.e-10*_ONE_)
     aw(i) = max(aw(i),1.e-10*_ONE_)
     at(i) = max(at(i),1.e-10*_ONE_)
  end do

  return
end subroutine alpha_mvwnb

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
