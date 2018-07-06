#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update turbulence production\label{sec:production}
!
! !INTERFACE:
!RRH:vvv
!  subroutine production(nlev,NN,SS,CSSTK,xP)
   subroutine production(nlev,NN,SS,SSSTK,CSSTK,xP)
!RRH:^^^
!
! !DESCRIPTION:
!  This subroutine calculates the production terms of turbulent kinetic
!  energy as defined in \eq{PandG} and the production of buoayancy
!  variance as defined in \eq{Pbvertical}.
!  The shear-production is computed according to
!  \begin{equation}
!    \label{computeP}
!     P = \nu_t (M^2 + \alpha_w N^2) + X_P
!    \comma
!  \end{equation}
!  with the turbulent diffusivity of momentum, $\nu_t$, defined in
!  \eq{nu}. The shear-frequency, $M$, is discretised as described
!  in \sect{sec:shear}.
!   The term multiplied by $\alpha_w$ traces back to
!  a parameterisation of breaking internal waves suggested by
!  \cite{Mellor89}. $X_P$ is an extra production term, connected for
!  example with turbulence production caused by sea-grass, see
!  \eq{sgProduction} in  \sect{sec:seagrass}. {\tt xP} is an {\tt optional}
!  argument in the FORTRAN code.
!
!  Similarly, according to \eq{PeVertical}, the buoyancy production
!  is computed from the expression
!  \begin{equation}
!   \label{computeG}
!    G=-\nu^B_t N^2 + \tilde{\Gamma}_B
!    \comma
!  \end{equation}
!  with the turbulent diffusivity, $\nu^B_t$, defined in
!  \eq{nu}. The second term in \eq{computeG} represents the non-local
!  buoyancy flux. The buoyancy-frequency, $N$, is discretised as described
!  in \sect{sec:stratification}.
!
!  The production of buoyancy variance by vertical meanflow gradients follows
!  from \eq{PeVertical} and \eq{computeG}
!  \begin{equation}
!   \label{computePb}
!    P_b = -G N^2
!    \point
!  \end{equation}
!  Thus, according to the definition of the potential energy \eq{defkb},
!  the buoyancy production $G$ describes the conversion between turbulent
!  kinetic and potential energy in \eq{tkeA} and \eq{kbeq}, respectively.
!
! !USES:
   use turbulence, only: P,B,Pb
   use turbulence, only: PS
   use turbulence, only: num,nuh
!RRH:vvv
   use turbulence, only: nucl
!RRH:^^^
   use turbulence, only: alpha,iw_model
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  boyancy frequency squared (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev)

!  shear-frequency squared (1/s^2)
   REALTYPE, intent(in)                :: SS(0:nlev)

!RRH:vvv
!  Stokes shear squared (1/s^2)
   REALTYPE, intent(in)                :: SSSTK(0:nlev)

!RRH:^^^
!  Stokes-Eulerian cross-shear (1/s^2)
   REALTYPE, intent(in)                :: CSSTK(0:nlev)

!  TKE production due to seagrass
!  friction (m^2/s^3)
   REALTYPE, intent(in), optional      :: xP(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!     Qing Li, 20180418, update Stokes production, PS
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   REALTYPE                      :: alpha_eff
   integer                       :: i
!-----------------------------------------------------------------------
!BOC
   alpha_eff=_ZERO_
   if (iw_model.eq.1) then
      alpha_eff=alpha
   end if

   if ( PRESENT(xP) ) then
      do i=0,nlev
         P(i)    =  num(i)*( SS(i)+alpha_eff*NN(i) ) + xP(i)
!RRH:vvv
#if defined(STOKESFLUX)
         P(i)    =  P(i) + nucl(i)*CSSTK(i)
#endif
!RRH:^^^
         B(i)    = -nuh(i)*NN(i)
         Pb(i)   = -  B(i)*NN(i)
         PS(i)   =  num(i)*CSSTK(i)
!RRH:vvv
#if defined(STOKESFLUX)
         PS(i)    =  PS(i) + nucl(i)*SSSTK(i)
#endif
!RRH:^^^
      enddo
   else
      do i=0,nlev
         P(i)    =  num(i)*( SS(i)+alpha_eff*NN(i) )
!RRH:vvv
#if defined(STOKESFLUX)
         P(i)    =  P(i) + nucl(i)*CSSTK(i)
#endif
!RRH:^^^
         B(i)    = -nuh(i)*NN(i)
         Pb(i)   = -  B(i)*NN(i)
         PS(i)   =  num(i)*CSSTK(i)
!RRH:vvv
#if defined(STOKESFLUX)
         PS(i)    =  PS(i) + nucl(i)*SSSTK(i)
#endif
!RRH:^^^
      enddo
   endif

   return
   end subroutine production
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
