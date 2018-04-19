#include"cppdefs.h"
!-----------------------------------------------------------------------
!
! !ROUTINE: The Langmuir turbulence quasi-equilibrium stability functions  after Harcourt (2015)\label{sec:cmueDH15}
!
! !INTERFACE:
   subroutine cmue_d_h15(nlev)
!
! !DESCRIPTION:
!
!  This subroutine updates the explicit solution of
!  \eq{bijVertical} and \eq{giVertical} under the same assumptions
!  as those discussed in \sect{sec:cmueC}. Now, however, an additional
!  equilibrium assumption is invoked. With the help of \eq{PeVertical},
!  one can write the equilibrium condition for the TKE as
! \begin{equation}
!  \label{quasiEquilibrium}
!     \dfrac{P+G}{\epsilon} =
!    \hat{c}_\mu(\alpha_M,\alpha_N) \alpha_M
!    - \hat{c}'_\mu(\alpha_M,\alpha_N) \alpha_N = 1
!   \comma
! \end{equation}
! where \eq{alphaIdentities} has been used. This is an implicit relation
! to determine $\alpha_M$ as a function of $\alpha_N$.
! With the definitions given in \sect{sec:cmueC}, it turns out that
! $\alpha_M(\alpha_N)$ is a quadratic polynomial that is easily solved.
! The resulting value for $\alpha_M$ is substituted into the stability
! functions described in \sect{sec:cmueC}. For negative $\alpha_N$
! (convection) the shear number $\alpha_M$ computed in this way may
! become negative. The value of $\alpha_N$ is limited such that this
! does not happen, see \cite{UmlaufBurchard2005a}.
! This Langmuir turbulence version includes the CL Vortex force in
! the algebraic models as well as in the vortex production of TKE and L or epsilon
!
! !USES:
   use turbulence, only: an,as,at
! nondimensional forcing functions for Eulerian shear dot Stokes shear, and Stokes shear squared:
! Also, surface proximity function SPF=(1-fzs), goes to zero at surface as tanh(0.25*z/l_S) where l_S
! the vortex-production-weighted dissipation length scale
   use turbulence, only: av, aw, SPF
   use turbulence, only: tke, L
   use turbulence, only: cmue1,cmue2
! Using 'nonlocal' fluxes to carry flux down the Stokes gradient.  
! This is not a nonlocal flux, so don't call it that outside of the code!
! Any additional nonlocal fluxes in future models will simply add to these, 
! which all get zeroed out in turbulence.F90 before calls to cmue_XX
  use turbulence,   only: gamu,gamv,gamh,gams

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer, intent(in)       :: nlev

! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: small       = 1.0D-10

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!  Converted to nefarious uses by Ramsey Harcourt
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
!
     integer                 ::   i
     REALTYPE            ::   Gv, Gs, Gh
     REALTYPE            ::   Sm, Ss, Sh

     REALTYPE, parameter :: my_A1 = 0.92D0
     REALTYPE, parameter :: my_A2 = 0.74D0
     REALTYPE, parameter :: my_B1 = 16.6D0
     REALTYPE, parameter :: my_B2 = 10.1D0
     REALTYPE, parameter :: my_C1 = 0.08D0
     REALTYPE, parameter :: my_C2 = 0.7D0
     REALTYPE, parameter :: my_C3 = 0.2D0
     REALTYPE, parameter :: h15_Ghmax = 0.029D0
     REALTYPE, parameter :: h15_Ghmin = -0.28D0
     REALTYPE, parameter :: h15_Gvmax = 0.024D0
     REALTYPE, parameter :: h15_Gsmax = 0.152D0

     REALTYPE :: h15_Shn0, h15_Shnh, h15_Shns, h15_Shnv
     REALTYPE :: h15_Shdah, h15_Shdav, h15_Shdbh
     REALTYPE :: h15_Shdv, h15_Shdvh, h15_Shdvv
     REALTYPE :: h15_Ssn0, h15_Ssdh, h15_Ssdv
     REALTYPE :: h15_Smn0, h15_SmnhSh h15_SmnsSs
     REALTYPE :: h15_Smdh h15_Smdv

     REALTYPE :: tmp0,tmp1,tmp2

!-----------------------------------------------------------------------
! These constants  above & below should all be computed elsewhere in advance, 
! subject to adjustments in A's, B's & C's. Just sticking them all in here for now.

     h15_Shn0=my_A2*(1-6.D0*my_A1/my_B1)
     h15_Shnh=-9.D0*my_A1*my_A2*( my_A2*(1-6.D0*my_A1/my_B1))
     h15_Shns=9.D0*my_A1*my_A2*(1-6.D0*my_A1/my_B1)*                                    &
              (2.D0*my_A1+my_A2)
     h15_Shnv=9.D0*my_A1*my_A2*                                                         &
              (my_A2*(1-6.D0*my_A1/my_B1-3.D0*my_C1)                                    &
               -2.D0*my_A1*(1-6.D0*my_A1/my_B1+3.D0*my_C1))
     h15_Shdah=-9.D0*my_A1*my_A2
     h15_Shdav=-36.D0*my_A1*my_A1
     h15_Shdbh=-3.D0*my_A2*(6.D0*my_A1+my_B2)
     h15_Shdv=-9.D0*my_A2*my_A2
     h15_Shdvh=-324.D0*my_A1*my_A1*my_A2*(my_A1+my_A2)
     h15_Shdvv=648.D0*my_A1*my_A1*my_A2*my_A2
     h15_Ssn0=my_A1*(1-6.D0*my_A1/my_B1)
     h15_Ssdh=-9.D0*my_A1*my_A2
     h15_Ssdv=-9.D0*my_A1*my_A1
     h15_Smn0 = my_A1*(1-6.D0*my_A1/my_B1-3.D0*my_C1)
     h15_SmnhSh = 9.D0*my_A1*(2.D0*my_A1+my_A2)
     h15_SmnsSs = 27.D0*my_A1*my_A1
     h15_Smdh = -9.D0*my_A1*my_A2
     h15_Smdv = -36.D0*my_A1*my_A1

     tmp0 = 4.D0/my_B1
     do i=1,nlev-1
! convert nondimensional forcing functions to q2-q2l formulation, 
! at leastuntil this is rederived in k-epsilon formulation
        Gh = -tmp0*an(i)
        Gv =  tmp0*av(i)
        Gs =  tmp0*aw(i)

! crude Harcourt(2015) limits. 
        Gh=max(Gh,h15_Ghmin)
        Gv=max(Gv,h15_Ghmin)
        Gh=min(Gh,h15_Ghmax)
        Gv=min(Gv,h15_Gvmax)
        Gs=min(Gs,h15_Gsmax)

        Sh=(h15_Shn0+h15_Shnh*Gh+h15_Shns*Gs+h15_Shnv*Gv)/                              &
               ((1.D0+h15_Shdah*Gh+h15_Shdav*Gv)*                                       &
                (1.D0+h15_Shdbh*Gh)+                                                    &
                (h15_Shdv+h15_Shdvh*Gh+h15_Shdvv*Gv)*Gv)
        Sh=max(0.D0,Sh)
        Ss=h15_Ssn0/(1.D0+h15_Ssdh*Gh+h15_Ssdv*Gv)
        Ss=max(0.D0,Ss)
        Sm = (h15_Smn0+h15_SmnhSh*Gh*Sh+h15_SmnsSs*Gs*Ss)/                              &
                 (1.D0+h15_Smdh*Gh+h15_Smdv*Gv)
        Sm=max(0.D0,Sm)

        Ss=Ss*SPF

        cmue1(i) =  sqrt(2.D0)*Sm
        cmue2(i) =  sqrt(2.D0)*Sh

        tmp1 = sqrt(2.D0*tke(i))*L(i)*Ss
        gamu(i) = gamu(i) - tmp1*dusdz
        gamv(i) = gamv(i) - tmp1*dvsdz

     end do

     return
   end subroutine cmue_d_h15


!-----------------------------------------------------------------------
