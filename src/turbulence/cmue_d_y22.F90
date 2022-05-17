#include"cppdefs.h"
!-----------------------------------------------------------------------
!
! !ROUTINE: The Langmuir turbulence quasi-equilibrium stability functions  after Harcourt (2015)\label{sec:cmueDH15}
!
! !INTERFACE:
   subroutine cmue_d_y22(nlev)
!
! !DESCRIPTION:
!
!  Old Description from GTOM: This subroutine updates the explicit solution of
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
   use turbulence, only: cmue3
   use turbulence, only: length_lim
! YCC vvv
   use turbulence, only: ct1, cc1, ctt, at3, a1, a2, a3, a5, at1, at2, at4, at5
   use turbulence, only: cm0
! YCC ^^^
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer, intent(in)       :: nlev

! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: small       = 1.0D-8
! cmue_c
   REALTYPE, parameter       :: asLimitFact=1.0d0
   REALTYPE, parameter       :: anLimitFact=0.5d0
! cmue_c

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!  Converted by Ramsey Harcourt, last updated 31 July 2018.  This version uses the Harcourt(2015)
!  stability functions from the quasi-equilibrium Second Moment closure (SMC) with Craik-Leibovich
!  terms, but it has been further modified by replacing the crude limiters applied
!  individually to Gh, Gv and Gs in Harcourt(2015) under unstable/positive production conditions
!  with a combinations of limitations on (L/q)**2 applied consistently across Gm, Gs, Gh, Gv,
!  as a function of Gh and Gv input to the ARSM. This ARSM also applies the Galperin limit to
!  L going into the ARSM (algebraic) diagnosis of Sm, Sh, Ss, regardless of whether it/s
!  being enforced within the dyanamic model.
!
!  Recomend running with e3=5 & length.lim=false, but e3=1.2 & length.lim=true is similar
!  When length.lim=false, length scale or at least L/q is still limited within the ARSM for
!  Stability fcns cmue1,cmue2,cmue3. length.lim=false allows the elevated length scale within
!  the mixed layer to impact the transition zone, while restraining the Stability functions
!  to the stability-limited length scale regime.
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
!
     integer                 ::   i
     REALTYPE            ::   Gv, Gs, Gh, Gm
     REALTYPE            ::   Sm, Ss, Sh

     REALTYPE, parameter :: my_A1 = 0.92D0
     REALTYPE, parameter :: my_A2 = 0.74D0
     REALTYPE, parameter :: my_B1 = 16.6D0
     REALTYPE, parameter :: my_B2 = 10.1D0
     REALTYPE, parameter :: my_C1 = 0.08D0
     REALTYPE, parameter :: my_C2 = 0.7D0
     REALTYPE, parameter :: my_C3 = 0.2D0
     REALTYPE, parameter :: h15_Ghmin = -0.28D0
     REALTYPE, parameter :: h15_Ghoff = 0.003D0
     REALTYPE, parameter :: h15_Gvoff = 0.006D0
     REALTYPE, parameter :: h15_Sxmax = 2.12D0

     REALTYPE :: h15_Shn0, h15_Shnh, h15_Shns, h15_Shnv
     REALTYPE :: h15_Shdah, h15_Shdav, h15_Shdbh
     REALTYPE :: h15_Shdv, h15_Shdvh, h15_Shdvv
     REALTYPE :: h15_Ssn0, h15_Ssdh, h15_Ssdv
     REALTYPE :: h15_Smn0, h15_SmnhSh, h15_SmnsSs
     REALTYPE :: h15_Smdh, h15_Smdv

     REALTYPE :: tmp0,tmp1,tmp2,tmp3,tmp4
     REALTYPE :: Ghcrit, Gvcrit
! yucc VVV
		 REALTYPE :: l1
		 REALTYPE :: l2
		 REALTYPE :: l3
		 REALTYPE :: l4
		 REALTYPE :: l5
		 REALTYPE :: l6
		 REALTYPE :: l7
		 REALTYPE :: l8
     REALTYPE :: T1, T2, T3, T4, T5
     REALTYPE :: H1, H2, H3, H4, H5, H6
     REALTYPE :: F, F1, F2
     REALTYPE :: ab3, N, Nb
     REALTYPE, parameter :: Gh_off = 0.06D0
     REALTYPE, parameter :: Gv_off = 0.12D0
     REALTYPE, parameter :: Gh_min = -19.31D0
! yucc ^^^
!-----------------------------------------------------------------------
! yucc VVV
     N    =   0.5*cc1
     Nb   =   ct1
     ab3   =   at3
     l1=2*a1/N
     l2=a2/N
     l3=a3/N
     l4=2*a5/N
     l6=at1/Nb
     l7=at2/Nb
     l8=at4*ctt/Nb
		 T1=2*l2/3+2*l3+l6+l7
		 T2=2*l2/3-2*l3+l6-l7
		 T3=l2**2/3-3*l3**2
		 T4=l2**2/3+6*l2*l3+3*l3**2
		 T5=l2**2/3-6*l2*l3+3*l3**2
! yucc new


		do i=1,nlev-1
      	
				Gh = -an(i)
				Gh = max(Gh,-1.4045) ! -0.28
				Gh = min(Gh,0.145) ! 0.029
				Gm = as(i)
				Gm = min(Gm,5+1.875*Gh)
				Gv = av(i)
				Gs = aw(i)

				F=1.D0-1.D0/3.D0*(l2**2-3*l3**2)*(Gm+Gs)-2.D0/3.D0*(l2**2+3*l3**2)*Gv-1.D0/2.D0*l4*ab3/Nb*Gh
				F1=F**2+2*l2**2*l3**2*(Gm*Gs-Gv**2)+1.D0/4.D0*F*((1.D0/3.D0*l2**2-3*l3**2)*(Gm+Gs)+2*(1.D0/3.D0*l2**2+3*l3**2)*Gv)
				F2=1.D0-1.D0/4.D0*(l6**2-l7**2)*(Gm+Gs)-1.D0/2.D0*(l6**2+l7**2)*Gv-(2*ab3/3.D0/Nb*l4+ctt*l8)*Gh
								
				H1=1.D0/8.D0*l1*(4*F+T3*Gm+T4*Gs+(T3+T4)*Gv)
				H2=1.D0/8.D0*l4*(4*F*T1+T3*T1*Gm+T4*T2*Gs+(T3*T2+T4*T1)*Gv)*Gh
				H3=1.D0/8.D0*l1*(4*F+T5*Gm+T3*Gs+(T3+T5)*Gv)
				H4=1.D0/8.D0*l4*(4*F*T2+T5*T1*Gm+T3*T2*Gs+(T3*T1+T5*T2)*Gv)*Gh
				H5=ab3/4.D0/Nb*(T2*Gm+T1*Gv)
				H6=ab3/4.D0/Nb*(T2*Gv+T1*Gs)
				
				tmp1=1.D0/4.D0*l1*l4*(T2-T1)*H6*Gh-ab3/3.D0/Nb*H2-H1*F2
				tmp2=-1.D0/4.D0*l1*l4*(T2-T1)*H5*Gh-ab3/3.D0/Nb*H4-H3*F2
				tmp3=-ab3/3.D0/Nb*F1-H1*H5-H3*H6
				tmp0=H2*H5+H4*H6-F1*F2

				Sm=tmp1/tmp0
				Ss=tmp2/tmp0
				Sh=tmp3/tmp0
				Sm=min(max(small,Sm),h15_Sxmax)
				Ss=min(max(small,Ss),h15_Sxmax)
				Sh=min(max(small,Sh),h15_Sxmax)
				cmue1(i) =  1.D0/((0.3)**(3.D0/2.D0))*Sm
				cmue2(i) =  1.D0/((0.3)**(3.D0/2.D0))*Sh
				cmue3(i) =  1.D0/((0.3)**(3.D0/2.D0))*Ss

		enddo
		
     return
! yucc ******

     return
   end subroutine cmue_d_y22


!-----------------------------------------------------------------------
