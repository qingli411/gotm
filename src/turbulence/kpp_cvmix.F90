#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: kpp: the KPP-turbulence model \label{sec:kpp}
!
! !INTERFACE:
    module kpp_cvmix
!
! !DESCRIPTION:
!  This is a wrapper of CVMix to do KPP using CVMix routines

    use cvmix_kpp

    IMPLICIT NONE

    private

! !PUBLIC MEMBER FUNCTIONS:
!
    public init_kpp_cvmix, do_kpp_cvmix, clean_kpp_cvmix
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
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!BOC
    continue
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
