#include"cppdefs.h"
!/ ------------------------------------------------------------------- /
    MODULE Langmuir
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
      Use Stokes, only :&!
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
      private

      ! Right now these parameters are only used in KPP to set for output
      ! they could be set here.
      REALTYPE, public :: Mixing_Efactor !< Output enhancement to mixing
                                         !! coefficient
      REALTYPE, public :: Entrainment_Efactor !< Output enhancement to Vt2
                                              !! if LF16
      REALTYPE, public :: LA_to_Vt2 !< Output Langmuir number used to 
                                    !! get Vt2 if RWHGK16 or LF17

      public Get_LaNum
!/ ------------------------------------------------------------------- /
    CONTAINS
!/
      SUBROUTINE Get_LaNum(AVGDPTH,USTin,LAOUT)
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
        USE TIME, only  : julianday,secondsofday
        use meanflow, only: u, v
!/
        implicit none
        REALTYPE, INTENT(IN)  :: AVGDPTH,USTin
        REALTYPE, INTENT(OUT) :: LAOUT
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
           LAdir_L=max(LAdir_L,1.E-8)
           LAOUT = sqrt (USTin / LA_LOWER) * min(10.0,(sqrt(_ONE_/LAdir_L)))
        ELSE
           LAOUT=100.0!Arbitrary large number for 'no-wave' condition.
        ENDIF
      END SUBROUTINE Get_LaNum
    ENDMODULE LANGMUIR
