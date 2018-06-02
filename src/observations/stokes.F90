#include"cppdefs.h"
!/ ------------------------------------------------------------------- /
    MODULE STOKES
!/  ------------------------------------------------------------------- /
!/ Contains routines for reading in Stokes drift and calculating
!/ Langmuir numbers for use in GOTM
!/  August - 04 - 2015
!/
!/   List of subroutines:
!/    1. Init Stokes
!/    2. Calc Stokes
!/    3. Get_LaNum
!/    4. Lagrangian_Profile
!/ ------------------------------------------------------------------- /
      Private
!/ ------------------------------------------------------------------- /
        REALTYPE, ALLOCATABLE, DIMENSION(:) :: JDATESPC
! JDATESPC - Julian Date for spectra
        REALTYPE, ALLOCATABLE, DIMENSION(:,:,:) :: SPC
! SPC - Spectra (f'n of wavenum, dir, and time)
        REALTYPE, ALLOCATABLE, DIMENSION(:,:) :: STOKESx,STOKESy
! STOKESx(/y) - Stokes drift in x(/y) direction (f'n of depth, time)
        REALTYPE, ALLOCATABLE, DIMENSION(:) :: WD, WN,Z_CEN
! WD - Wave directions
! WN - Wavenumbers
! Z_CEN - Depth of center
        REALTYPE, ALLOCATABLE, DIMENSION(:) :: LVLTHCK,Z_FACE
        REALTYPE, public, allocatable, dimension(:) :: US, VS
! US - Stokes U velocity
! VS - Stokes V velocity
! LVLTHCK - Thickness of level
! Z_FACE - Depth of face
        INTEGER :: NSPC, NDIR, NK, NZ
! NSPC - number spectra (total times)
! NDIR - number of spectral directions
! NK   - numer of spectral wavenumbers
! NZ   - number of levels
        REALTYPE, PARAMETER :: TPI=2*3.141592653589793

        public init_stokes, InterpStokesProfile, Get_LaNum
!/ ------------------------------------------------------------------- /
        CONTAINS
!/ ------------------------------------------------------------------- /
        SUBROUTINE INIT_STOKES(Stokes_Switch, Stokes_WL)
!/ ------------------------------------------------------------------- /
!/ Name    : INIT_STOKES
!/ Purpose : Read inputs and calculate time series of Stokes drift
!/           Profiles for use in GOTM.
!/ Author  : Brandon Reichl (URI-GSO)
!/ Intent(in) : Stokes_Switch - Switch for calculating Stokes drift
!/ Intent(in) : Stokes_WL - Upper wavelength limit for Stokes drift
!/ ------------------------------------------------------------------- /
          USE time, ONLY: read_time_string
!/ ------------------------------------------------------------------- /
          IMPLICIT NONE
!/ ------------------------------------------------------------------- /
          LOGICAL, intent(in)  :: Stokes_Switch
          REALTYPE, intent(in) :: Stokes_WL
!/ ------------------------------------------------------------------- /
          REALTYPE    :: DATA
          INTEGER :: T, I, D
          INTEGER :: YR,MO,DA,HR,MI,SC,Day,Sec
          CHARACTER(19) :: ST
!/
!/ ------------------------------------------------------------------- /
!/
!/ - Open Depth used in GOTM and read in (should simply access via GOTM)
!/
          OPEN(82,FILE='Depth.txt',status='old')
          READ(82,*) NZ
          ALLOCATE(Z_CEN(NZ),LVLTHCK(NZ),Z_FACE(NZ))
          READ(82,*) DATA
          Z_CEN(NZ)=-DATA*.5
          Z_FACE(NZ)=0.0
          LVLTHCK(NZ) = DATA
          DO I=1, NZ-1
             READ(82,*) DATA
             LVLTHCK(NZ-I)=DATA
             Z_FACE(NZ-I)=Z_FACE(NZ-I+1)-DATA
             Z_CEN(NZ-I)=Z_CEN(NZ-I+1)-DATA
          ENDDO
!/
!/ - Open Spectrum
!/
          OPEN(81,FILE='SPC.txt',status='old')
          REWIND(81)
          READ(81,*)NK,NDIR,NSPC
!/
!/ - Allocate Spectral variables based on spectral header
!/
          IF (.NOT. ALLOCATED(WD))      ALLOCATE(WD(NDIR))
          IF (.NOT. ALLOCATED(WN))      ALLOCATE(WN(NK))
          IF (.NOT. ALLOCATED(SPC))    ALLOCATE(SPC(NSPC,NK,NDIR))
          IF (.NOT. ALLOCATED(JDATESPC)) ALLOCATE(JDATESPC(NSPC))
          IF (.NOT. ALLOCATED(STOKESx)) ALLOCATE(STOKESx(NSPC,0:NZ))
          STOKESx = _ZERO_
          IF (.NOT. ALLOCATED(STOKESy)) ALLOCATE(STOKESy(NSPC,0:NZ))
          STOKESy = _ZERO_
          IF (.NOT. ALLOCATED(US)) ALLOCATE(US(0:NZ));US = _ZERO_
          IF (.NOT. ALLOCATED(VS)) ALLOCATE(VS(0:NZ));VS = _ZERO_
!/
!/ - Read in spectra and calculate Stokes drift
!/
          READ(81,*)WN
          READ(81,*)WD
          DO D=1,NSPC
             READ(81,'(i4.4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)') &
                  YR,MO,DA,HR,MI,SC
             WRITE(ST,'(i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1, &
                  i2.2,a1,i2.2)')YR,'-',MO,'-',DA,' ',HR,':',MI,':',SC
             CALL READ_TIME_STRING(ST,Day,Sec)
             JDATESPC(D)=real(sec)/86400.
             JDATESPC(D)=real(day)+JDATESPC(D)
             if (STOKES_SWITCH) then
                DO I=1,NK
                   READ(81,*)SPC(D,I,:)
                ENDDO
             else
                SPC=0.0 !Waves set to 0.
             end if

             CALL CalcStokes(SPC(D,:,:),STOKESx(D,:),STOKESy(D,:),&
                  TPI/Stokes_WL)

          ENDDO
!/
          CLOSE(81)
!/
          RETURN
!/
        ENDSUBROUTINE init_STOKES
!/
!/ ------------------------------------------------------------------- /
        SUBROUTINE CalcStokes(SPCin,STKx,STKy,WNUL)
!/ ------------------------------------------------------------------- /
!/ Name    : CalcStokes
!/ Purpose : Receives Wave spectrum and returns Stokes drift profile.
!/           Stokes drift is integrated and averaged over each layer.
!/ Author  : Brandon Reichl (URI-GSO)
!/ Intent(in) : SPCin - Wave spectrum
!/ Intent(in) : WNUL - Upper wavelength limit for Stokes drift
!/ Intent(out): STKx - X-component Stokes drift profile
!/ Intent(out): STKy - Y-component Stokes drift profile
!/ ------------------------------------------------------------------- /
          IMPLICIT NONE
!/
          REALTYPE, INTENT(IN), DIMENSION(NK,NDIR) :: SPCin
          REALTYPE, INTENT(OUT), DIMENSION(0:NZ) :: STKx, STKy
          REALTYPE, INTENT(IN) :: WNUL
!/
          REAL, PARAMETER :: GRAV=9.806
          REALTYPE    :: DTH, DKTH,SQG,HS
          INTEGER :: I,J,Z,K,T
!/
          STKx = _ZERO_;
          STKy = _ZERO_;
          SQG = sqrt(GRAV)
          DTH=TPI/REAL(NDIR)
          HS=0.0
          STKx(1) = 0.0
          STKy(1) = 0.0
          DO Z=2,NZ
             STKx(Z)=0.0
             STKy(Z)=0.0
             DO K=1,NK
                if (K.GT.1 .AND. K.LT.NK) THEN
                   DKTH=(WN(K+1)-WN(K-1))/2.*DTH
                elseif (K.EQ.1) then
                   DKTH=(WN(K+1)-WN(K))*DTH
                else
                   DKTH=(WN(K)-WN(K-1))*DTH
                end if
                   IF (2*3.1415/WN(k).gt.1.0) then
                   DO T=1,NDIR
                      STKx(Z) = STKx(Z) + SPCin(K,T) * SQRT(GRAV*WN(K))*COS(WD(T)) &
                           * (EXP(2.*WN(K)*(Z_FACE(Z))) - EXP(2.*WN(K)*Z_FACE(Z-1))) &
                           * DKTH
                      STKy(Z) = STKy(Z) + SPCin(K,T) * SQRT(GRAV*WN(K))*SIN(WD(T)) &
                           * (EXP(2.*WN(K)*(Z_FACE(Z))) - EXP(2.*WN(K)*Z_FACE(Z-1))) &
                           * DKTH
                   ENDDO
                ENDIF
             ENDDO
             STKx(Z) = STKx(Z) / LVLTHCK(z)
             STKy(Z) = STKy(Z) / LVLTHCK(z)
          ENDDO
          HS = 4. * sqrt(HS)

          RETURN
        END SUBROUTINE CalcStokes
!/
!/ ------------------------------------------------------------------- /
!/
        SUBROUTINE InterpStokesProfile()
!/ ------------------------------------------------------------------- /
!/ Name    : LagrangianProfile
!/ Purpose : Computes Lagrangian currents by interpolating Stokes
!/           drift profile in time.
!/ Author  : Brandon Reichl (URI-GSO)
!/ ------------------------------------------------------------------- /
          USE TIME, only  : julianday,secondsofday
!/
          implicit none
          REALTYPE              :: JDAYIn
          REALTYPE              :: DT,DZ
          INTEGER               :: I,Z
!/
          JDAYIn=real(secondsofday/86400.)
          JDAYIN=JDAYIN+REAL(JULIANDAY)
          US(:) = 0
          VS(:) = 0
          IF (JDayIn.GT.JDATESPC(1).AND.JDayIn.LT.JDATESPC(Nspc)) THEN
             I=1
             DO WHILE(JDayIN.GT.JDATESPC(I+1))
                I=I+1
             ENDDO
             DT=JDATESPC(I+1)-JDATESPC(I)
             DO Z=1,NZ
                US(Z) = Stokesx(I+1,Z) * (JDAYIN-JDATESPC(I))/DT &
                     +Stokesx(I,Z)* (JDATESPC(I+1)-JDAYIN)/DT
                VS(Z) = Stokesy(I+1,Z) * (JDAYIN-JDATESPC(I))/DT &
                     +Stokesy(I,Z)* (JDATESPC(I+1)-JDAYIN)/DT
             ENDDO
          ELSE
             US(:) = 0.0
             VS(:) = 0.0
          ENDIF
!/
          RETURN
!/
        END SUBROUTINE InterpStokesProfile
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
    ENDMODULE STOKES
