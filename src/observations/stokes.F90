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
!/    3. Get_LaNum (Relocated June, 2018 to turbulence/Langmuir)
!/    4. Lagrangian_Profile
!/ ------------------------------------------------------------------- /
      Private
!/ ------------------------------------------------------------------- /
        REALTYPE, ALLOCATABLE, DIMENSION(:), public :: JDATESPC
! JDATESPC - Julian Date for spectra
        REALTYPE, ALLOCATABLE, DIMENSION(:,:,:), public :: SPC
! SPC - Spectra (f'n of wavenum, dir, and time)
        REALTYPE, ALLOCATABLE, DIMENSION(:,:), public :: STOKESx,STOKESy
! STOKESx(/y) - Stokes drift in x(/y) direction (f'n of depth, time)
        REALTYPE, ALLOCATABLE, DIMENSION(:), public :: WD, WN,Z_CEN
! WD - Wave directions
! WN - Wavenumbers
! Z_CEN - Depth of center
        REALTYPE, ALLOCATABLE, DIMENSION(:), public :: LVLTHCK,Z_FACE
        REALTYPE, allocatable, dimension(:), public :: US, VS
! US - Stokes U velocity
! VS - Stokes V velocity
! LVLTHCK - Thickness of level
! Z_FACE - Depth of face
        INTEGER, public :: NSPC, NDIR, NK, NZ
! NSPC - number spectra (total times)
! NDIR - number of spectral directions
! NK   - numer of spectral wavenumbers
! NZ   - number of levels
        REALTYPE, PARAMETER :: TPI=2*3.141592653589793

        public init_stokes, InterpStokesProfile
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
    ENDMODULE STOKES
