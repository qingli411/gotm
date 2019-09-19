#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stokes_drift --- Stokes drift \label{sec:stokes_drfit}
!
! !INTERFACE:
   module stokes_drift
!
! !DESCRIPTION:
!
! !USES:
   use observations, only: us_prof_method, nfreq, wave_age, pi
   use airsea, only: u10, v10

   use Stokes, only: init_stokes, InterpStokesProfile
   use Stokes, only: US,VS

   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public do_stokes_drift

! !DEFINED PARAMETERS:

!  pre-defined parameters
   integer, parameter        :: NOTHING=0
   integer, parameter        :: CONSTANT=1
   integer, parameter        :: FROMFILE=2
   integer, parameter        :: FROMSPEC=3
   integer, parameter        :: FROMUSP=4
   integer, parameter        :: EXPONENTIAL=5
   integer, parameter        :: THEORYWAVE=6
   integer, parameter        :: HURRSPEC=7
   integer, parameter        :: DHH85SPEC=8
   REALTYPE, parameter       :: gravity = 9.81

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
! !IROUTINE: do_stokes_drift
!
! !INTERFACE:
   subroutine do_stokes_drift(freq,spec,xcmp,ycmp,nlev,z,zi, &
                           us_x,us_y,delta,ustokes,vstokes, &
                           dusdz,dvsdz)
!
! !DESCRIPTION:
!  A wrapper for all the subroutines to calculate the Stokes drift profile.
! TODO: Documentation  <11-03-18, Qing Li> !
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: freq(:), spec(:), xcmp(:), ycmp(:)
   REALTYPE, intent(in)                :: z(0:nlev), zi(0:nlev)
   REALTYPE, intent(inout)             :: us_x, us_y, delta
!
! !OUTPUT PARAMETERS:
   ! grid cell averaged Stokes drift
   REALTYPE, intent(out)               :: ustokes(0:nlev), vstokes(0:nlev)
   ! vertical shear of Stokes drift at cell interface
   REALTYPE, intent(out)               :: dusdz(0:nlev), dvsdz(0:nlev)

! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: ustran
!-----------------------------------------------------------------------
!BOC

   select case (us_prof_method)
      case (NOTHING)
         ustokes = _ZERO_
         vstokes = _ZERO_
         us_x = _ZERO_
         us_y = _ZERO_
         delta = _ZERO_
         dusdz = _ZERO_
         dvsdz = _ZERO_
      case (CONSTANT)
!        Empirical spectrum
         LEVEL1 'The following us_prof_method has yet been supported: ', us_prof_method
         stop 'stokes_drift()'
      case (FROMFILE)
         ! ustokes and vstokes already read from file
         call surface_stokes_drift(ustokes,vstokes,nlev,zi,us_x,us_y,delta)
      case (FROMSPEC)
         call stokes_drift_spec(freq,spec,xcmp,ycmp,nlev,z,zi,us_x,us_y,delta,ustokes,vstokes)
      case (FROMUSP)
         call stokes_drift_usp(freq,xcmp,ycmp,nlev,z,zi,us_x,us_y,delta,ustokes,vstokes)
      case (EXPONENTIAL)
         call stokes_drift_exp(nlev,z,zi,us_x,us_y,delta,ustokes,vstokes)
      case (THEORYWAVE)
         ! TODO: wrap in a subroutine, compute surface Stokes drift and decay depth <20190116, Qing Li> !
         ustokes(:) = 0.0
         vstokes(:) = 0.0
         do i=1,nlev
            call stokes_drift_theory(u10,zi(i),zi(i-1),ustokes(i))
            call stokes_drift_theory(v10,zi(i),zi(i-1),vstokes(i))
         enddo
         call surface_stokes_drift(ustokes,vstokes,nlev,zi,us_x,us_y,delta)
      case (HURRSPEC)
         ! TODO: wrap in a subroutine <20190116, Qing Li> !
         ustokes(:) = 0.0
         vstokes(:) = 0.0
         call interpstokesprofile()
         ustokes(:) = US(:)
         vstokes(:) = VS(:)
         call surface_stokes_drift(ustokes,vstokes,nlev,zi,us_x,us_y,delta)
      case (DHH85SPEC)
         call stokes_drift_dhh85(nlev,z,zi,u10,v10,wave_age,us_x,us_y,delta,ustokes,vstokes)
      case default
         LEVEL1 'A non-valid us_prof_method has been given ', us_prof_method
         stop 'stokes_drift()'
   end select

   if (us_prof_method .ne. NOTHING) then
      do i=1,nlev-1
         dusdz(i) = (ustokes(i+1)-ustokes(i))/(z(i+1)-z(i))
         dvsdz(i) = (vstokes(i+1)-vstokes(i))/(z(i+1)-z(i))
      end do
      dusdz(0   ) = dusdz(1     )
      dusdz(nlev) = dusdz(nlev-1)
      dvsdz(0   ) = dvsdz(1     )
      dvsdz(nlev) = dvsdz(nlev-1)
   end if

   end subroutine do_stokes_drift
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surface_stokes_drift
!
! !INTERFACE:
   subroutine surface_stokes_drift(ustokes,vstokes,nlev,zi,us_x,us_y,delta)
!
! !DESCRIPTION:
!  Get surface Stokes drift and Stokes drift decay depth from Stokes drift profiles
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)      :: nlev
   REALTYPE, intent(in)     :: ustokes(nlev), vstokes(nlev)
   REALTYPE, intent(in)     :: zi(0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)    :: us_x, us_y, delta

! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: ustran
!-----------------------------------------------------------------------
!BOC

   ! set the surface Stokes drift to be the value at the upper most level
   us_x = ustokes(nlev)
   us_y = vstokes(nlev)
   ! Stokes transport
   ustran = _ZERO_
   do i=1,nlev
      ustran = ustran + sqrt(ustokes(i)**2+vstokes(i)**2)*(zi(i)-zi(i-1))
   end do
   delta = ustran/max(SMALL, sqrt(us_x**2.+us_y**2.))

   end subroutine surface_stokes_drift
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_drift_spec
!
! !INTERFACE:
   subroutine stokes_drift_spec(freq,spec,xcmp,ycmp,nlev,z,zi,&
                                us_x,us_y,delta,ustokes,vstokes)
!
! !DESCRIPTION:
!  Calculate the Stokes drift profile from wave spectrum. The wave spectrum
!   here is the spectrum weighted by the frequency bin width,
!   i.e., spec = S(f) df
! TODO: Documentation  <20-12-17, Qing Li> !
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: freq(:), spec(:), xcmp(:), ycmp(:)
   REALTYPE, intent(in)                :: z(0:nlev), zi(0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: us_x, us_y, delta
   REALTYPE, intent(out)               :: ustokes(0:nlev), vstokes(0:nlev)

! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   integer                             :: i, k
   REALTYPE                            :: const, tmp, ustran
   REALTYPE                            :: dz, kdz, freqc, dfreqc
   REALTYPE                            :: aplus, aminus, iplus, iminus
   REALTYPE                            :: factor(nfreq), factor2(nfreq)
!-----------------------------------------------------------------------
!BOC
!  initialization
   ustokes = _ZERO_
   vstokes = _ZERO_
   us_x = _ZERO_
   us_y = _ZERO_
   delta = _ZERO_
!  some factors
   const = 8.*pi**2./gravity
   do i=1,nfreq
      factor2(i) = const*freq(i)**2.
      factor(i) = 2.*pi*freq(i)*factor2(i)
   end do
! cutoff frequency
   freqc = 1.5*freq(nfreq)-0.5*freq(nfreq-1)
   dfreqc = freq(nfreq)-freq(nfreq-1)
!  Stokes drift averaged over the grid cell, z(0) is not used
   do k=1,nlev
      dz = zi(k)-zi(k-1)
      do i=1,nfreq
         kdz = factor2(i)*dz/2.
         if (kdz .lt. 100.) then
             tmp = sinh(kdz)/kdz*factor(i)*spec(i)*exp(factor2(i)*z(k))
         else
             tmp = factor(i)*spec(i)*exp(factor2(i)*z(k))
         end if
         ustokes(k) = ustokes(k)+tmp*xcmp(i)
         vstokes(k) = vstokes(k)+tmp*ycmp(i)
      end do
!     add contribution from a f^-5 tail (aplus > 0)
      aplus  = max(SMALL, -const*freqc**2.*zi(k))
      aminus = -const*freqc**2.*zi(k-1)
      iplus  = 2.*aplus/3.*(sqrt(pi*aplus)*erfc(sqrt(aplus)) \
              -(_ONE_-0.5/aplus)*exp(-aplus))
      iminus = 2.*aminus/3.*(sqrt(pi*aminus)*erfc(sqrt(aminus)) \
              -(_ONE_-0.5/aminus)*exp(-aminus))
      tmp = 2.*pi*freqc**2./dz*spec(nfreq)/dfreqc*(iplus-iminus)
      ustokes(k) = ustokes(k)+tmp*xcmp(nfreq)
      vstokes(k) = vstokes(k)+tmp*ycmp(nfreq)
   end do
!  Surface Stokes drift and penetration depth
   ustran = _ZERO_
   do i=1,nfreq
      us_x = us_x+factor(i)*spec(i)*xcmp(i)
      us_y = us_y+factor(i)*spec(i)*ycmp(i)
      ! Stokes transport
      ustran = ustran+2.*pi*freq(i)*spec(i)
   end do
!  add contribution from a f^-5 tail
   tmp = 2.*pi*const*freqc**4*spec(nfreq)/dfreqc
   us_x = us_x+tmp*xcmp(nfreq)
   us_y = us_y+tmp*ycmp(nfreq)
   ustran = ustran+4./3.*pi*freqc**2.*spec(nfreq)/dfreqc
   delta = ustran/max(SMALL, sqrt(us_x**2.+us_y**2.))

   end subroutine stokes_drift_spec
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_drift_usp
!
! !INTERFACE:
   subroutine stokes_drift_usp(freq,ussp,vssp,nlev,z,zi,&
                               us_x,us_y,delta,ustokes,vstokes)
!
! !DESCRIPTION:
!  Calculate the Stokes drift profile from partitioned Stokes drift.
! TODO: Documentation  <11-03-18, Qing Li> !
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: freq(:), ussp(:), vssp(:)
   REALTYPE, intent(in)                :: z(0:nlev), zi(0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: us_x, us_y, delta
   REALTYPE, intent(out)               :: ustokes(0:nlev), vstokes(0:nlev)

! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   integer                             :: i, k
   REALTYPE                            :: const, tmp, ustran, vstran
   REALTYPE                            :: kdz, dz
   REALTYPE                            :: factor(nfreq)
!-----------------------------------------------------------------------
!BOC
!  initialization
   ustokes = _ZERO_
   vstokes = _ZERO_
   us_x = _ZERO_
   us_y = _ZERO_
   delta = _ZERO_
   ustran = _ZERO_
   vstran = _ZERO_
!  some factors
   const = 8.*pi**2./gravity
   do i=1,nfreq
      ! 2k
      factor(i) = const*freq(i)**2.
   end do
!  Stokes drift averaged over the grid cell, z(0) is not used
   do k=1,nlev
      dz = zi(k)-zi(k-1)
      do i=1,nfreq
         kdz = factor(i)*dz/2.
         if (kdz .lt. 100.) then
             tmp = sinh(kdz)/kdz*exp(factor(i)*z(k))
         else
             tmp = exp(factor(i)*z(k))
         end if
         ustokes(k) = ustokes(k)+tmp*ussp(i)
         vstokes(k) = vstokes(k)+tmp*vssp(i)
      end do
   end do
!  Surface Stokes drift and penetration depth
   do i=1,nfreq
      us_x = us_x+ussp(i)
      us_y = us_y+vssp(i)
      ustran = ustran+ussp(i)/factor(i)
      vstran = vstran+vssp(i)/factor(i)
   end do
   delta = sqrt(ustran**2.+vstran**2.)/max(SMALL, sqrt(us_x**2.+us_y**2.))

   end subroutine stokes_drift_usp
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_drift_exp
!
! !INTERFACE:
   subroutine stokes_drift_exp(nlev,z,zi,us_x,us_y,delta,ustokes,vstokes)
!
! !DESCRIPTION:
!  Calculate the Stokes drift profile from surface Stokes drift and the
!  Stokes penetration depth, assuming exponential profile.
! TODO: Documentation  <05-04-18, Qing Li> !
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: z(0:nlev), zi(0:nlev)
   REALTYPE, intent(in)                :: us_x, us_y, delta
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: ustokes(0:nlev), vstokes(0:nlev)

! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   integer                             :: i, k
   REALTYPE                            :: tmp
   REALTYPE                            :: kdz, dz
!-----------------------------------------------------------------------
!BOC
!  initialization
   ustokes = _ZERO_
   vstokes = _ZERO_
!  Stokes drift averaged over the grid cell, z(0) is not used
   do k=1,nlev
      dz = zi(k)-zi(k-1)
      kdz = 0.5*dz/delta
      if (kdz .lt. 100.) then
          tmp = sinh(kdz)/kdz*exp(z(k)/delta)
      else
          tmp = exp(z(k)/delta)
      end if
      ustokes(k) = tmp*us_x
      vstokes(k) = tmp*us_y
   end do

   end subroutine stokes_drift_exp
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_drift_theory
!
! !INTERFACE:
   subroutine stokes_drift_theory(u10, zup, zdn, Stokes_Avg)
! !DESCRIPTION:
!   This function returns the surface layer averaged Stokes drift, given
!    the 10-meter wind (m/s) and the boundary layer depth (m).
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
! 10 meter wind (m/s)
   REALTYPE, intent(in)                 :: u10
! Top of Layer for averaging (m)
   REALTYPE, intent(in)                 :: zup
! Bottom of Layer for averaging (m)
   REALTYPE, intent(in)                 :: zdn
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: Stokes_Avg

! !REVISION HISTORY:
!  Original author(s): Qing Li
!  (Re)added to GOTM: Brandon Reichl
!
!EOP
! !LOCAL VARIABLES:
   REALTYPE, parameter :: &
   ! ratio of U19.5 to U10 (Holthuijsen, 2007)
   u19p5_to_u10 = 1.075, &
   ! ratio of mean frequency to peak frequency for
   ! Pierson-Moskowitz spectrum (Webb, 2011)
   fm_to_fp = 1.296, &
   ! ratio of surface Stokes drift to U10
   us0_to_u10 = 0.0162, &
   ! loss ratio of Stokes transport
   r_loss = 0.667
   REALTYPE :: z_store(2), stokes_store(2) !Store the depth/integrals

   REALTYPE :: us0, hm0, fm, fp, stokes_trans, kphil, kstar
   REALTYPE :: z0, z0i, r1, r2, r3, r4, tmp
   integer  :: i

   z_store(1) = zup
   z_store(2) = zdn
   stokes_store(:) = _ZERO_
   stokes_avg = _ZERO_

   if (u10 .gt. _ZERO_) then

      ! surface Stokes drift
      us0 = us0_to_u10 * u10
      !
      ! significant wave height from Pierson-Moskowitz
      ! spectrum (Bouws, 1998)
      hm0 = 0.0246 * u10**2
      !
      ! peak frequency (PM, Bouws, 1998)
      tmp = 2.0 *  pi * u19p5_to_u10 * u10
      fp = 0.877 * gravity/tmp
      !
      ! mean frequency
      fm = fm_to_fp * fp
      !
      ! total Stokes transport (a factor r_loss is applied to account
      !  for the effect of directional spreading, multidirectional waves
      !  and the use of PM peak frequency and PM significant wave height
      !  on estimating the Stokes transport)
      stokes_trans = 0.125 * pi * r_loss * fm * hm0**2
      !
      ! the general peak wavenumber for Phillips' spectrum
      ! (Breivik et al., 2016) with correction of directional spreading
      kphil = 0.176 * us0 / stokes_trans
      !
      ! surface layer averaged Stokes dirft with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! the directional spreading effect from Webb and Fox-Kemper, 2015
      ! is also included
      kstar = kphil * 2.56
      ! Stokes integral to ZUP
      do i=1,2
         z0 = abs(z_store(i))
         if (z0>1.e-4) then
            z0i = abs(_ONE_/z0)
            ! term 1 to 4
            r1 = (0.151 / kphil * z0i - 0.84) &
                 *(_ONE_ - exp(-2.0 * kphil * z0))
            r2 = -(0.84 + 0.0591 / kphil * z0i) &
                 *sqrt(2.0 * pi * kphil * z0) &
                 *erfc(sqrt(2.0 * kphil * z0))
            r3 = (0.0632 / kstar * z0i + 0.125) &
                 *(_ONE_ - exp(-2.0 * kstar * z0))
            r4 = (0.125 + 0.0946 / kstar * z0i) &
                 *sqrt(2.0 * pi * kstar * z0) &
                 *erfc(sqrt(2.0 * kstar * z0))
            stokes_store(i) = us0*(0.715 + r1 + r2 + r3 + r4)
         else
            stokes_store(i) = _ZERO_
         endif
      enddo
      stokes_avg = (stokes_store(2)*z_store(2) - stokes_store(1)*z_store(1) ) &
           /(z_store(2)-z_store(1))
   endif

   return

   end subroutine stokes_drift_theory
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_drift_dhh85
!
! !INTERFACE:
   subroutine stokes_drift_dhh85(nlev,z,zi,u10,v10,wave_age,           &
                                 us_x,us_y,delta,ustokes,vstokes)
!
! !DESCRIPTION:
!  Compute grid cell-averaged Stokes drift profile from the empirical wave
!  spectrum of Donelan et al., 1985
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: z(0:nlev), zi(0:nlev)
   REALTYPE, intent(in)                :: u10, v10, wave_age
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: us_x, us_y, delta
   REALTYPE, intent(out)               :: ustokes(0:nlev), vstokes(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   REALTYPE, parameter     :: min_omega = 0.1, max_omega = 10.0
   integer, parameter      :: nomega = 1000

   integer                 :: i, k
   REALTYPE                :: xcomp, ycomp, wind_speed
   REALTYPE                :: domega, sd_omega, tmp, dz, ustran
!-----------------------------------------------------------------------
!BOC

!  initialization
   ustokes = _ZERO_
   vstokes = _ZERO_
   ustran = _ZERO_

!  wind direction
   wind_speed = sqrt(u10**2+v10**2)
   xcomp = u10 / wind_speed
   ycomp = v10 / wind_speed
!
!  stokes drift at u-levels
   domega = ( max_omega - min_omega ) / real(nomega)
   do  k = 1, nlev
      tmp = _ZERO_
      sd_omega = min_omega + 0.5 * domega
      dz = zi(k)-zi(k-1)
      do  i = 1, nomega
         tmp = tmp +                                                       &
            domega * stokes_drift_kernel_dhh85(sd_omega, z(k), dz,         &
                                               wind_speed, wave_age)
         sd_omega = sd_omega + domega
      enddo
      ustokes(k) = xcomp * tmp
      vstokes(k) = ycomp * tmp
      ustran = ustran + dz * tmp
   enddo
   us_x = ustokes(nlev)
   us_y = vstokes(nlev)
   delta = ustran/max(SMALL, sqrt(us_x**2.+us_y**2.))

   end subroutine stokes_drift_dhh85
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_drift_kernel_dhh85
!
! !INTERFACE:
   function stokes_drift_kernel_dhh85(sd_omega, sd_z, sd_dz, wind_speed, wave_age)
!
! !DESCRIPTION:
!  Evaluate kernel of the Stokes integral for the Donelan et al., 1985 spectrum
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)               :: sd_omega, sd_z, sd_dz
   REALTYPE, intent(in)               :: wind_speed, wave_age
!
! !REVISION HISTORY:
!  Original author(s): Qing Li
!
!EOP
! !LOCAL VARIABLES:
   REALTYPE  :: dhh_omega_p, dhh_alpha, dhh_sigma, dhh_gamma1, dhh_gamma2
   REALTYPE  :: wave_spec, sd_filter, kdz, iwa
   REALTYPE  :: stokes_drift_kernel_dhh85
!-----------------------------------------------------------------------
!BOC

!  DHH 85 spectrum
   iwa = _ONE_ / wave_age !< inverse wave age
   dhh_omega_p = gravity * iwa / wind_speed !< peak frequency
   dhh_alpha  = 0.006 * iwa**(0.55)
   dhh_sigma  = 0.08 * ( _ONE_ + 4.0 * wave_age**3 )
   if ( iwa .le. _ONE_) then
      dhh_gamma1 = 1.7
   else
      dhh_gamma1 = 1.7 + 6.0 * log10( iwa )
   endif
   dhh_gamma2 = exp( -0.5 * ( sd_omega - dhh_omega_p )**2 /              &
                    dhh_sigma**2 / dhh_omega_p**2 )
   wave_spec  = dhh_alpha * gravity**2 / (dhh_omega_p * sd_omega**4 ) *  &
                exp( -( dhh_omega_p / sd_omega )**4 ) *                  &
                dhh_gamma1**dhh_gamma2
!
!--Stokes drift integral kernel
   kdz = sd_omega**2 * sd_dz / gravity
   if ( kdz .lt. 10.0 ) then
      sd_filter = sinh(kdz) / kdz
   else
      sd_filter = _ONE_
   endif
   stokes_drift_kernel_dhh85 = 2.0 * ( wave_spec * sd_omega**3 ) *       &
          sd_filter * exp( 2.0 * sd_omega**2 * sd_z / gravity ) / gravity
   return

   end function stokes_drift_kernel_dhh85
!EOC

!-----------------------------------------------------------------------

   end module stokes_drift

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
