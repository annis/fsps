SUBROUTINE GETMAGS(zred,spec,mags,mag_compute)

  !routine to calculate magnitudes in the Vega or AB systems,
  !given an input spectrum and redshift.
  !see parameter compute_vega_mags in sps_vars.f90
  !magnitudes defined in accordance with Fukugita et al. 1996, Eqn 7
  !This routine also redshifts the spectrum, if necessary.

  USE sps_vars; USE sps_utils, ONLY : linterp, tsum
  IMPLICIT NONE

  INTEGER  :: i
  REAL(SP), INTENT(in) :: zred
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(inout), DIMENSION(nbands) :: mags
  INTEGER, DIMENSION(nbands), INTENT(in), OPTIONAL  :: mag_compute
  INTEGER, DIMENSION(nbands) :: magflag
  REAL(SP), DIMENSION(nspec)  :: tspec

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  mags = 99.
  
  !set up the flags determining which mags are computed
  IF (PRESENT(mag_compute)) THEN
     magflag = mag_compute
     magflag(1) = 1  !always compute V band mag
  ELSE
     magflag = 1
  ENDIF

  ! JTA changing alpha/Fe via spline fit to MILES spectra ratio
  ! spec = spec * ( &
  !   + (spec_lambda**0)*(9.95400175e-01)  &
  !   + (spec_lambda**1)*(5.64599641e-05) &
  !   + (spec_lambda**2)*(-2.68732708e-08) &
  !   + (spec_lambda**3)*(4.23254854e-12) & 
  !   + (spec_lambda**4)*( -2.16644332e-16) &
  !   ) 
  ! spec = spec * ( &
  !   + (spec_lambda**0)*(-6.60525884e-01)  &
  !   + (spec_lambda**1)*(1.68594678e-03) &
  !   + (spec_lambda**2)*(-6.62213048e-07) &
  !   + (spec_lambda**3)*(1.26620043e-10) & 
  !   + (spec_lambda**4)*(-1.18475548e-14) &
  !   + (spec_lambda**5)*(4.35957547e-19) &
  !   ) 
  ! the following is the one used for the BMA paper
  ! spec = spec * ( &
  !   + (spec_lambda**2)*(1.66454871e-08 ) &
  !   + (spec_lambda**1)*(-2.01350117e-04) & 
  !   + (spec_lambda**0)*(1.63908679e+00) &
  !   ) 

  !redshift the spectrum
  IF (ABS(zred).GT.tiny_number) THEN
     DO i=1,nspec
        tspec(i) = MAX(linterp(spec_lambda*(1+zred),spec,&
        spec_lambda(i)),0.0)
     ENDDO
     spec = tspec
  ELSE
     tspec = spec  
  ENDIF

  !the units of the spectra are Lsun/Hz; convert to
  !erg/s/cm^2/Hz, at 10pc for absolute mags
  tspec = tspec*lsun/4.0/mypi/(pc2cm*pc2cm)/100.0

  !integrate over each filter
  DO i=1,nbands
     IF (magflag(i).EQ.0) CYCLE
     mags(i) = TSUM(spec_lambda,tspec*bands(:,i)/spec_lambda)
     mags(i) = MAX(mags(i),tiny_number)
  ENDDO

  !convert to magnitudes in AB system
  DO i=1,nbands
     IF (magflag(i).EQ.0) CYCLE
     IF (mags(i).LE.tiny_number) THEN
        mags(i) = 99.0 
     ELSE
        mags(i) = -2.5*LOG10(mags(i))-48.60
     ENDIF
  ENDDO

  !put magnitudes in the Vega system if keyword is set
  !(V-band is the first element in the array)
  IF (compute_vega_mags.EQ.1) &
       mags(2:nbands) = (mags(2:nbands)-mags(1)) - &
       (magvega(2:nbands)-magvega(1)) + mags(1)


END SUBROUTINE GETMAGS
