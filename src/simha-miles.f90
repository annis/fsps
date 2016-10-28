PROGRAM SIMHA

  !set up modules
  USE sps_vars; USE sps_utils
  
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90
  INTEGER :: i,j,k,l,m
  !define variable for SSP spectrum
  REAL(SP), DIMENSION(ntfull,nspec)  :: spec_pz
  !define variables for Mass and Lbol info
  REAL(SP), DIMENSION(ntfull)    :: mass_pz,lbol_pz
  CHARACTER(100) :: filename=''
  CHARACTER(100) :: format=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp
  REAL(SP) :: zave

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Now we're going to show you how to use full  
  ! metallicity-dependent info                   
  !
  ! modified from lesssimple.f90 to implement the simha-2014 models
  !
  ! This version is to be used with the Miles+Padova spectra+isochrones
  ! These have 5 metallicities, not 20.
  
  pset%sfh = 5    ! compute simha model
  pset%sf_theta = -0.52  ! assuming radians, this is 30 degrees
  pset%sf_trunc = 11.0
  pset%sf_start = 1.0
  pset%tau = 1.0


  ! part of the get LRG colors right
  pset%fbhb = 0.2

  !here we have to read in all the libraries
  CALL SPS_SETUP(-1)

  !compute all SSPs (i.e. at all Zs)
  !nz and the various *ssp_zz arrays are stored 
  !in the common block set up in sps_vars.f90
  DO i=1,nz
     pset%zmet = i
     CALL SSP_GEN(pset,mass_ssp_zz(:,i),&
          lbol_ssp_zz(:,i),spec_ssp_zz(:,:,i))
  ENDDO

  !run compsp 


  do i=1,5
     if (i == 1 ) pset%zmet = 5  ! Z = 0.030
     if (i == 2 ) pset%zmet = 4  ! solar metalicity, Z  = 0.0190
     if (i == 3 ) pset%zmet = 2  ! Z = 0.0031  ! log Z = -0.78
     if (i == 4 ) pset%zmet = 3  ! Z = 0.0096
     if (i == 5 ) pset%zmet = 1  ! Z = 0.0008  low metallicty pop
     do j=1,4
     	if (j == 1 ) pset%sf_start = 2.0
     	if (j == 2 ) pset%sf_start = 1.5
     	if (j == 3 ) pset%sf_start = 1.0
     	if (j == 4 ) pset%sf_start = 0.7
     	do k=1,4
     	   if (k == 1 ) pset%sf_trunc = 7
     	   if (k == 2 ) pset%sf_trunc = 9
     	   if (k == 3 ) pset%sf_trunc = 11
     	   if (k == 4 ) pset%sf_trunc = 13
     	   do l=1,7
     	      if (l == 1 ) pset%tau = 0.3
     	      if (l == 2 ) pset%tau = 0.7
     	      if (l == 3 ) pset%tau = 1.0
     	      if (l == 4 ) pset%tau = 1.3
     	      if (l == 5 ) pset%tau = 2.0
     	      if (l == 6 ) pset%tau = 9.0
     	      if (l == 7 ) pset%tau = 13.0
     	      do m=1,5
     	         if (m == 1 ) pset%sf_theta = -1.047  ! 60 degrees
     	         if (m == 2 ) pset%sf_theta = -0.175  ! 10 degrees
     	         if (m == 3 ) pset%sf_theta = -0.524  ! 30 degrees
     	         if (m == 4 ) pset%sf_theta = -0.785  ! 45 degrees
     	         if (m == 5 ) pset%sf_theta = -1.396  ! 80 degrees

                 ! define redshift
                 ! redshift_colors
                 ! 0 = colors redshifted to a fixed redshift, specified in parameter set
                 ! 1 = colors redshifted according to the age of the SSP or CSP
                 pset%redshift_colors=1
                 pset%zred=0.0

         ! if (pset%zmet == 10 ) continue

		 ! '(A2,    I2,    A1,  F3.1,         A1,    I1,             A1, F3.1,         F6.3,    A2)' 
	         ! 's-',pset%zmet,'-', pset%sf_start,"-",int(pset%sf_trunc),"-",pset%tau,pset%sf_theta,'-z'
                 !         123456789112345678921234567893123456
		 format = '(A2,I2,A1,F3.1,A1,I1,A1,F3.1,F6.3,A2)'
		 IF (pset%zmet >= 10)     WRITE(format(5:6),'(A2)') 'I2'
		 IF (pset%zmet < 10)      WRITE(format(5:6),'(A2)') 'I1'
		 IF (pset%sf_trunc >= 10) WRITE(format(19:20),'(A2)') 'I2'
		 IF (pset%sf_trunc < 10)  WRITE(format(19:20),'(A2)') 'I1'
		 IF (pset%tau >= 10)      WRITE(format(24:26),'(A2)') 'F4.1'
		 IF (pset%tau < 10)       WRITE(format(24:26),'(A2)') 'F3.1'
		 WRITE(filename, format) &
		   's-',pset%zmet,'-', pset%sf_start,"-",int(pset%sf_trunc),"-",pset%tau,pset%sf_theta,'-z'

                 CALL COMPSP(1,nz,filename,mass_ssp_zz(:,pset%zmet),lbol_ssp_zz(:,pset%zmet),&
                      spec_ssp_zz(:,:,pset%zmet),pset,ocompsp)

                 ! define redshift
                 pset%redshift_colors=0
                 pset%zred=0.0

		 WRITE(filename, format) &
		   's-',pset%zmet,'-', pset%sf_start,"-",int(pset%sf_trunc),"-",pset%tau,pset%sf_theta,'-0'
            
                 CALL COMPSP(1,nz,filename,mass_ssp_zz(:,pset%zmet),lbol_ssp_zz(:,pset%zmet),&
                      spec_ssp_zz(:,:,pset%zmet),pset,ocompsp)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END PROGRAM SIMHA
