SUBROUTINE angle_sat2toa ( earth_curv, ers2_alt, atm_hght, n_arr_pix, zen0, zen, relazm, sza, vza, aza, sca )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg, pi
  IMPLICIT NONE

  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN) :: n_arr_pix
  REAL (KIND=dp), INTENT (IN) :: earth_curv, ers2_alt, atm_hght

  ! Modified variables
  ! ==================
  REAL (KIND=dp), DIMENSION (n_arr_pix), INTENT (INOUT) :: zen0, zen, relazm

  ! Output variables
  ! ==================
  REAL (KIND=dp), INTENT (OUT)          :: sza, vza, aza, sca

  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (n_arr_pix) :: earth_zen, locazm

  ! ==================================
  ! Compute GOME viewing angles at TOA
  ! ==================================
  ! (Adopted in part from R.Spurr's geometry module)

  ! ----------------------------------
  ! Step 1: Compute Earth center angle
  ! ----------------------------------
  earth_zen = ASIN( SIN (zen) * (earth_curv+ers2_alt) / (earth_curv+atm_hght) )
  earth_zen = earth_zen - zen

  ! ------------------------------------------------------------
  ! Step 2: Correction of relative azimuth angle, save as LOCAZM
  ! ------------------------------------------------------------
  locazm = COS(earth_zen) * SIN(zen0) * COS(relazm) - SIN(ABS(earth_zen)) * COS(zen0)

  ! -----------------------------------
  ! Step 3: Correction of zenith angles
  ! -----------------------------------
  zen0 = ACOS(COS(zen0) * COS(earth_zen) + SIN(ABS(earth_zen)) * SIN(zen0) * COS(relazm) )
  zen  = zen + earth_zen

  ! Step 4: Final computation of relative azimuth angle
  ! ---------------------------------------------------
  relazm = locazm / SIN(zen0)
  WHERE ( relazm > 1.0 )    ! should be impossible
     relazm = 1.0
  ENDWHERE
  WHERE ( relazm < -1.0 )   ! should be impossible
     relazm = -1.0
  ENDWHERE
  relazm = ACOS ( relazm )
   
  !xliu: compute effective viewing geometry

  !vza: integrate from A to C over the path length (1/cos(vza))
  !use formula: ln(tan(x/2)+45)
  IF (zen(1) == zen(3)) THEN
     vza = zen(1)
  ELSE 
     vza = ACOS( 1.0 / ((LOG(TAN(zen(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen(3)/2.0+pi/4.0))) / (zen(1)-zen(3)))) 
  ENDIF

  ! sza: similar to vza
  IF (zen0(1) == zen0(3)) THEN
     sza = zen0(1)
  ELSE 
     sza = ACOS( 1.0 / ((LOG(TAN(zen0(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen0(3)/2.0+pi/4.0))) / (zen0(1)-zen0(3)))) 
  ENDIF

  ! aza: average between A and C (flip flop at B for nadir and back-scan pixels)
  aza = (180.0 - (relazm(1) + relazm(3)) / 2.0 * rad2deg) / rad2deg

  ! compute scattering angle for polarization correction
  sca = ACOS(COS(sza) * COS(vza) + SIN(sza) * SIN(vza) * COS(aza))
  sca = 180.0 - sca * rad2deg

  aza = 180.0 - aza * rad2deg 
  vza = vza * rad2deg
  sza = sza * rad2deg; 

  RETURN
END SUBROUTINE angle_sat2toa

! For OMI, viewing gometry is provided at surface
! In lidort, viewing geometry is defined at the bottom of surface. 
SUBROUTINE omi_angle_sat2toa ( n_arr_pix, zen0, zen, relazm, sza, vza, aza, sca )


  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg, pi, rearth
  USE OMSAO_variables_module,    ONLY: zatmos
  IMPLICIT NONE

  ! Input/Output variables
  ! =======================
  INTEGER,        INTENT (IN)                           :: n_arr_pix
  REAL (KIND=dp), DIMENSION (n_arr_pix), INTENT (INOUT) :: zen0, zen, relazm
  REAL (KIND=dp), INTENT (OUT)                          :: sza, vza, aza, sca

  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (n_arr_pix) :: earth_zen, locazm

  ! ==================================
  ! Compute GOME viewing angles at TOA
  ! ==================================
  ! (Adopted in part from R.Spurr's geometry module)
  zen0   = zen0   * deg2rad
  zen    = zen    * deg2rad
  relazm = relazm * deg2rad 
 
  IF (zatmos > 0.) THEN  ! No adjustment
     ! ----------------------------------
     ! Step 1: Compute Earth center angle
     ! ----------------------------------
     earth_zen = ASIN( SIN (zen) * (rearth) / (rearth + zatmos ) )
     earth_zen = earth_zen - zen
     
     ! ------------------------------------------------------------
     ! Step 2: Correction of relative azimuth angle, save as LOCAZM
     ! ------------------------------------------------------------
     locazm = COS(earth_zen) * SIN(zen0) * COS(relazm) - SIN(ABS(earth_zen)) * COS(zen0)
     
     ! -----------------------------------
     ! Step 3: Correction of zenith angles
     ! -----------------------------------
     zen0 = ACOS(COS(zen0) * COS(earth_zen) + SIN(ABS(earth_zen)) * SIN(zen0) * COS(relazm) )
     zen  = zen + earth_zen
     
     ! Step 4: Final computation of relative azimuth angle
     ! ---------------------------------------------------
     relazm = locazm / SIN(zen0)
     WHERE ( relazm > 1.0 )    ! should be impossible
        relazm = 1.0
     ENDWHERE
     WHERE ( relazm < -1.0 )   ! should be impossible
        relazm = -1.0
     ENDWHERE
     relazm = ACOS ( relazm )   
  ENDIF

  ! xliu: compute effective viewing geometry

  !vza: integrate from A to C over the path length (1/cos(vza))
  !use formula: ln(tan(x/2)+45)
  IF (zen(1) == zen(3)) THEN
     vza = zen(1)
  ELSE 
     vza = ACOS( 1.0 / ((LOG(TAN(zen(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen(3)/2.0+pi/4.0))) / (zen(1)-zen(3)))) 
  ENDIF
  !vza = ABS(zen(2))
  
  ! sza: similar to vza
  IF (zen0(1) == zen0(3)) THEN
     sza = zen0(1)
  ELSE 
     sza = ACOS( 1.0 / ((LOG(TAN(zen0(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen0(3)/2.0+pi/4.0))) / (zen0(1)-zen0(3)))) 
  ENDIF
  !sza = zen0(2)
  
  ! aza: average between A and C (flip flop at B for nadir and back-scan pixels)
  aza = (relazm(1) + relazm(3)) / 2.0
  !aza = relazm(2)

  ! compute scattering angle for polarization correction
  sca = ACOS(COS(sza) * COS(vza) + SIN(sza) * SIN(vza) * COS(aza))  * rad2deg
  sca = 180.0 - sca 

  ! Use center viewing geometry
  !sza = zen0(2)
  !vza = ABS(zen(2))
  !aza = relazm(2)
  aza = aza * rad2deg 
  vza = vza * rad2deg
  sza = sza * rad2deg

  ! This is the correct relative azmithmual angle for radiative transfer model
  aza = 180.0 - aza

  RETURN
END SUBROUTINE omi_angle_sat2toa


! Adjust the input viewing geometry from surface to those at 
! sea surface level (zero km), which is required in LIDORT 

SUBROUTINE adjust_angle (thght, zen0, zen, relazm, sca)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg, pi, rearth
  IMPLICIT NONE

  ! Input/Output variables
  ! =======================
  REAL (KIND=dp), INTENT ( IN) :: thght
  REAL (KIND=dp), INTENT (INOUT) :: zen0, zen, relazm
  REAL (KIND=dp), INTENT (  OUT) :: sca

  ! Local variables
  ! ===============
  REAL (KIND=dp)                 :: earth_zen, locazm

  zen0   = zen0   * deg2rad
  zen    = zen    * deg2rad
  relazm = relazm * deg2rad
 
  IF ( thght /= 0.0 ) THEN  ! Perform adjustment when thght is no zero.
     ! ----------------------------------
     ! Step 1: Compute Earth center angle
     ! ----------------------------------
     earth_zen = ASIN( SIN (zen) * (rearth) / (rearth - thght ) )
     earth_zen = earth_zen - zen
     
     ! ------------------------------------------------------------
     ! Step 2: Correction of relative azimuth angle, save as LOCAZM
     ! ------------------------------------------------------------
     locazm = COS(earth_zen) * SIN(zen0) * COS(relazm) - SIN(ABS(earth_zen)) * COS(zen0)
     
     ! -----------------------------------
     ! Step 3: Correction of zenith angles
     ! -----------------------------------
     zen0 = ACOS(COS(zen0) * COS(earth_zen) + SIN(ABS(earth_zen)) * SIN(zen0) * COS(relazm) )
     zen  = zen + earth_zen
     
     ! Step 4: Final computation of relative azimuth angle
     ! ---------------------------------------------------
     relazm = locazm / SIN(zen0)
     IF ( relazm >  1.0 )  relazm =  1.0    ! should be impossible
     IF ( relazm < -1.0 )  relazm = -1.0    ! should be impossible
     relazm = ACOS ( relazm )   
  ENDIF
  
  ! Recompute the scattering angle
  ! compute scattering angle for polarization correction
  relazm = pi - relazm
  sca    = 180.0 - ACOS(COS(zen0) * COS(zen) + SIN(zen0) * SIN(zen) * COS(relazm))  * rad2deg
  
  relazm = 180.0 - relazm * rad2deg 
  zen    = zen   * rad2deg
  zen0   = zen0  * rad2deg
  
  RETURN

END SUBROUTINE adjust_angle

! Calculate the probability of Sun Glint
SUBROUTINE SUNGLINT_PROBABILITY (sza, vza, azm, prob)
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  REAL (KIND=dp), INTENT (INOUT) :: sza, vza, azm
  REAL (KIND=dp), INTENT (  OUT) :: prob

  ! =======================
  ! Local Variables
  ! ======================= 
  REAL (KIND=dp) :: big, half, a1, a2, a3, b1, b2, b3, &
       c1, c2, c3, denom, n3, theta

  sza = sza * deg2rad
  vza = vza * deg2rad
  azm = (180. - azm) * deg2rad

  big  = COS(sza) * COS(vza) + SIN(sza) * SIN(vza) * COS(azm)
  half = SQRT(( big + 1.0) / 2.0)

  a1 = SIN(sza); a2 = 0.; a3 = COS(sza)
  b1 = SIN(vza) * COS(azm);  b2 = SIN(vza) * SIN(azm); b3 = COS(vza)
  c1 = -a3 * b2; c2 = a3 * b1 - a1 * b3; c3 = a1 * b2

  denom = ((b1 * a3)**2 + b2 ** 2 + (b3 * a1) ** 2 &
       - 2 * a1 * a3 * b1 * b3)
  n3    = half * (b1 * c2 - b2 * c1 - a1 * c2) / denom

  theta = ACOS(n3)
  prob  = EXP( -0.5 * (theta / (8.4932 * deg2rad)) ** 2)  
 
  sza = sza / deg2rad
  vza = vza / deg2rad
  azm = 180. - azm / deg2rad 

  RETURN
END SUBROUTINE SUNGLINT_PROBABILITY


SUBROUTINE scia_angle_sat2toa ( earth_curv, ers2_alt, atm_hght, n_arr_pix, zen0, zen, relazm, sza, vza, aza, sca )

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg, pi
  IMPLICIT NONE

  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN) :: n_arr_pix
  REAL (KIND=dp), INTENT (IN) :: earth_curv, ers2_alt, atm_hght

  ! Modified variables
  ! ==================
  REAL (KIND=dp), DIMENSION (n_arr_pix), INTENT (INOUT) :: zen0, zen, relazm

  ! Output variables
  ! ==================
  REAL (KIND=dp), INTENT (OUT)          :: sza, vza, aza, sca

  ! Local variables
  ! ===============
  REAL (KIND=dp), PARAMETER             :: toa_height = 100.0
  REAL (KIND=dp), DIMENSION (n_arr_pix) :: earth_zen, locazm

  ! ==================================
  ! Compute GOME viewing angles at TOA
  ! ==================================
  ! (Adopted in part from R.Spurr's geometry module)

  ! ----------------------------------
  ! Step 1: Compute Earth center angle
  ! ----------------------------------
  earth_zen = ASIN( SIN (zen) * (earth_curv + toa_height) / (earth_curv + atm_hght) )
  earth_zen = earth_zen - zen

  ! ------------------------------------------------------------
  ! Step 2: Correction of relative azimuth angle, save as LOCAZM
  ! ------------------------------------------------------------
  locazm = COS(earth_zen) * SIN(zen0) * COS(relazm) - SIN(ABS(earth_zen)) * COS(zen0)

  ! -----------------------------------
  ! Step 3: Correction of zenith angles
  ! -----------------------------------
  zen0 = ACOS(COS(zen0) * COS(earth_zen) + SIN(ABS(earth_zen)) * SIN(zen0) * COS(relazm) )
  zen  = zen + earth_zen

  ! Step 4: Final computation of relative azimuth angle
  ! ---------------------------------------------------
  relazm = locazm / SIN(zen0)
  WHERE ( relazm > 1.0 )    ! should be impossible
     relazm = 1.0
  ENDWHERE
  WHERE ( relazm < -1.0 )   ! should be impossible
     relazm = -1.0
  ENDWHERE
  relazm = ACOS ( relazm )

   
  !xliu: compute effective viewing geometry

  !vza: integrate from A to C over the path length (1/cos(vza))
  !use formula: ln(tan(x/2)+45)
  IF (zen(1) == zen(3)) THEN
     vza = zen(1)
  ELSE 
     vza = ACOS( 1.0 / ((LOG(TAN(zen(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen(3)/2.0+pi/4.0))) / (zen(1)-zen(3)))) 
  ENDIF

  ! sza: similar to vza
  IF (zen0(1) == zen0(3)) THEN
     sza = zen0(1)
  ELSE 
     sza = ACOS( 1.0 / ((LOG(TAN(zen0(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen0(3)/2.0+pi/4.0))) / (zen0(1)-zen0(3)))) 
  ENDIF
  
  ! aza: average between A and C (flip flop at B for nadir and back-scan pixels)
  aza = (180.0 - (relazm(1) + relazm(3)) / 2.0 * rad2deg) / rad2deg

  ! compute scattering angle for polarization correction
  sca = ACOS(COS(sza) * COS(vza) + SIN(sza) * SIN(vza) * COS(aza))
  sca = 180.0 - sca * rad2deg

  aza = 180.0 - aza * rad2deg 
  vza = vza * rad2deg
  sza = sza * rad2deg; 

  RETURN

END SUBROUTINE scia_angle_sat2toa

SUBROUTINE gome2_angle_sat2toa ( earth_curv, ers2_alt, atm_hght, n_arr_pix, zen0, zen, relazm, sza, vza, aza, sca )

  ! CRN: The angles given for GOME-2 are referenced to h=0.0 km, so corrections are not needed to put
  ! at the surface coordinate. This subroutine needs to be modified if atm_hght used is not equal to
  ! zero.


  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: deg2rad, rad2deg, pi
  IMPLICIT NONE

  ! Input variables
  ! ===============
  INTEGER,        INTENT (IN) :: n_arr_pix
  REAL (KIND=dp), INTENT (IN) :: earth_curv, ers2_alt, atm_hght

  ! Modified variables
  ! ==================
  REAL (KIND=dp), DIMENSION (n_arr_pix), INTENT (INOUT) :: zen0, zen, relazm

  ! Output variables
  ! ==================
  REAL (KIND=dp), INTENT (OUT)          :: sza, vza, aza, sca

  ! Local variables
  ! ===============
  REAL (KIND=dp), DIMENSION (n_arr_pix) :: earth_zen, locazm

  ! ==================================
  ! Compute GOME viewing angles at TOA
  ! ==================================
  ! (Adopted in part from R.Spurr's geometry module)

  ! ----------------------------------
  ! Step 1: Compute Earth center angle
  ! ----------------------------------
  !earth_zen = ASIN( SIN (zen) * (earth_curv+ers2_alt) / (earth_curv+atm_hght) )
  !earth_zen = earth_zen - zen

  ! ------------------------------------------------------------
  ! Step 2: Correction of relative azimuth angle, save as LOCAZM
  ! ------------------------------------------------------------
  !locazm = COS(earth_zen) * SIN(zen0) * COS(relazm) - SIN(ABS(earth_zen)) * COS(zen0)

  ! -----------------------------------
  ! Step 3: Correction of zenith angles
  ! -----------------------------------
!!$  zen0 = ACOS(COS(zen0) * COS(earth_zen) + SIN(ABS(earth_zen)) * SIN(zen0) * COS(relazm) )
!!$  zen  = zen + earth_zen


  ! Step 4: Final computation of relative azimuth angle
  ! ---------------------------------------------------
!!$  relazm = locazm / SIN(zen0)
!!$  WHERE ( relazm > 1.0 )    ! should be impossible
!!$     relazm = 1.0
!!$  ENDWHERE
!!$  WHERE ( relazm < -1.0 )   ! should be impossible
!!$     relazm = -1.0
!!$  ENDWHERE
!!$  relazm = ACOS ( relazm )


  !xliu: compute effective viewing geometry

  !CRN: move zen to TOA height (default for GOME-2 is 0 km)
  !zen = ASIN( earth_curv / (earth_curv + atm_hght) * SIN(pi - zen) 


  !vza: integrate from A to C over the path length (1/cos(vza))
  !use formula: ln(tan(x/2)+45)
  IF (zen(1) == zen(3)) THEN
     vza = zen(1)
  ELSE 
     vza = ACOS( 1.0 / ((LOG(TAN(zen(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen(3)/2.0+pi/4.0))) / (zen(1)-zen(3)))) 
  ENDIF

  ! sza: similar to vza
  IF (zen0(1) == zen0(3)) THEN
     sza = zen0(1)
  ELSE 
     sza = ACOS( 1.0 / ((LOG(TAN(zen0(1)/2.0+pi/4.0)) - &
          LOG(TAN(zen0(3)/2.0+pi/4.0))) / (zen0(1)-zen0(3)))) 
  ENDIF
  
  ! aza: average between A and C (flip flop at B for nadir and back-scan pixels)
  aza = (180.0 - (relazm(1) + relazm(3)) / 2.0 * rad2deg) / rad2deg

  ! compute scattering angle for polarization correction
  sca = ACOS(COS(sza) * COS(vza) + SIN(sza) * SIN(vza) * COS(aza))
  sca = 180.0 - sca * rad2deg

  aza = 180. - aza * rad2deg 
  vza = vza * rad2deg
  sza = sza * rad2deg 

  RETURN

END SUBROUTINE gome2_angle_sat2toa


