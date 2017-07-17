!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE OMO3PROF_l2_fs
! 
! This module contains the definitions for the HE5 output of OMO3PROF PGE. 
! Each geo or data field listed in this module contians the the name, 
! data type, dimeninsion names, rank, valid range, scale, offset, title, unit,
! and UniqueFieldDefinition for the field to be written out to the output file.
! To change the output geo and data field, one only need to change the 
! definition in this file. 
!
!!Revision History:
! Modified from Kai Yang's O3T_omto3_fs.f90 by Xiong Liu 
!
!!END
!!****************************************************************************

!Define data fields

!1. List of possible dimensions (will be part of swath attributes)
!   nTimes:    number of scan lines (may not be the same as that in L1B data)
!   nTimesp1:  number of scan lines + 1
!   nXtrack:   number of cross track position
!   nXtrack1:  number of cross track position + 1
!   nLayer:    number of retrieval layers
!   nLayerp1:  number of levels 
!   nFitvar:   number of fitted variables
!   mWavel:    maximum number of wavelengths (since # of wavelengths varies with x-track position)
!   nGas:      number of trace gases 
!   nAer:      number of aerosol wavelength   
!   nOth:      number of other variables
!   nCh:       number of channels

!   mTimes:    maximum number of scan lines in L1B data
!   mXtrack:   maximum number of cross-track positions in L1B data
!   xOffset:   offset in Across-Track Position
!   yOffset:   offset in Along-Track  Position
!   xBin:      number of Across-Track Coadding
!   yBin:      number of Along-Track  Coadding

!   GasNames:   name list of trace gases                       : NGas
!   Gasunits:   units of trace gases
!   OthNames:   name list of other variables                   : nOth
!   OthUnits:   units of other variables           
!   VarName:    list of retrieved variables                    :nFitvar
!   VarUnits:   units of retrieved variables                    

!2. Geolocation Fields
!   GroundPixelQualityFlags: nXtrack x nTimes
!   Latitude               : nXtrack x nTimes 
!   Longitude              : nXtrack x nTimes
!   cLatitude              : nXtrackp1 x nTimesp1
!   cLongitude             : nXtrackp1 x nTimesp1
!   SolarZenithAngle       : nXtrack x nTimes
!   ViewingZenithAngle     : nXtrack x nTimes
!   RelativeAzimuthalAngle : nXtrack x nTimes
!   Time                   : nTimes
!   SecondsInDay           : nTimes
!   SpacecraftLatitude     : nTimes
!   SpacecraftLongitude    : nTimes
!   SpacecraftAltitude     : nTimes
!   
!3. Data Fields
!   MeasurementQualityFlags: copy from level 1b data           : nTimes
!   NumberSmallPixelsColumns                                   : number of small pixels
!   PresProf:    pressure                                      : nLayerp1 x nXtrack x nTimes
!   AltProf:     altitude                                      : nLayerp1 x nXtrack x nTimes
!   TempProf:    temperature                                   : nLayerp1 x nXtrack x nTimes
  
!   OzAp:        A Priori ozone                                : nLayer x nXtrack x nTimes
!   OzApErr:     A Priori error                                : nLayer x nXtrack x nTimes
!   OzRet:       Retrieved ozone                               : nLayer x nXtrack x nTimes
!   OzNErr:      Retrieved ozone precision                     : nLayer x nXtrack x nTimes
!   OzSMErr:     Retrieved ozone solution error                : nLayer x nXtrack x nTimes

!   TOZ:         Total ozone column                            : nXtrack x nTimes     
!   TOZNErr:     Precision in total ozone column               : nXtrack x nTimes
!   TOZSnErr:    Solution Errors in total ozone column         : nXtrack x nTimes
!   SOC:         Stratospheric Ozone Column                    : nXtrack x nTimes     
!   SOCNErr:     Precision in Stratospheric ozone column       : nXtrack x nTimes
!   SOCSNErr:    Solution Errors in Stratospheric ozone column : nXtrack x nTimes
!   TOC:         Tropospheric Ozone Column                     : nXtrack x nTimes     
!   TOCNErr:     Precision in Tropospheric ozone column        : nXtrack x nTimes
!   TOCSNerr:    Solution Errors in Tropospheric ozone column  : nXtrack x nTimes
!   TropopauseIdx: index for tropopause                        : nXtrack x nTimes
!   SurfaceIdx:  Index for surface                           : nXtrack x nTimes

!   GasAp:       A Priori trace gas VCD                        : nGas x nXtrack x nTimes
!   GasApErr:    A Priori error                                : nGas
!   GasRet:      retrieved trace gases                         : nGas x nXtrack x nTimes
!   GasNerr:     retrieved trace gas precision                 : nGas x nXtrack x nTimes
!   GasSnerr:    retrieved trace gas solution error            : nGas x nXtrack x nTimes
!   GasProf:     Trace Gas Profile                             : nGas x nLayer x nXtrack x nTimes
!   GasAvgk:     Trace Gas Averaging Kernels                   : nGas x nLayer x nXtrack x nTimes

!   Other variablesm including common mode, albedo polynomial, cloud fraction polynomial, 
!   I/F shift, I/o3 shift, Ring, SO2 altitude, aerosol, surface pressure 
!   OthAp:       a priori for other variables                  : nOth x nXtrack x nTimes
!   OthApErr:    a priori error for other variables            : nOth
!   OthRet:      Retrieved for variables                       : nOth x nXtrack x nTimes
!   OthNErr:     Retrieval precision for variables             : nOth x nXtrack x nTimes
!   OthSNErr:    Retrieval solution errors for other varibles  : nOth x nXtrack x nTimes
   
!   Wavel:       Wavelength                                    : mWavel x nXtrack x nTimes
!   ComWavel:    CommonWavelength                              : mWavel
!   NumWavel:    Number of Wavelengths                         : nXtrack x nTimes
!   NumChWavel:  Number of Wavelengths in each channel         : nCh x nXtrack x nTimes
!   SimRad:      Simulated Radiances                           : mWavel x nXtrack x nTimes
!   FitSpec:     Measured Radiance                             : mWavel x nXtrack x nTimes
!   Avgk:        Averaging kernel                              : nLayer x nLayer x nXtrack x nTimes
!   Contri:      Contribution function                         : mWavel x nFitvar x nXtrack x nTimes
!   Covar:       Covariance Maxtrix                            : nLayer x nLayer x nXtrack x nTimes            
!   Correl:      Correlation Maxtrix                           : nFitvar x nFitvar x nXtrack x nTimes
!   Merr:        Meausrement Error                             : mWavel x nXtrack x nTimes
!   WeightFunc:  Weighting function                            : nFitvar x mWavel x nXtrack x nTimes
!   RingSpec:    Ring Spectrum                                 : mWavel x nXtrack x nTimes
!
!   ExitStatus:  Exit Status                                   : nXtrack x nTimes                           
!   NumIter:     number of iterations                          : nXtrack x nTimes
!   RMS:         fitting rms                                   : nCh x nXtrack x nTimes
!   AvgRes:      fitting residuals                             : nCh x nXtrack x nTimes
!   OzInfo:      Ozone Information Content                     : nXtrack x nTimes
!   CldFrac      CLoud fraction                                : nXtrack x nTimes
!   CldFlg:      Cloud flag                                    : nXtrack x nTimes
!   SfcAlb:      Surface Albedo                                : nXtrack x nTimes
!   CTP:         Cloud top pressure                            : nXtrack x nTimes
!   COD:         Cloud optical depth                           : nXtrack x nTimes    
!   GlintProb:   Sun glint probability                         : nXtrack x nTimes
!   AerWavel:    Wavelengths for inputing aerosol properties   : nAer x nXtrack x nTimes
!   AOD:         Aerosol optical thickness                     : nAer x nXtrack x nTimes
!   ASOD:        Aerosol scattering optical thickness          : nAer x nXtrack x nTimes 
!   nSaaSpike:   number of SAA large spikes                    : nXtrack x nTimes

MODULE he5_l2_fs

  USE OMSAO_he5_module
  IMPLICIT NONE

! Geolocation Fields
    TYPE (DFHE5_T) :: gf_GroundPixelQualityFlags  =                       &
           DFHE5_T( 0.0, 65534.0D0, 1.0D0, 0.0D0,                         &
                    "GroundPixelQualityFlags      ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Ground Pixel Quality Flags   ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_UINT16, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_Latitude =                                       &
           DFHE5_T( -90.0, 90.0, 1.0D0, 0.0D0,                            &
                    "Latitude                     ",                      &
                    "nXtrack,nTimes               ",                      & 
                    "deg                          ",                      &
                    "Geodetic Latitude            ",                      &
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_Longitude =                                      &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "Longitude                    ",                      &
                    "nXtrack,nTimes               ",                      & 
                    "deg                          ",                      &
                    "Geodetic Longitude           ",                      &
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_cLatitude =                                      &
           DFHE5_T( -90.0, 90.0, 1.0D0, 0.0D0,                            &
                    "cLatitude                     ",                     &
                    "nXtrackp1,nTimesp1            ",                     & 
                    "deg                           ",                     &
                    "Geodetic Corner Latitude      ",                     &
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_cLongitude =                                     &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "cLongitude                    ",                     &
                    "nXtrackp1,nTimesp1            ",                     & 
                    "deg                           ",                     &
                    "Geodetic Corner Longitude     ",                     &
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

!    TYPE (DFHE5_T) :: gf_SolarAzimuthAngle =                              &
!           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
!                    "SolarAzimuthAngle                              ",    &
!                    "nXtrack,nTimes                                 ",    & 
!                    "deg(EastofNorth)                               ",    &
!                    "Solar Azimuth Angle                            ",    &
!                    "TOMS-OMI-Shared                                ",    & 
!                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))
!
!    TYPE (DFHE5_T) :: gf_ViewingAzimuthAngle =                            &
!           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
!                    "ViewingAzimuthAngle                            ",    &
!                    "nXtrack,nTimes                                 ",    & 
!                    "deg(EastofNorth)                               ",    &
!                    "Viewing Azimuth Angle                          ",    &
!                    "TOMS-OMI-Shared                                ",    & 
!                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))
!
    TYPE (DFHE5_T) :: gf_RelativeAzimuthAngle =                           &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "RelativeAzimuthAngle                           ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg(EastofNorth)                               ",    &
                    "Relative Azimuth Angle (sun + 180 - view)      ",    &
                    "TOMS-OMI-Shared                                ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SolarZenithAngle =                               &
           DFHE5_T( 0.0, 180.0, 1.0D0, 0.0D0,                             &
                    "SolarZenithAngle                               ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg                                            ",    &
                    "Solar Zenith Angle                             ",    &
                    "TOMS-OMI-Shared                                ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_ViewingZenithAngle =                             &
           DFHE5_T( 0.0, 70.0, 1.0D0, 0.0D0,                              &
                    "ViewingZenithAngle                             ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg                                            ",    &
                    "Viewing Zenith Angle                           ",    &
                    "TOMS-OMI-Shared                                ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_Time  =                                          &
           DFHE5_T( -5.0D09, 1.0D10, 1.0D0, 0.0D0,                        & 
                    "Time                            ",                   &
                    "nTimes                          ",                   &
                    "s                               ",                   &
                    "Time at Start of Scan (TAI93)   ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_DOUBLE, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SecondsInDay =                                   &
           DFHE5_T( 0.0, 86401.0, 1.0D0, 0.0D0,                           & 
                    "SecondsInDay                    ",                   &
                    "nTimes                          ",                   &
                    "s                               ",                   &
                    "Seconds after UTC midnight      ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SpacecraftLatitude =                             &
           DFHE5_T( -90.0, 90.0, 1.0D0, 0.0D0,                            & 
                    "SpacecraftLatitude              ",                   &
                    "nTimes                          ",                   &
                    "deg                             ",                   &
                    "Spacecraft Latitude             ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SpacecraftLongitude =                            &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "SpacecraftLongitude             ",                   &
                    "nTimes                          ",                   &
                    "deg                             ",                   &
                    "Spacecraft Longitude            ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SpacecraftAltitude =                             &
           DFHE5_T( 4.0D05, 9.0D05, 1.0D0, 0.0D0,                         & 
                    "SpacecraftAltitude              ",                   &
                    "nTimes                          ",                   &
                    "m                               ",                   &
                    "Spacecraft Altitude             ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1, -1 /))

! Data Fields
    TYPE (DFHE5_T) :: df_MeasurementQualityFlags  =                       &
           DFHE5_T( 0, 254, 1.0D0, 0.0D0,                                 & 
                    "MeasurementQualityFlags         ",                   &
                    "nTimes                          ",                   &
                    "NoUnits                         ",                   &
                    "Measurement Quality Flags       ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_UINT8, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_NumberSmallPixelColumns  =                       &
           DFHE5_T( 0, 5, 1.0D0, 0.0D0,                                   & 
                    "NumberSmallPixelColumns         ",                   &
                    "nTimes                          ",                   &
                    "NoUnits                         ",                   &
                    "Number of Small Pixel Columns   ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_UINT8, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_PresProf =                                       &
           DFHE5_T( 0.01, 1100.0, 1.0D0, 0.0D0,                           & 
                    "PresProf                     ",                      &
                    "nLayerp1,nXtrack,nTimes      ",                      &
                    "hPa                          ",                      &
                    "Pressure Profile at Each Level",                     & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_AltProf =                                        &
           DFHE5_T( -1.0, 100.0, 1.0D0, 0.0D0,                            & 
                    "AltProf                      ",                      &
                    "nLayerp1,nXtrack,nTimes      ",                      &
                    "km                           ",                      &
                    "Altitude Profile at Each Level",                     & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))             

    TYPE (DFHE5_T) :: df_TempProf =                                       &
           DFHE5_T( 150, 350.0, 1.0D0, 0.0D0,                             & 
                    "TempProf                     ",                      &
                    "nLayerp1,nXtrack,nTimes      ",                      &
                    "K                            ",                      &
                    "Temperature Profile at Each Level",                  & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_OzAp =                                           &
           DFHE5_T( 0., 100.0, 1.0D0, 0.0D0,                              & 
                    "OzAp                         ",                      &
                    "nLayer,nXtrack,nTimes        ",                      &
                    "DU                           ",                      &
                    "A Priori Ozone Profile       ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OzApErr =                                        &
           DFHE5_T( 0., 100.0, 1.0D0, 0.0D0,                              & 
                    "OzApErr                      ",                      &
                    "nLayer,nXtrack,nTimes        ",                      &
                    "DU                           ",                      &
                    "A Priori Error in Ozone Profile",                    & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OzRet =                                          &
           DFHE5_T( -50., 100.0, 1.0D0, 0.0D0,                            & 
                    "OzRet                        ",                      &
                    "nLayer,nXtrack,nTimes        ",                      &
                    "DU                           ",                      &
                    "Retrieved Ozone Profile      ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OzNErr =                                         &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "OzNErr                       ",                      &
                    "nLayer,nXtrack,nTimes        ",                      &
                    "DU                           ",                      &
                    "Retrieval Precision in Ozone Profile",               & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OzSNErr =                                        &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "OzSNErr                      ",                      &
                    "nLayer,nXtrack,nTimes        ",                      &
                    "DU                           ",                      &
                    "Retrieval Solution Errors in Ozone Profile",         & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TOZ =                                            &
           DFHE5_T( 0., 600.0, 1.0D0, 0.0D0,                              & 
                    "TOZ                          ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Total Ozone Column           ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TOZNErr =                                        &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "TOZNErr                      ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Precision in Total Ozone Column",                    & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TOZSNErr =                                       &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "TOZSNErr                     ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Solution Error in Total Ozone Column",               & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SOC =                                            &
           DFHE5_T( 0., 600.0, 1.0D0, 0.0D0,                              & 
                    "SOC                          ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Stratospheric Ozone Column   ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SOCNErr =                                        &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "SOCNErr                      ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Precision in Stratospheric Ozone Column",            & 
                    "OMI-Specific                ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SOCSNErr =                                       &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "SOCSNErr                     ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Solution Error in Stratospheric Ozone Column",       & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TOC =                                            &
           DFHE5_T( -20., 150.0, 1.0D0, 0.0D0,                            & 
                    "TOC                          ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Tropspheric Ozone Column     ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TOCNErr =                                        &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "TOCNErr                      ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Precision in Tropospheric Ozone Column",             & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TOCSNErr =                                       &
           DFHE5_T( 0., 50.0, 1.0D0, 0.0D0,                               & 
                    "TOCSNErr                     ",                      &
                    "nXtrack,nTimes               ",                      &
                    "DU                           ",                      &
                    "Solution Error in Tropospheric Ozone Column",        & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_TropopauseIdx =                                  &
           DFHE5_T( 0, 100, 1.0D0, 0.0D0,                                 & 
                    "TropopauseIdx                ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Index for Tropopause         ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SurfaceIdx =                                     &
           DFHE5_T( 0, 100, 1.0D0, 0.0D0,                                 & 
                    "SurfaceIdx                   ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Index for Surface            ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasAp =                                           &
           DFHE5_T( 0.0, 1.0D20, 1.0D0, 0.0D0,                             & 
                    "GasAp                         ",                      &
                    "nGas,nXtrack,nTimes           ",                      &
                    "molecules cm^-2               ",                      &
                    "Trace Gas A Priori Vertical Column Density",          & 
                    "OMI-Specific                  ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasApErr =                                        &
           DFHE5_T( 0., 1.0D20, 1.0D0, 0.0D0,                              & 
                    "GasApErr                      ",                      &
                    "nGas,nXtrack,nTimes           ",                      &
                    "molecules cm^-2               ",                      &
                    "Trace Gas A Priori Error      ",                      & 
                    "OMI-Specific                  ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasRet =                                         &
           DFHE5_T( -1.0D20, 1.0D20, 1.0D0, 0.0D0,                        & 
                    "GasRet                        ",                     &
                    "nGas,nXtrack,nTimes           ",                     &
                    "molecules cm^-2               ",                     &
                    "Retrieved Trace Gas Vertical Column Density",        & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasNErr =                                        &
           DFHE5_T( 0.0, 1.0D20, 1.0D0, 0.0D0,                            & 
                    "GasNEr                      ",                      &
                    "nGas,nXtrack,nTimes          ",                      &
                    "molecules cm^-2              ",                      &
                    "Retrieval Precision in Trace Gas",                   & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasSNErr =                                       &
           DFHE5_T( 0.0, 1.0D20, 1.0D0, 0.0D0,                            & 
                    "GasSNErr                      ",                     &
                    "nGas,nXtrack,nTimes           ",                     &
                    "molecules cm^-2               ",                     &
                    "Retrieval Solution Errors in Trace Gas",             & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasProf =                                        &
           DFHE5_T( 0.0, 1.0D20, 1.0D0, 0.0D0,                            & 
                    "GasProf                       ",                     &
                    "nGas,nLayer,nXtrack,nTimes    ",                     &
                    "molecules cm^-2               ",                     &
                    "Trace Gas Profile             ",                     & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GasAvgk =                                        &
           DFHE5_T( -10.D0, 10.D0, 1.0D0, 0.0D0,                          & 
                    "GasAvgk                       ",                     &
                    "nGas,nLayer,nXtrack,nTimes    ",                     &
                    "NoUnits                       ",                     &
                    "Trace Gas Retrieval Averaging Kernel",               & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OthAp =                                           &
           DFHE5_T( -1.0D20, 1.0D20, 1.0D0, 0.0D0,                         & 
                    "OthAp                         ",                      &
                    "nOth,nXtrack,nTimes           ",                      &
                    "Varying                       ",                      &
                    "A Priori for Other Variables  ",                      & 
                    "OMI-Specific                  ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OthApErr =                                        &
           DFHE5_T( -1.0D20, 1.0D20, 1.0D0, 0.0D0,                         & 
                    "OthApErr                      ",                      &
                    "nOth,nXtrack,nTimes           ",                      &
                    "Varying                       ",                      &
                    "A Priori Error for Other Variables",                  & 
                    "OMI-Specific                  ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OthRet =                                          &
           DFHE5_T( -1.0D20, 1.0D20, 1.0D0, 0.0D0,                         & 
                    "OthRet                         ",                     &
                    "nOth,nXtrack,nTimes            ",                     &
                    "Varying                        ",                     &
                    "Retrievals for Other Variables ",                     & 
                    "OMI-Specific                   ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OthNErr =                                        &
           DFHE5_T( -1.0D20, 1.0D20, 1.0D0, 0.0D0,                        & 
                    "OthNerr                      ",                      &
                    "nOth,nXtrack,nTimes          ",                      &
                    "Varying                      ",                      &
                    "Retrieval Precision in Other Variables",             & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OthSNErr =                                       &
           DFHE5_T( -1.0D20, 1.0D20, 1.0D0, 0.0D0,                        & 
                    "OthSNErr                      ",                     &
                    "nOth,nXtrack,nTimes           ",                     &
                    "Varying                       ",                     &
                    "Retrieval Solution Errors in Other Variables",       & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_NumWavel  =                                      &
           DFHE5_T( 0, 1000, 1.0D0, 0.0D0,                                & 
                    "NumWavel                        ",                   &
                    "nXtrack,nTimes                  ",                   &
                    "NoUnits                         ",                   &
                    "Number of Wavelengths in the Fitting",               & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_UINT16, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_NumChWavel  =                                    &
           DFHE5_T( 0, 1000, 1.0D0, 0.0D0,                                & 
                    "NumChWavel                      ",                   &
                    "nCh,nXtrack,nTimes              ",                   &
                    "NoUnits                         ",                   &
                    "Number of Wavelengths in Each Fitting Window",       & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_UINT16, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_ComWavel  =                                      &
           DFHE5_T( 250.0, 400.0, 1.0D0, 0.0D0,                           & 
                    "ComWavel                        ",                   &
                    "mWavel                          ",                   &
                    "nm                              ",                   &
                    "Common Wavelengths Used in the Fitting",             & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_Wavel  =                                         &
           DFHE5_T( 250.0, 400.0, 1.0D0, 0.0D0,                           & 
                    "Wavel                           ",                   &
                    "mWavel,nXtrack,nTimes           ",                   &
                    "nm                              ",                   &
                    "Wavelengths Used in the Fitting ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SimRad  =                                        &
           DFHE5_T( 0.D0, 1.D0, 1.0D0, 0.0D0,                             & 
                    "SimRad                          ",                   &
                    "mWavel,nXtrack,nTimes           ",                   &
                    "I/F, NoUnits                    ",                   &
                    "Simulated Normalized Radiance   ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_FitSpec  =                                       &
           DFHE5_T( 0.D0, 1.D0, 1.0D0, 0.0D0,                             & 
                    "FitSpec                         ",                   &
                    "mWavel,nXtrack,nTimes           ",                   &
                    "I/F, NoUnits                    ",                   &
                    "Observed Normalized Radiance    ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_MErr  =                                          &
           DFHE5_T( 0.D0, 1.0D2, 1.0D0, 0.0D0,                            & 
                    "MErr                            ",                   &
                    "mWavel,nXtrack,nTimes           ",                   &
                    "relative, NoUnits               ",                   &
                    "Relative Measurement Error      ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_RingSpec  =                                      &
           DFHE5_T( -1.D0, 1.D0, 1.0D0, 0.0D0,                            & 
                    "RingSpec                        ",                   &
                    "mWavel,nXtrack,nTimes           ",                   &
                    "relative, NoUnits               ",                   &
                    "Ring Spectrum                   ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_AvgK  =                                          &
           DFHE5_T( -1.D2, 1.D2, 1.0D0, 0.0D0,                            & 
                    "AvgK                            ",                   &
                    "nLayer,nLayer,nXtrack,nTimes    ",                   &
                    "DU/DU                           ",                   &
                    "Retrieval Ozone Averaging Kernels",                  & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_Correl  =                                        &
           DFHE5_T( -1.D0, 1.D0, 1.0D0, 0.0D0,                            & 
                    "Correl                            ",                 &
                    "nFitvar,nFitvar,nXtrack,nTimes    ",                 &
                    "NoUnits                           ",                 &
                    "Correlation Matrix (all variables)",                 & 
                    "OMI-Specific                      ",                 & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_Covar  =                                         &
           DFHE5_T( -1.D20, 1.D20, 1.0D0, 0.0D0,                          & 
                    "Covar                             ",                 &
                    "nLayer,nLayer,nXtrack,nTimes      ",                 &
                    "DU                                ",                 &
                    "Covariance Matrix (ozone)         ",                 & 
                    "OMI-Specific                      ",                 & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_Contri  =                                        &
           DFHE5_T( -1.D20, 1.D20, 1.0D0, 0.0D0,                          & 
                    "Contri                            ",                 &
                    "mWavel,nFitvar,nXtrack,nTimes     ",                 &
                    "Varying                           ",                 &
                    "Contribution Function Matrix (all variables)",       & 
                    "OMI-Specific                      ",                 & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_WeightFunc  =                                    &
           DFHE5_T( -1.D29, 1.D29, 1.0D0, 0.0D0,                          & 
                    "WeightFunc                        ",                 &
                    "nFitvar,mWavel,nXtrack,nTimes     ",                 &
                    "Varying                           ",                 &
                    "Weighting Function Matrix (all variables)",          & 
                    "OMI-Specific                      ",                 & 
                     -1, HE5T_NATIVE_FLOAT, 4, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_ExitStatus  =                                    &
           DFHE5_T( -10, 110, 1.0D0, 0.0D0,                                & 
                    "ExitStatus                      ",                   &
                    "nXtrack,nTimes                  ",                   &
                    "NoUnits                         ",                   &
                    "Retrieval Exit Status           ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_INT8, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_NumIter  =                                       &
           DFHE5_T( 0, 20, 1.0D0, 0.0D0,                                  & 
                    "NumIter                         ",                   &
                    "nXtrack,nTimes                  ",                   &
                    "NoUnits                         ",                   &
                    "Number of Iterations            ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_RMS  =                                           &
           DFHE5_T( 0, 1.0D4, 1.0D0, 0.0D0,                               & 
                    "RMS                             ",                   &
                    "nCh,nXtrack,nTimes              ",                   &
                    "NoUnits                         ",                   &
                    "Ratio of Fitting Residual to Measurement Error",     & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_AvgRes  =                                        &
           DFHE5_T( 0, 1.0D4, 1.0D0, 0.0D0,                               & 
                    "AvgRes                          ",                   &
                    "nCh,nXtrack,nTimes              ",                   &
                    "%                               ",                   &
                    "Average Fitting Residuals       ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_OzInfo    =                                      &
           DFHE5_T( 0.0, 100., 1.0D0, 0.0D0,                              & 
                    "OzInfo                                   ",          &
                    "nXtrack,nTimes                           ",          &
                    "NoUnits                                  ",          &
                    "Ozone Information Content                ",          & 
                    "OMI-Specific                             ",          & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_CldFrac    =                                     &
           DFHE5_T( -0.1, 1.1, 1.0D0, 0.0D0,                              & 
                    "CldFrac                                  ",          &
                    "nXtrack,nTimes                           ",          &
                    "NoUnits                                  ",          &
                    "Effective Cloud Fraction                 ",          & 
                    "OMI-Specific                             ",          & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_CTP  =                                           &
           DFHE5_T( 90.0, 1013.25, 1.0D0, 0.0D0,                          & 
                    "CTP                          ",                      &
                    "nXtrack,nTimes               ",                      &
                    "hPa                          ",                      &
                    "Effective Cloud Top Pressure ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_COD  =                                           &
           DFHE5_T( 0.0, 1000., 1.0D0, 0.0D0,                             & 
                    "COD                          ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Effective Cloud Optical Depth",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_CldFlg  =                                        &
           DFHE5_T( 0.0, 4., 1.0D0, 0.0D0,                                & 
                    "CldFlg                       ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Cloud Flag                   ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SfcAlb  =                                        &
           DFHE5_T( -0.2, 1.2, 1.0D0, 0.0D0,                              & 
                    "SfcAlb                       ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Surface Albedo (Zeroth UV2)  ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_GlintProb  =                                     &
           DFHE5_T( 0.0, 1., 1.0D0, 0.0D0,                                & 
                    "GlintProb                    ",                      &
                    "nXtrack,nTimes               ",                      &
                    "Relative, NoUnits            ",                      &
                    "Glint Probability            ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_AerWavel  =                                      &
           DFHE5_T( 250., 400., 1.0D0, 0.0D0,                             & 
                    "AerWavel                     ",                      &
                    "nAer,nXtrack,nTimes          ",                      &
                    "nm                           ",                      &
                    "Input Aerosol Wavelengths    ",                      & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_AOD  =                                           &
           DFHE5_T( 0., 20., 1.0D0, 0.0D0,                                & 
                    "AOD                          ",                      &
                    "nAer,nXtrack,nTimes          ",                      &
                    "NoUnits                      ",                      &
                    "Input Aerosol Optical Thickness",                    & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_ASOD  =                                          &
           DFHE5_T( 0., 20., 1.0D0, 0.0D0,                                & 
                    "AOD                          ",                      &
                    "nAer,nXtrack,nTimes          ",                      &
                    "NoUnits                      ",                      &
                    "Input Aerosol Scattering Optical Thickness",         & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SAANumSpike  =                                   &
           DFHE5_T( 0.0, 100., 1.0D0, 0.0D0,                              & 
                    "SAANumSpike                  ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Number of Spikes in a SAA Region",                   & 
                    "OMI-Specific                  ",                     & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1, -1 /))

   INTEGER (KIND=4), PARAMETER       :: nGeoF = 11, nGeoF1 = 2, mDataF = 58
   INTEGER (KIND=4)                  :: nDataF
   TYPE (DFHE5_T), DIMENSION(nGeoF)  :: gf_o3prof
   TYPE (DFHE5_T), DIMENSION(nGeoF1) :: gf1_o3prof
   TYPE (DFHE5_T), DIMENSION(mDataF) :: df_o3prof
   TYPE (DFHE5_T), DIMENSION(1)      :: df_wav
   TYPE (DFHE5_T), DIMENSION(2)      :: df_l1b

   PUBLIC :: DefineFields
 CONTAINS 

   ! Define Geolocation and Data Fields
   SUBROUTINE DefineFields (nGas, nOth, wrtavgk, wrtcorr, wrtcovar, wrtcontri, wrtres, wrtwf, &
        wrtsnr, wrtring, gaswrt, wrtvar, do_lambcld, aerosol, saa_flag, reduce_resolution )

     IMPLICIT NONE
     INTEGER, INTENT (IN) :: nGas, nOth
     LOGICAL, INTENT (IN) :: wrtavgk, wrtcorr, wrtcovar, wrtcontri, wrtres, wrtwf, &
          wrtsnr, wrtring, do_lambcld, aerosol, saa_flag, reduce_resolution, gaswrt, wrtvar
     INTEGER              :: sidx, eidx
   
     gf_o3prof(1:nGeoF) = (/ &
          gf_GroundPixelQualityFlags, &       ! nXtrack,nTimes
          gf_Latitude               , &
          gf_Longitude              , &
          gf_SolarZenithAngle       , &
!          gf_SolarAzimuthAngle      , &
          gf_ViewingZenithAngle     , &
!          gf_ViewingAzimuthAngle    , &
          gf_RelativeAzimuthAngle   , &
          gf_Time                   , &       ! nTimes
          gf_SecondsInDay           , &
          gf_SpacecraftLatitude     , &
          gf_SpacecraftLongitude    , &
          gf_SpacecraftAltitude     /)

     gf1_o3prof(1:nGeoF1) = (/ &
          gf_cLatitude              , &       ! nXtrackp1 x nTimesp1
          gf_cLongitude             /)

     ! Necessary data fields
     nDataF = 31
     df_o3prof(1:31) = (/ &
          df_PresProf,                   &    ! nLayerp1 x nXtrack x nTimes
          df_AltProf,                    &
          df_TempProf,                   &
          df_OzAp,                       &    ! nLayerp1 x nXtrack x nTimes  
          df_OzApErr,                    &
          df_OzRet,                      &
          df_OzNErr,                     &
          df_OzSNErr,                    &    
          df_TOZ,                        &    ! nXtrack x nTimes
          df_TOZNErr,                    &
          df_TOZSNErr,                   &
          df_SOC,                        &
          df_SOCNErr,                    &
          df_SOCSNErr,                   &
          df_TOC,                        &
          df_TOCNErr,                    &
          df_TOCSNErr,                   &
          df_TropopauseIdx,              &
          df_SurfaceIdx,                 &
          df_ExitStatus,                 &
          df_NumIter,                    &
          df_OzInfo,                     &
          df_CldFrac,                    &
          df_CTP,                        &
          df_CldFlg,                     &
          df_SfcAlb,                     &
          df_GlintProb,                  &
          df_RMS,                        &    ! nCh x nXtrack x nTimes
          df_AvgRes,                     &
          df_NumWavel,                   &    ! nXtrack x nTimes
          df_NumChWavel                  /)   ! nCh x nXtrack x nTimes 
     
     IF ( nGas > 0 .AND. (gaswrt .OR. wrtvar)) THEN
        sidx = nDataF + 1
        eidx = nDataF + 5       ! GasAvgk/GasProf not implemented
        nDataF = eidx
        df_o3prof(sidx:eidx) = (/        &
             df_GasAp,                   &    ! nXtrack x nTimes
             df_GasApErr,                &   
             df_GasRet,                  &   
             df_GasNErr,                 &
             df_GasSNErr                /)   
             !df_GasProf                 /)   ! nLayer x nXtrack x nTimes
             !df_GasAvgk,                &        
     ENDIF
     
     IF ( nOth > 0 .AND. wrtvar) THEN
        sidx = nDataF + 1
        eidx = nDataF + 5   
        nDataF = eidx    
        df_o3prof(sidx:eidx) = (/        &    ! nXtrack x nTimes
             df_OthAp,                   &    
             df_OthApErr,                &
             df_OthRet,                  &
             df_OthNErr,                 &
             df_OthSNErr                /)
     ENDIF

     IF (wrtavgk) THEN                        ! nLayer x nLayer x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_AvgK
        nDataF = nDataF + 1
     ENDIF

     IF (wrtcorr) THEN                        ! nFitvar x nFitvar x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_Correl
        nDataF = nDataF + 1
     ENDIF

     IF (wrtcovar) THEN                       ! nLayer x nLayer x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_Covar    
        nDataF = nDataF + 1
     ENDIF

     IF (wrtcontri) THEN                      ! mWavel x nFitvar x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_Contri
        nDataF = nDataF + 1
     ENDIF

     IF (wrtwf) THEN                          ! nFitvar x mWavel x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_WeightFunc
        nDataF = nDataF + 1
     ENDIF
        
     IF ( (.NOT. reduce_resolution) .AND. &
          (wrtcontri .OR. wrtres .OR. wrtwf .OR. wrtsnr .OR. wrtring)) THEN
        nDataF = nDataF + 1
        df_o3prof(nDataF) = df_Wavel          ! mWavel x nXtrack x nTimes
     ENDIF

     IF (wrtres) THEN                         ! mWavel x nFitvar x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_SimRad
        df_o3prof(nDataF + 2) = df_FitSpec
        nDataF = nDataF + 2
     ENDIF
        
     IF (wrtsnr) THEN                         ! mWavel x nFitvar x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_MErr
        nDataF = nDataF + 1
     ENDIF

     IF (wrtring) THEN                        ! mWavel x nFitvar x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_RingSpec
        nDataF = nDataF + 1
     ENDIF

     IF (.NOT. do_lambcld) THEN               ! nXtrack x nTimes            
        df_o3prof(nDataF + 1) = df_COD
        nDataF = nDataF + 1
     ENDIF  

     IF (saa_flag) THEN                       ! nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_SAANumSpike
        nDataF = nDataF + 1
     ENDIF

     IF (aerosol) THEN                        ! nAer x nXtrack x nTimes
        df_o3prof(nDataF + 1) = df_AerWavel
        df_o3prof(nDataF + 2) = df_AOD
        df_o3prof(nDataF + 3) = df_ASOD
        nDataF = nDataF + 3
     ENDIF

     df_wav(:) = (/df_ComWavel/)

     df_l1b(1:2) = (/df_MeasurementQualityFlags, df_NumberSmallPixelColumns/)

   RETURN
 END SUBROUTINE DefineFields
   
END MODULE he5_l2_fs
