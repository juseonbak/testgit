  ! *******************************************************
  ! Author:  xiong liu
  ! Date  :  July 25, 2003
  ! Purpose: Variables used in various SUBROUTINE/FUNCTION
  !			 such as fcurv, scurv, find_lcorner
  ! *******************************************************

MODULE gsvd_data_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : max_fit_pts
  USE OMSAO_indices_module,    ONLY : n_max_fitpars
  USE OMSAO_variables_module,  ONLY : n_fitvar_rad, n_rad_wvl

  IMPLICIT NONE

  INTEGER                                                :: ptr_l
  REAL (KIND=dp), DIMENSION(max_fit_pts)                 :: utg
  REAL (KIND=dp), DIMENSION(n_max_fitpars)               :: gamma, alpha, xfa
  REAL (KIND=dp), DIMENSION(n_max_fitpars,n_max_fitpars) :: x
  LOGICAL                                                :: gcv

END  MODULE gsvd_data_module

 
