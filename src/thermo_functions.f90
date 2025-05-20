!----------------------------------------
 function esat(t)
!****************************************
! Saturation water vapour pressure
!
! Input  : t in K
! Output : esat in Pa
!
! Jean-Francois Mahfouf (03/07/2001)
!
!****************************************
 implicit none
 real :: esat, lnesat
 real, intent(in) :: t
 if (t>=273.15) then
   esat=611.2*exp(17.67*(t-273.15)/(t-27.65)) ! bolton (1980)
 else
   lnesat=23.33086-6111.72784/t+0.15215*log(t)
   esat=100.0*exp(lnesat)                     ! emanuel (1994)
 endif
 end function esat
!---------------------------------------
 function qsat(p,t)
!****************************************
! Specific humidity at saturation
!
! Inputs : p in Pa
!          t in K
! Output : qsat in kg/kg
!
! Jean-Francois Mahfouf (03/07/2001)
!
!*****************************************
 implicit none
 real :: qsat,eps,esat
 real, intent(in) :: p,t
 real, parameter :: rv=461.5,rd=287.04
 eps=rd/rv
 qsat=eps*esat(t)/(p-esat(t)*(1.0-eps))
 end function qsat
!------------------------------------------------------------
 function desat(t)
!************************************************************
! Derivative of saturation water vapour pressure w.r.t.to T
!
! Input : t in K
! Output : desat in Pa/K
!
! Jean-Francois Mahfouf (03/07/2001)
!
!************************************************************
 implicit none
 real :: lnesat, desat
 real, intent(in) :: t
 if (t>=273.15) then
   desat=2651376.432*exp(17.67*(t-273.15)/(t-27.65))/ &
&        ((t-27.65)*(t-27.65))          
 else
   lnesat=23.33086-6111.72784/t+0.15215*log(t)
   desat=100.0*exp(lnesat)* &
&        (6111.72784/(t*t)-0.15215/t)
 endif
 end function desat 
!--------------------------------------------------------
 function dqsat(p,t)
!********************************************************
! Derivative of specific humidity at saturation w.r.t T
!
! Inputs : p in Pa
!          t in K
! Output : qsat in kg/kg/K
!
! Jean-Francois Mahfouf (03/07/2001)
!
!*********************************************************
 implicit none
 real :: dqsat,eps,desat,esat
 real, intent(in) :: p,t
 real, parameter :: rv=461.5,rd=287.04
 eps=rd/rv
 dqsat=eps*desat(t)*p/(p-esat(t)*(1.0-eps))**2
 end function dqsat
