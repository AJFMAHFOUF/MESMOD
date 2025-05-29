subroutine pressure_gradx (k,pi,phi,tv,press,d6)
 use consts
 use vgrid, only  : sigma
 use hgrid, only  : nx, dx
! 
 implicit none
!
 integer,                intent(in)  :: k 
 real, dimension (nx),   intent(in)  :: pi,phi,tv,press
 real, dimension (nx+1), intent(out) :: d6
! 
 integer :: i, im, ip
 real, dimension (nx) :: c
 real                 :: cm
! 
 d6(:) = 0.0
 c(:) = phi(:) - rd*tv(:)*sigma(k)*pi(:)/press(:)
 do i=1,nx+1
   im = max(i-1,1)
   ip = min(i,nx)
   cm = 0.50*(c(ip) + c(im))
   d6(i) = cm*(pi(ip) - pi(im))/dx - (pi(ip)*phi(ip) - pi(im)*phi(im))/dx 
 enddo  
!
 return   
end subroutine pressure_gradx 
