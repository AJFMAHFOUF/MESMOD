subroutine diff_horiz_h (kh,h,d10)
 use hgrid, only : nx, dx
! 
 implicit none
! 
 real,                 intent(in)  :: kh
 real, dimension (nx), intent(in)  :: h
 real, dimension (nx), intent(out) :: d10
! 
 integer :: i, im, ip
! 
 d10(:) = 0.0
! 
 do i=1,nx
   im = max(1,i-1)
   ip = min(i+1,nx)
   d10(i) = kh*(h(ip) + h(im) - 2.0*h(i))/(dx*dx)  
 enddo
!
 return     
end subroutine diff_horiz_h
