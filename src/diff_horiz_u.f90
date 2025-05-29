subroutine diff_horiz_u (kh,u,d11)
 use hgrid, only : nx, dx
! 
 implicit none
! 
 real,                   intent(in)  :: kh
 real, dimension (nx+1), intent(in)  :: u
 real, dimension (nx+1), intent(out) :: d11
! 
 integer :: i, im, ip
! 
 d11(:) = 0.0
! 
 do i=1,nx+1
   im = max(1,i-1)
   ip = min(i+1,nx+1)
   d11(i) = kh*(u(ip) + u(im) - 2.0*u(i))/(dx*dx)  
 enddo
! 
 return    
end subroutine diff_horiz_u
