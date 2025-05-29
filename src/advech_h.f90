subroutine advech_h (u,h,d4)
 use hgrid, only : nx, dx
!
 implicit none
!
 real, dimension (nx+1), intent(in)   :: u
 real, dimension (nx),   intent(in)   :: h
 real, dimension (nx),   intent(out)  :: d4
!
 integer :: i
!
 d4(:) = 0.0
!
 do i=2,nx-1
   d4(i) = -(u(i+1)*(h(i+1) + h(i)) - u(i)*(h(i) + h(i-1)))/(2.0*dx) 
 enddo 
!   
 return 
end subroutine advech_h
