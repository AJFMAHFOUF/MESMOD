subroutine advech_u (u,pu,d1)
 use hgrid, only : nx, dx
! 
 implicit none
! 
 real, dimension (nx+1), intent(in)  :: u, pu
 real, dimension (nx+1), intent(out) :: d1
! 
 integer :: i
 real :: up, um
! 
 d1(:) = 0.0
! 
 do i=2,nx
   up = 0.50*(pu(i+1) + pu(i))
   um = 0.50*(pu(i) + pu(i-1))
   d1(i) = -((u(i+1) + u(i))*up - (u(i) + u(i-1))*um)/(2.0*dx)
 enddo  
!   
 return
end subroutine advech_u
