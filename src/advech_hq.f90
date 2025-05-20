subroutine advech_hq (u,v,x,d4)
 use hgrid, only : nx, ny, dx, dy 
!
 implicit none
!
 real, dimension (nx+1,ny+1), intent(in)  :: u,v
 real, dimension (nx,ny), intent(in)      :: x
 real, dimension (nx,ny), intent(out)     :: d4
!
 integer :: i, j
!
 d4(:,:) = 0.0
!
 do i=2,nx-1
   do j=2,ny-1
     d4(i,j) = -((u(i+1,j+1) + u(i+1,j))*(x(i+1,j) + x(i,j)) - &
   &             (u(i,j+1) + u(i,j))*(x(i,j) + x(i-1,j)))/(4.0*dx) - &
   &            ((v(i+1,j+1) + v(i,j+1))*(x(i,j+1) + x(i,j)) - &
   &             (v(i+1,j) + v(i,j))*(x(i,j) + x(i,j-1)))/(4.0*dy) 
   enddo
 enddo 
!   
 return 
end subroutine advech_hq
