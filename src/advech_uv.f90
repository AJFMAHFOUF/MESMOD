subroutine advech_uv (u,v,x,d1)
 use hgrid, only : nx, ny, dx, dy
! 
 implicit none
! 
 real, dimension (nx+1,ny+1), intent(in)  :: u,v,x
 real, dimension (nx+1,ny+1), intent(out) :: d1
! 
 integer :: i, j
 real :: xpp, xpm, xmp, xmm
! 
 d1(:,:) = 0.0
! 
 do i=2,nx
   do j=2,ny
     xpp = 0.25*(x(i+1,j+1) + x(i+1,j) + x(i,j+1) + x(i,j))
     xpm = 0.25*(x(i+1,j) + x(i+1,j-1) + x(i,j) + x(i,j-1))
     xmp = 0.25*(x(i,j+1) + x(i,j) + x(i-1,j+1) + x(i-1,j))
     xmm = 0.25*(x(i,j) + x(i,j-1) + x(i-1,j) + x(i-1,j-1)) 
!
     d1(i,j) = -((u(i+1,j) + u(i,j))*(xpp + xpm) - &
   &             (u(i,j) + u(i-1,j))*(xmp + xmm))/(4.0*dx) - &
   &            ((v(i,j+1) + v(i,j))*(xpp + xmp) - &
   &             (v(i,j) + v(i,j-1))*(xpm + xmm))/(4.0*dy) 
   enddo
 enddo  
!   
 return
end subroutine advech_uv
