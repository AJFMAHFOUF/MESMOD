subroutine diff_horiz_uv (kh,x,d11)
 use hgrid, only : nx, ny, dx, dy
! 
 implicit none
! 
 real, intent(in)                         :: kh
 real, dimension (nx+1,nx+1), intent(in)  :: x
 real, dimension (nx+1,ny+1), intent(out) :: d11
! 
 integer :: i, j, im, jm, ip, jp
! 
 d11(:,:) = 0.0
! 
 do i=1,nx+1
   do j=1,ny+1
     im = max(1,i-1)
     ip = min(i+1,nx+1)
     jm = max(1,j-1)
     jp = min(j+1,ny+1) 
     d11(i,j) = kh*((x(ip,j) + x(im,j) - 2.0*x(i,j))/(dx*dx)  + &
              &     (x(i,jp) + x(i,jm) - 2.0*x(i,j))/(dy*dy))
   enddo
 enddo
! 
 return    
end subroutine diff_horiz_uv
