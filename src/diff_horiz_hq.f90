subroutine diff_horiz_hq (kh,x,d10)
 use hgrid, only : nx, ny, dx, dy
! 
 implicit none
! 
 real,                    intent(in)  :: kh
 real, dimension (nx,nx), intent(in)  :: x
 real, dimension (nx,ny), intent(out) :: d10
! 
 integer :: i, j, im, jm, ip, jp
! 
 d10(:,:) = 0.0
! 
 do i=1,nx
   do j=1,ny
     im = max(1,i-1)
     ip = min(i+1,nx)
     jm = max(1,j-1)
     jp = min(j+1,ny)
     d10(i,j) = kh*((x(ip,j) + x(im,j) - 2.0*x(i,j))/(dx*dx)  + &
              &     (x(i,jp) + x(i,jm) - 2.0*x(i,j))/(dy*dy))
   enddo
 enddo
!
 return     
end subroutine diff_horiz_hq
