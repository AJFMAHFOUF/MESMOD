subroutine pressure_grady (k,pi,phi,tv,press,d7)
 use consts
 use vgrid, only  : sigma
 use hgrid, only  : nx, ny, dy 
! 
 implicit none
! 
 integer,                     intent(in)  :: k 
 real, dimension (nx,ny),     intent(in)  :: pi,phi,tv,press
 real, dimension (nx+1,ny+1), intent(out) :: d7
! 
 integer :: i, j, im, ip, jm, jp
 real, dimension (nx,ny) :: c
 real                    :: cm
!
 d7(:,:) = 0.0
 c(:,:) = phi(:,:) - rd*tv(:,:)*sigma(k)*pi(:,:)/press(:,:)
! 
 do i=1,nx+1
   do j=1,ny+1
     im = max(i-1,1)
     jm = max(j-1,1)
     ip = min(i,nx)
     jp = min(j,ny)
     cm = 0.25*(c(ip,jp) + c(ip,jm) + c(im,jp) + c(im,jm))
     d7(i,j) = cm*(pi(ip,jp) - pi(ip,jm) + pi(im,jp) - pi(im,jm))/(2.0*dy) - &
              &   (pi(ip,jp)*phi(ip,jp) - pi(ip,jm)*phi(ip,jm) +             &
              &    pi(im,jp)*phi(im,jp) - pi(im,jm)*phi(im,jm))/(2.0*dy)
   enddo
 enddo   
 return 
end subroutine pressure_grady 
