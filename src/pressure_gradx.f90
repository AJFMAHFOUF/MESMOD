subroutine pressure_gradx (k,pi,phi,tv,press,d6)
 use consts
 use vgrid, only  : sigma
 use hgrid, only  : nx, ny, dx
! 
 implicit none
!
 integer,                     intent(in)  :: k 
 real, dimension (nx,ny),     intent(in)  :: pi,phi,tv,press
 real, dimension (nx+1,ny+1), intent(out) :: d6
! 
 integer :: i, j, im, ip, jm, jp
 real, dimension (nx,ny) :: c
 real                    :: cm
! 
 d6(:,:) = 0.0
 c(:,:) = phi(:,:) - rd*tv(:,:)*sigma(k)*pi(:,:)/press(:,:)
 do i=1,nx+1
   do j=1,ny+1
     im = max(i-1,1)
     jm = max(j-1,1)
     ip = min(i,nx)
     jp = min(j,ny)
     cm = 0.25*(c(ip,jp) + c(ip,jm) + c(im,jp) + c(im,jm))
     d6(i,j) = cm*(pi(ip,jp) - pi(im,jp) + pi(ip,jm) - pi(im,jm))/(2.0*dx) - &
              &   (pi(ip,jp)*phi(ip,jp) - pi(im,jp)*phi(im,jp) +             &
              &    pi(ip,jm)*phi(ip,jm) - pi(im,jm)*phi(im,jm))/(2.0*dx) 
   enddo
 enddo  
!
 return   
end subroutine pressure_gradx 
