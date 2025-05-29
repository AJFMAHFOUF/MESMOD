subroutine compute_geopotential (z_surf,pi,ts,h,phi)
 use consts
 use vgrid, only  : sigma, nz
 use hgrid, only  : nx
! 
 implicit none
!
 real, dimension (nx),      intent(in)  :: z_surf, pi, ts
 real, dimension (nx,nz),   intent(in)  :: h
 real, dimension (nx,nz+1), intent(out) :: phi
! 
 integer :: i, k
 real    :: pm, pp, phatm, phatp, theta
! 
 phi(:,nz+1) = rg*z_surf(:)
 do i=1,nx
   pm = sigma(nz)*pi(i) + ptop
   phatm = (pm/p00)**rscp
   phatp = (pi(i)/p00)**rscp
   theta = 0.5*(h(i,nz)/pi(i) + ts(i)/phatp)    
   phi(i,nz) = phi(i,nz+1) + cp*(theta*(phatp - phatm))
   phatp = phatm
   pp = pm
   do k=nz-1,1,-1
     pm = sigma(k)*pi(i) + ptop
     phatm = (pm/p00)**rscp
     theta = 0.5*(h(i,k) + h(i,k+1))/pi(i)
     phi(i,k) = phi(i,k+1) + cp*(theta*(phatp - phatm))
     phatp = phatm
     pp = pm
   enddo
 enddo 
!
 return   
end subroutine compute_geopotential
