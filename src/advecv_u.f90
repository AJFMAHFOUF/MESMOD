subroutine advecv_u (dnu,nudot,sigp,pu,d2)
 use hgrid, only : nx
 use vgrid, only : nz
! 
 implicit none
! 
 real, dimension (nz)     , intent(in)  :: dnu
 real, dimension (nx,nz+1), intent(in)  :: nudot
 real, dimension (0:nz+1) , intent(in)  :: sigp
 real, dimension (nx+1,nz), intent(in)  :: pu
 real, dimension (nx+1,nz), intent(out) :: d2
! 
 integer :: i, k
 real :: kp, km, nudotp1, nudotp2, sigp1, sigp2, zzpup, zzpum
! 
 d2(:,:) = 0.0
! 
 do i=2,nx
   do k=2,nz-1
!
       kp  = min(k+1,nz)  
       km  = max(k-1,1)
!           
       nudotp1 = 0.50*(nudot(i,k+1) + nudot(i-1,k+1)) 
       nudotp2 = 0.50*(nudot(i,k)   + nudot(i-1,k))     
       sigp1 = 0.5*(sigp(k)+sigp(k+1))
       sigp2 = 0.5*(sigp(k)+sigp(k-1))
       zzpup = 0.5*(pu(i,kp) + pu(i,k))
       zzpum = 0.5*(pu(i,k) + pu(i,km))
!       
       d2(i,k) = -1.0/(sigp(k)*dnu(k))*(nudotp1*sigp1*zzpup - nudotp2*sigp2*zzpum) 
     
    enddo
 enddo      
!   
 return
end subroutine advecv_u
