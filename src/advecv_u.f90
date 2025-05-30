subroutine advecv_u (dnu,nudot,sigp,sigpm,pu,d2)
 use hgrid, only : nx
 use vgrid, only : nz
! 
 implicit none
! 
 real, dimension (nz)     , intent(in)  :: dnu
 real, dimension (nx,nz+1), intent(in)  :: nudot
 real, dimension (nz+1)   , intent(in)  :: sigp
 real, dimension (nz)     , intent(in)  :: sigpm
 real, dimension (nx+1,nz), intent(in)  :: pu
 real, dimension (nx+1,nz), intent(out) :: d2
! 
 integer :: i, k
 real    :: kp, km, zzpup, zzpum
! 
 d2(:,:) = 0.0
! 
 do i=2,nx
   do k=1,nz
   
       kp  = min(k+1,nz)  
       km  = max(k-1,1)
!      
       zzpup = 0.5*(pu(i,kp) + pu(i,k))
       zzpum = 0.5*(pu(i,k) + pu(i,km))
!       
       d2(i,k) = -1.0/(sigp(k)*dnu(k))*(nudot(i,k+1)*sigpm(kp)*zzpup - nudot(i,k)*sigpm(k)*zzpum) 
    
    enddo
 enddo      
!   
 return
end subroutine advecv_u
