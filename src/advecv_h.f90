subroutine advecv_h (dnu,nudot,sigp,sigpm,h,d2)
 use hgrid, only : nx
 use vgrid, only : nz
! 
 implicit none
! 
 real, dimension (nz)     , intent(in)  :: dnu
 real, dimension (nx,nz+1), intent(in)  :: nudot
 real, dimension (nz+1)   , intent(in)  :: sigp
 real, dimension (nz)     , intent(in)  :: sigpm
 real, dimension (nx,nz)  , intent(in)  :: h
 real, dimension (nx+1,nz), intent(out) :: d2
! 
 integer :: i, k
 real    :: kp, km, zzpup, zzpum  
! 
 d2(:,:) = 0.0
!  
  do i=2,nx
   do k=2,nz
   
       kp  = min(k+1,nz)  
       km  = max(k-1,1)
!       
       zzpup = 0.5*(h(i,kp) + h(i,k))
       zzpum = 0.5*(h(i,k) + h(i,km))
!       
       d2(i,k) = -1.0/(sigp(k)*dnu(k))*(nudot(i,k+1)*sigpm(kp)*zzpup - nudot(i,k)*sigpm(k)*zzpum) 
    
    enddo
 enddo    
!   
 return
end subroutine advecv_h
