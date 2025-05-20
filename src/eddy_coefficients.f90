subroutine eddy_coefficients (phi,ustar,ff,fg,kdifh,kdifm)
 use consts
 use vgrid, only : nz
 use hgrid, only : nx, ny
! 
 implicit none
! 
 real, dimension (nx,ny,nz+1), intent(in)  :: phi
 real, dimension (nx,ny),      intent(in)  :: ustar, ff, fg
 real, dimension (nx,ny,nz+1), intent(out) :: kdifh, kdifm
! 
 integer                      :: i, j, k 
 real                         :: zi, kzi, kzh, dkzh, zalt, zalt_nz
!
!  Eddy exchange coefficients are defined at mass points
!
 zi  = 1000.0 ! top of the planetary boundary layer
 kzi = 0.0    ! exchange coefficient at top of planetary boundary layer
!
 kdifh(:,:,:) = 0.0
 kdifm(:,:,:) = 0.0
 do i=1,nx
   do j=1,ny
     zalt_nz = (phi(i,j,nz) - phi(i,j,nz+1))/rg 
     do k=1,nz
!
!  Define altitude above surface
!      
       zalt = (phi(i,j,k) - phi(i,j,nz+1))/rg
!
!  Exchange coefficient for heat and moisture (formula from O'Brien, 1970)
!
       kzh  = ustar(i,j)*fg(i,j)*zalt_nz
       dkzh = ustar(i,j)*fg(i,j)
       kdifh(i,j,k) =  kzi + ((zi - zalt)/(zi - zalt_nz))**2 * &
                  &   (kzh - kzi + (zalt - zalt_nz) * &
                  &   (dkzh + 2.0*(kzh - kzi)/(zi - zalt_nz)))
       if (zalt > zi) kdifh(i,j,k) = kzi
!
!  Exchange coefficient for momentum (formula from O'Brien, 1970)
!
       kzh  = ustar(i,j)*ff(i,j)*zalt_nz
       dkzh = ustar(i,j)*ff(i,j)
       kdifm(i,j,k) =  kzi + ((zi - zalt)/(zi - zalt_nz))**2 * &
                  &   (kzh - kzi + (zalt - zalt_nz) * &
                  &   (dkzh + 2.0*(kzh - kzi)/(zi - zalt_nz)))
       if (zalt > zi) kdifm(i,j,k) = kzi
     enddo
   enddo
 enddo
!
 return      
end subroutine eddy_coefficients
