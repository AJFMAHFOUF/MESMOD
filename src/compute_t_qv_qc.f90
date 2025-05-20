subroutine compute_t_qv_qc (h,pq,pi,qv,ql,t)
 use consts
 use vgrid, only  : sigma, nz
 use hgrid, only  : nx, ny
! 
 implicit none
!  
 real, dimension (nx,ny,nz), intent(in)   :: h,  pq 
 real, dimension (nx,ny,nz), intent(out)  :: qv, ql, t 
 real, dimension (nx,ny),    intent(in)   :: pi
! 
 real, dimension (nx,ny,nz+1) :: p  
 integer                      :: i, j, k, niter
 real                         :: dzfdt, hspi, qw, tguess, tw, xlscp, za, zb, zf, &
                             &   zzz, qsat, dqsat
!
! Compute the pressure field 
!
 do i=1,nx
   do j=1,ny
     do k=1,nz+1
       p(i,j,k) = sigma(k)*pi(i,j) + ptop
     enddo
   enddo
 enddo
!
! Compute the temperature from the enthalpy and total humidity
! 
 do i=1,nx
   do j=1,ny
     do k=1,nz
        hspi = h(i,j,k)/pi(i,j)
        za = (p00/p(i,j,k))**rscp
        t(i,j,k)=exp(hspi)/za
       !tguess = t(i,j,k) ! temperature from the previous time step
       !hspi = h(i,j,k)/pi(i,j)
       !za = (p00/p(i,j,k))**rscp
       !zb = lv*pq(i,j,k)/(cp*pi(i,j))
       !do niter=1,5
       !  zf = hspi - log(tguess*za) - zb/tguess  
       !  dzfdt = - za/tguess + zb/(tguess*tguess)
       !  tguess = tguess - zf/dzfdt
       !enddo 
       !if ( abs(t(i,j,k) - tguess) > 50.) print *,' difference t',i,j,k,tguess - t(i,j,k)
       !t(i,j,k) = min(tguess,350.)
       !if (tguess > 4000.) print*,'ca va peter',i,j,k,tguess,h(i,j,k),p(i,j,k),pq(i,j,k)/pi(i,j)*1000.0
     enddo
   enddo
 enddo
!
! Compute dew point temperature and saturated specific humidity
!
 xlscp = lv/cp
 do i=1,nx
   do j=1,ny
     do k=1,nz
       !zzz = xlscp*pq(i,j,k)/pi(i,j) + t(i,j,k)
       !tw = t(i,j,k) 
       !print *,'temperature dans compute t',t(i,j,k),p(i,j,k),i,j,k
       !do niter=1,3
       !  zf = zzz - xlscp*qsat(p(i,j,k),tw) - tw
       !  dzfdt = - xlscp*dqsat(p(i,j,k),tw) - 1.0
       !  tw = tw - zf/dzfdt
       !enddo
       !qw = qsat(p(i,j,k),tw)
       !ql(i,j,k) = max(0.0,pq(i,j,k)/pi(i,j) - qw) ! liquid water
        qv(i,j,k) = pq(i,j,k)/pi(i,j) !- ql(i,j,k)   ! water vapour
!
! Recompute entropy, temperature and water vapour when supersaturation is diagnosed
!
       !if (ql(i,j,k) > 0.0) then
       !  qv(i,j,k) = qw
       !  t(i,j,k) = tw
       !  h(i,j,k)  = pi(i,j)*(log(tw*(p00/p(i,j,k))**rscp) + lv*qw/(cp*tw))
       !endif
     enddo
   enddo
 enddo
!
 return 
end subroutine compute_t_qv_qc
