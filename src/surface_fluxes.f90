subroutine surface_fluxes (t,qv,pu,pv,phi,pi,z0,ts,qs,ustar,tstar,qstar,ff,fg)
 use consts
 use vgrid, only  : sigma, nz
 use hgrid, only  : nx, ny
!
 implicit none
!
 real, dimension (nx,ny,nz),   intent(in)  :: t, qv, pu, pv  
 real, dimension (nx,ny,nz+1), intent(in)  :: phi
 real, dimension (nx,ny),      intent(in)  :: z0, pi, ts, qs
 real, dimension (nx,ny),      intent(out) :: ustar, tstar, qstar, ff, fg
!
 integer                      :: i, j
 real                         :: pnz, thetanz, thetas, thetam, phinz, unz, cdn, cm, ch, rib
!
!  Formulation proposed by Louis (1979) 
!
!  The expressions below assume that winds have been interpolated at mass points 
!  before calling the routine 
!
 do i=1,nx
   do j=1,ny
     pnz = sigma(nz)*pi(i,j) + ptop
     thetanz = t(i,j,nz)*(1. + 0.61*qv(i,j,nz))*(p00/pnz)**rscp
     thetas = ts(i,j)*(1. + 0.61*qs(i,j))*(p00/pi(i,j))**rscp
     thetam = 0.5*(thetanz + thetas)
     phinz = phi(i,j,nz) - phi(i,j,nz+1) + rg*z0(i,j)
     unz = max(0.01,sqrt(pu(i,j,nz)*pu(i,j,nz) + pv(i,j,nz)*pv(i,j,nz))/pi(i,j))
     rib = phinz*(thetanz - thetas)/(thetam*unz*unz)
     cdn = karman/(log(phinz/(rg*z0(i,j))))
     cm = 7.4*cdn*cdn*(9.4*sqrt(phinz/(rg*z0(i,j))))
     ch = 5.3*cdn*cdn*(9.4*sqrt(phinz/(rg*z0(i,j))))
     if (rib < 0) then      
       ff(i,j) = cdn*sqrt(1.0 - 9.4*rib/(1.0 + cm*sqrt(abs(rib))))
       fg(i,j) = cdn*sqrt(1.0 - 9.4*rib/(1.0 + ch*sqrt(abs(rib))))
     else
       ff(i,j) = cdn/(1.0 + 4.7*rib)
       fg(i,j) = ff(i,j)/0.74
     endif
     ustar(i,j) = sqrt(pu(i,j,nz)*pu(i,j,nz) + pv(i,j,nz)*pv(i,j,nz))/pi(i,j)*ff(i,j)
     tstar(i,j) = (thetanz - thetas)*fg(i,j)
     qstar(i,j) = (qv(i,j,nz) - qs(i,j))*fg(i,j)
   enddo
 enddo
!
 return 
end subroutine surface_fluxes
