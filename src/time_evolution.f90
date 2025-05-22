subroutine time_evolution
 use consts
 use vgrid , only  : sigma, nu, nz
 use hgrid , only  : nx, ny, dx, dy
! 
 implicit none
! 
 real, dimension (nx,ny,nz)       :: pu_init, pv_init, h_init, pq_init, t, qv, ql
 real, dimension (nx,ny)          :: pi_init
 real, dimension (nx,ny,nz+1)     :: phi, kdifh, kdifm, nudotsigpi, nudot
 real, dimension (nx+1,ny+1,nz,3) :: pu, pv
 real, dimension (nx,ny,nz,3)     :: h, pq
 real, dimension (nx,ny,3)        :: pi
 real, dimension (nx,ny)          :: z_surf, z0, ts, qs, dpidt
 real, dimension (nx,ny)          :: ustar, tstar, qstar, ff, fg
 real, dimension (nx,ny)          :: phik, tv, press
 real, dimension (nx,ny)          :: fluxu, fluxv, fluxh, fluxq
 real, dimension (nx+1,ny+1)      :: fluxum, fluxvm, pim
 real, dimension (nx+1,ny+1)      :: zpu, zpv, pgradx, pgrady
 real, dimension (nx+1,ny+1)      :: advhu, advhv, diffhu, diffhv
 real, dimension (nx,ny)          :: advhh, advhq, diffhh, diffhq
 real, dimension (0:nz+1)         :: sigp
 real, dimension (nz)             :: kh
 integer                          :: i, j, k, nstep_max, im, ip, jm, jp, kp, km, nstep, nt, nti, ntf, &
                                  &  ntm, kk 
 real                             :: advvh, advvq, advvu, advvv, alpha, diffvh, diffvq, diffvu, diffvv, &
                                  &  divu, dnu, dnum, dnup, dt, f_coriolis, phatm, &
                                  &  phatp, pm, pp, qvm, theta, sigp1, sigp2, wk, za, za1, &
                                  &  za2, za3, zabar, zakhp, zakhm, zakmm, zakmp, zam, zam1, zam2, zam3, &
                                  &  zambar, zap, zap1, zap2, zap3, zapbar, zkmbar, zkmbarm, zkmbarp, zp, &
                                  &  zp1, zp2, zp3, zpm, zpm1, zpm2, zpm3, zpp, zpp1, zpp2, zpp3, zzhm, zzhp, &
                                  &  zzpqm, zzpqp, zzpum, zzpup, zzpvm, zzpvp, p, xpi
 real                             :: kscale, filter
!
! Define input - output files
!
 open (unit=88,file='../data/windsurf.dat')
 open (unit=40,file='../data/initial_conditions_ideal1.dat')
!
! Set-up various parameters
!
 xpi = acos(-1.0)
 f_coriolis = 1.0e-4
 alpha = 0.5
 wk = 0.0 ! pure Asselin filter
 kscale = 0.5*((dx*dx) + (dy*dy))/(2.*dt0)
!
! Define horizontal diffusion coefficient + sponge layer @ model top
!
 do k=1,nz
   if (k < 5) then
     kh(k) = (5.0e-3 + 0.01*(sin(0.5*xpi*(sigma(6) - sigma(k+1))/sigma(6)))**2)*kscale
   else
     kh(k) =  5.0e-3*kscale
   endif 
 enddo
!
! Number of time steps
!
 nstep_max = 720
!
! vertical coordinate parameter : dsigma/dnu 
!
 sigp(:) = 4.0*(1. - nu(:)**3)/3.0
!
! read initial conditions
!
 read (40,*) h_init
 read (40,*) pu_init
 read (40,*) pv_init
 read (40,*) pq_init
 read (40,*) t
 read (40,*) phi
 read (40,*) pi_init
 read (40,*) z0
 read (40,*) ts
 read (40,*) qs
 read (40,*) z_surf
! 
 close (unit=40)
!
 do i=1,nx
   do j=1,ny
     do nt=1,3
       do k=1,nz
         h(i,j,k,nt)  = h_init(i,j,k)
         pq(i,j,k,nt) = pq_init(i,j,k)
         pu(i,j,k,nt) = pu_init(i,j,k)
         pv(i,j,k,nt) = pv_init(i,j,k)
       enddo
     pi(i,j,nt) = pi_init(i,j)
     enddo
   enddo
 enddo
!
! Set-up boundary values for winds (add one line and one column)
!
 do j=1,ny+1
   pu(nx+1,j,:,:) = pu(nx,j,:,:)
   pv(nx+1,j,:,:) = pv(nx,j,:,:)
 enddo
 do i=1,nx+1
   pu(i,ny+1,:,:) = pu(i,ny,:,:)
   pv(i,ny+1,:,:) = pv(i,ny,:,:)
 enddo
!
 print *,'Initial fields have been read and put in appropriate arrays'
!
! Compute surface fluxes (at mass points) assuming that pu_init and pv_init are at mass points (see below the other call)
!
 do k=1,nz
   qv(:,:,k) = pq_init(:,:,k)/pi_init(:,:)
 enddo
!
 call surface_fluxes (t,qv,pu_init,pv_init,phi,pi_init,z0,ts,qs,ustar,tstar,qstar,ff,fg)
!
 fluxu(:,:) = pi(:,:,1)*ustar(:,:)*ustar(:,:)*pu_init(:,:,nz)/(sqrt(pu_init(:,:,nz)**2 + pv_init(:,:,nz)**2))
 fluxv(:,:) = pi(:,:,1)*ustar(:,:)*ustar(:,:)*pv_init(:,:,nz)/(sqrt(pu_init(:,:,nz)**2 + pv_init(:,:,nz)**2))
 do i=1,nx+1
   do j=1,ny+1
     im = max(1,i-1)
     jm = max(1,j-1)
     ip = min(i,nx)
     jp = min(j,ny) 
     fluxum(i,j) = 0.25*(fluxu(ip,jp) + fluxu(ip,jm) + fluxu(im,jp) + fluxu(im,jm))
     fluxvm(i,j) = 0.25*(fluxv(ip,jp) + fluxv(ip,jm) + fluxv(im,jp) + fluxv(im,jm))
   enddo
 enddo
 fluxh(:,:) = pi(:,:,1)*(ustar(:,:)*tstar(:,:)*((pi(:,:,1)*sigma(nz) + ptop)/p00)**rscp/t(:,:,nz) + &
            &         lv*ustar(:,:)*qstar(:,:)/(cp*t(:,:,nz)))
 fluxq(:,:) = pi(:,:,1)*ustar(:,:)*qstar(:,:)
! 
print *,'Initial surface fluxes ok',fluxu(15,14),fluxv(15,14),fluxh(15,14),fluxq(15,14)
!
! Compute exchange coefficients (at mass points)
!
 call eddy_coefficients (phi,ustar,ff,fg,kdifh,kdifm)
 print *,'Exchange coefficients ok'
!
! Set-up vertical velocity to zero
!
 nudotsigpi(:,:,:) = 0.0
 dpidt(:,:) = 0.0
 nudot(:,:,:) = 0.0
!
! Start temporal loop
!
 do nstep=1,nstep_max
   print  *,'nstep',nstep
!
! Horizontal transfers (advection - pressure gradient - horizontal diffusion)
!
   if (nstep == 1) then
     nti = 1
     ntm = 1
     ntf = 2
     dt = dt0
   else
     nti = 1
     ntm = 2
     ntf = 3
     dt = 2.*dt0
   endif
!
   do k=1,nz ! Start vertical loop 
!
! Horizontal advection: leapfrog scheme
!
     do i=1,nx+1
       do j=1,ny+1
         im = max(1,i-1)
         jm = max(1,j-1)
         ip = min(i,nx)
         jp = min(j,ny)
         pim(i,j) = 0.25*(pi(ip,jp,ntm) + pi(im,jp,ntm) + pi(ip,jm,ntm) + pi(im,jm,ntm))
       enddo
     enddo 
     zpu(:,:) = pu(:,:,k,ntm)/pim(:,:)
     zpv(:,:) = pv(:,:,k,ntm)/pim(:,:)

     call advech_uv (zpu,zpv,pu(:,:,k,ntm),advhu)
     call advech_uv (zpu,zpv,pv(:,:,k,ntm),advhv)
     call advech_hq (zpu,zpv,h (:,:,k,ntm),advhh)
     call advech_hq (zpu,zpv,pq(:,:,k,ntm),advhq)
!
! Horizontal diffusion: Euler forward scheme 
!    
!     do i=1,nx+1
!       do j=1,ny+1
!         im = max(1,i-1)
!         jm = max(1,j-1)
!         ip = min(i,nx)
!         jp = min(j,ny)
!         pim(i,j) = 0.25*(pi(ip,jp,nti) + pi(im,jp,nti) + pi(ip,jm,nti) + pi(im,jm,nti))
!       enddo
!     enddo 
!     zpu(:,:) = pu(:,:,k,nti)/pim(:,:)
!     zpv(:,:) = pv(:,:,k,nti)/pim(:,:)

     call diff_horiz_uv (kh(k),pu(:,:,k,nti),diffhu)
     call diff_horiz_uv (kh(k),pv(:,:,k,nti),diffhv)
     call diff_horiz_hq (kh(k),h (:,:,k,nti),diffhh)
     call diff_horiz_hq (kh(k),pq(:,:,k,nti),diffhq)
! 
     phik(:,:) = phi(:,:,k)
     tv(:,:) = t(:,:,k)*(1.0 + 0.608*qv(:,:,k))
     press(:,:) = sigma(k)*pi(:,:,ntm) + ptop
!
     call pressure_gradx (k,pi(:,:,ntm),phik,tv,press,pgradx)
     call pressure_grady (k,pi(:,:,ntm),phik,tv,press,pgrady)
!
     pu(:,:,k,ntf) = pu(:,:,k,nti) + dt*diffhu(:,:)  ! + dt*(advhu(:,:) + f_coriolis*pv(:,:,k,ntm) + pgradx(:,:) + pim(:,:)*diffhu(:,:))
     pv(:,:,k,ntf) = pv(:,:,k,nti) + dt*diffhv(:,:) ! + dt*(advhv(:,:) - f_coriolis*pu(:,:,k,ntm) + pgrady(:,:) + pim(:,:)*diffhv(:,:)) 
     h (:,:,k,ntf) = h (:,:,k,nti) + dt*diffhh(:,:) ! + dt*(advhh(:,:) + diffhh(:,:))
     pq(:,:,k,ntf) = pq(:,:,k,nti) + dt*diffhq(:,:) ! + dt*(advhq(:,:) + diffhq(:,:)) 
!     
   enddo ! End vertical loop
!
!  Vertical transfers (vertical advection - vertical diffusion)
!
   dpidt (:,:) = 0.0
   do i=2,nx
     do j=2,ny
       do k=1,nz
!
         kp  = min(k+1,nz)
         km  = max(k-1,1)
!
         zp    = sigma(k)*pi(i,j,ntm) + ptop
         zp1   = sigma(k)*pi(i-1,j,ntm) + ptop
         zp2   = sigma(k)*pi(i,j-1,ntm) + ptop
         zp3   = sigma(k)*pi(i-1,j-1,ntm) + ptop
!
         zpp   = sigma(kp)*pi(i,j,ntm) + ptop
         zpp1  = sigma(kp)*pi(i-1,j,ntm) + ptop
         zpp2  = sigma(kp)*pi(i,j-1,ntm) + ptop
         zpp3  = sigma(kp)*pi(i-1,j-1,ntm) + ptop
!
         zpm   = sigma(k-1)*pi(i,j,ntm) + ptop
         zpm1  = sigma(k-1)*pi(i-1,j,ntm) + ptop
         zpm2  = sigma(k-1)*pi(i,j-1,ntm) + ptop
         zpm3  = sigma(k-1)*pi(i-1,j-1,ntm) + ptop
!
         za  = -rg*zp /(rd*t(i,j,k)*pi(i,j,ntm)*sigp(k))
         za1 = -rg*zp1/(rd*t(i-1,j,k)*pi(i-1,j,ntm)*sigp(k))
         za2 = -rg*zp2/(rd*t(i,j-1,k)*pi(i,j-1,ntm)*sigp(k))
         za3 = -rg*zp3/(rd*t(i-1,j-1,k)*pi(i-1,j-1,ntm)*sigp(k))
!
         zap  = -rg*zpp /(rd*t(i,j,kp)*pi(i,j,ntm)*sigp(kp))
         zap1 = -rg*zpp1/(rd*t(i-1,j,kp)*pi(i-1,j,ntm)*sigp(kp))
         zap2 = -rg*zpp2/(rd*t(i,j-1,kp)*pi(i,j-1,ntm)*sigp(kp))
         zap3 = -rg*zpp3/(rd*t(i-1,j-1,kp)*pi(i-1,j-1,ntm)*sigp(kp))

         zam  = -rg*zpm /(rd*t(i,j,km)*pi(i,j,ntm)*sigp(k-1))
         zam1 = -rg*zpm1/(rd*t(i-1,j,km)*pi(i-1,j,ntm)*sigp(k-1))
         zam2 = -rg*zpm2/(rd*t(i,j-1,km)*pi(i,j-1,ntm)*sigp(k-1))
         zam3 = -rg*zpm3/(rd*t(i-1,j-1,km)*pi(i-1,j-1,ntm)*sigp(k-1))
!
         zabar  = 0.25*(za  + za1  + za2  + za3 )
         zapbar = 0.25*(zap + zap1 + zap2 + zap3)
         zambar = 0.25*(zam + zam1 + zam2 + zam3)
!
         zkmbar  = 0.25*(kdifm(i,j,k)  + kdifm(i-1,j,k)  + kdifm(i,j-1,k)  + kdifm(i-1,j-1,k) )
         zkmbarp = 0.25*(kdifm(i,j,kp) + kdifm(i-1,j,kp) + kdifm(i,j-1,kp) + kdifm(i-1,j-1,kp))
         zkmbarm = 0.25*(kdifm(i,j,km) + kdifm(i-1,j,km) + kdifm(i,j-1,km) + kdifm(i-1,j-1,km))
!
!   Vertical advection ------------------------------------------------------------------
!
         if (k == 1) then
           dnu = 0.5*(nu(k) + nu(k+1)) - nu(k-1)
         elseif (k == nz) then
           dnu = nu(k+1) - 0.5*(nu(k)+nu(k-1))
         else
           dnu = 0.5*(nu(k+1) - nu(k-1))
         endif
!  
         sigp1 = 0.5*(sigp(k)+sigp(k+1))
         sigp2 = 0.5*(sigp(k)+sigp(k-1))
         zzpup = 0.5*(pu(i,j,kp,ntm) + pu(i,j,k,ntm))
         zzpum = 0.5*(pu(i,j,k,ntm) + pu(i,j,km,ntm))
         advvu = -1.0/(sigp(k)*dnu)*(nudot(i,j,k+1)*sigp1*zzpup - nudot(i,j,k)*sigp2*zzpum) 
         zzpvp = 0.5*(pv(i,j,kp,ntm) + pv(i,j,k,ntm))
         zzpvm = 0.5*(pv(i,j,k,ntm) + pv(i,j,km,ntm))
         advvv = -1.0/(sigp(k)*dnu)*(nudot(i,j,k+1)*sigp1*zzpvp - nudot(i,j,k)*sigp2*zzpvm) 
         zzhp  = 0.5*(h(i,j,kp,ntm) + h(i,j,k,ntm))
         zzhm  = 0.5*(h(i,j,k,ntm) + h(i,j,km,ntm))
         advvh = -1.0/(sigp(k)*dnu)*(nudot(i,j,k+1)*sigp1*zzhp - nudot(i,j,k)*sigp2*zzhm) 
         zzpqp = 0.5*(pq(i,j,kp,ntm) + pq(i,j,k,ntm))
         zzpqm = 0.5*(pq(i,j,k,ntm) + pq(i,j,km,ntm))
         advvq = -1.0/(sigp(k)*dnu)*(nudot(i,j,k+1)*sigp1*zzpqp - nudot(i,j,k)*sigp2*zzpqm) 
!
!   Vertical diffusion --------------------------------------------------------------------
!
         zakmp = 0.5*(zapbar*zkmbarp + zabar*zkmbar)
         zakmm = 0.5*(zambar*zkmbarm + zabar*zkmbar)
!
         zakhp = 0.5*(zap*kdifh(i,j,kp) + za*kdifh(i,j,k))
         zakhm = 0.5*(zam*kdifh(i,j,km) + za*kdifh(i,j,k))
!
         dnup = nu(k+1) - nu(k)
         dnum = nu(k) - nu(k-1)

         if (k /= nz) then 
     
           diffvu = (zabar/dnu)*((pu(i,j,kp,nti) - pu(i,j,k,nti))*zakmp/dnup - & 
                  &              (pu(i,j,k,nti) - pu(i,j,km,nti))*zakmm/dnum)
           diffvv = (zabar/dnu)*((pv(i,j,kp,nti) - pv(i,j,k,nti))*zakmp/dnup - & 
                  &              (pv(i,j,k,nti) - pv(i,j,km,nti))*zakmm/dnum)
           diffvh = (za/dnu)*((h (i,j,kp,nti) - h (i,j,k,nti))*zakhp/dnup - &  
                  &           (h (i,j,k,nti) - h (i,j,km,nti))*zakhm/dnum)
           diffvq = (za/dnu)*((pq(i,j,kp,nti) - pq(i,j,k,nti))*zakhp/dnup - &  
                  &           (pq(i,j,k,nti) - pq(i,j,km,nti))*zakhm/dnum)

         else ! introduction of surface fluxes
     
           diffvu = (zabar/dnu)*(fluxum(i,j) - & 
                  &             (pu(i,j,k,nti) - pu(i,j,km,nti))*zakmm/dnum)
           diffvv = (zabar/dnu)*(fluxvm(i,j) - & 
                  &             (pv(i,j,k,nti) - pv(i,j,km,nti))*zakmm/dnum)
           diffvh = (za/dnu)*(fluxh(i,j) - &  
                  &          (h (i,j,k,nti) - h (i,j,km,nti))*zakhm/dnum)
           diffvq = (za/dnu)*(fluxq(i,j) - &  
                  &          (pq(i,j,k,nti) - pq(i,j,km,nti))*zakhm/dnum)
         endif
!
         pu(i,j,k,ntf) =  pu(i,j,k,ntf) + dt*(diffvu + advvu)!(advvu + diffvu)
         pv(i,j,k,ntf) =  pv(i,j,k,ntf) + dt*(diffvv + advvv)!(advvv + diffvv)
         h (i,j,k,ntf) =  h (i,j,k,ntf) + dt*(diffvh + advvh)!(advvh + diffvh)
         pq(i,j,k,ntf) =  pq(i,j,k,ntf) + dt*(diffvq + advvq)!(advvq + diffvq)
!
         divu = (pu(i+1,j+1,k,ntm) - pu(i,j+1,k,ntm) + pu(i+1,j,k,ntm) - pu(i,j,k,ntm))/(2.0*dx) + &
             &  (pv(i+1,j+1,k,ntm) - pv(i+1,j,k,ntm) + pv(i,j+1,k,ntm) - pv(i,j,k,ntm))/(2.0*dy)
         dpidt(i,j) =  dpidt(i,j) - dnu*sigp(k)*divu
       enddo ! end of k loop
       pi(i,j,ntf) = pi(i,j,nti) + dt*dpidt(i,j)
     enddo ! end of j loop
   enddo ! end of i loop
!
! Write wind field in file
!
   if (nstep == 720) then
     do i = 1,nx
       do j = 1,ny
         kk = 11
         write (88,*) i,j,(pu(i,j,kk,ntm)/pim(i,j)),(pv(i,j,kk,ntm)/pim(i,j)) !,exp(h(i,j,kk,ntm)/pi(i,j,ntm)),ustar(i,j),cp*ustar(i,j)*tstar(i,j)
       enddo
     enddo
     close (unit=88)
   endif 
!
! Compute vertical velocity
!
   nudotsigpi(:,:,:) = 0.0
   do i=1,nx
     do j=1,ny
       do k=1,nz
         do kk=1,k
           if (kk == 1) then
             dnu = 0.5*(nu(kk) + nu(kk+1)) - nu(kk-1)
           elseif (kk == nz) then
             dnu = nu(kk+1) - 0.5*(nu(kk)+nu(kk-1))
           else
             dnu = 0.5*(nu(kk+1) - nu(kk-1))
           endif
           divu = (pu(i+1,j+1,kk,ntm) - pu(i,j+1,kk,ntm) + pu(i+1,j,kk,ntm) - pu(i,j,kk,ntm))/(2.0*dx) + &
               &  (pv(i+1,j+1,kk,ntm) - pv(i+1,j,kk,ntm) + pv(i,j+1,kk,ntm) - pv(i,j,kk,ntm))/(2.0*dy)
           nudotsigpi(i,j,k) = nudotsigpi(i,j,k) - dnu*sigp(kk)*(dpidt(i,j) + divu)      
         enddo
       enddo
     enddo
   enddo
!
   nudot(:,:,:) = 0.0
   do i=1,nx
     do j=1,ny
       do k=2,nz
         nudot(i,j,k) = nudotsigpi(i,j,k)/(pi(i,j,ntm)*sigp(k))
       enddo
     enddo
   enddo 
!
!  Apply an Asselin filter on prognostic variables (modified by Williams, 2009)
!
   if (nstep > 1) then
     do i=1,nx+1
       do j=1,ny+1
         do k=1,nz
           filter = pu(i,j,k,3) + pu(i,j,k,1) - 2.0*pu(i,j,k,2)
           pu(i,j,k,2) = pu(i,j,k,2) + 0.5*alpha*wk*filter
           pu(i,j,k,3) = pu(i,j,k,3) - 0.5*alpha*(1.0 - wk)*filter
           filter = pv(i,j,k,3) + pv(i,j,k,1) - 2.0*pv(i,j,k,2)
           pv(i,j,k,2) = pv(i,j,k,2) + 0.5*alpha*wk*filter
           pv(i,j,k,3) = pv(i,j,k,3) - 0.5*alpha*(1.0 - wk)*filter            
         enddo
       enddo
     enddo  
     do i=1,nx
       do j=1,ny
         do k=1,nz
           filter = h(i,j,k,3) + h(i,j,k,1) - 2.0*h(i,j,k,2)
           h(i,j,k,2) = h(i,j,k,2) + 0.5*alpha*wk*filter
           h(i,j,k,3) = h(i,j,k,3) - 0.5*alpha*(1.0 - wk)*filter
           filter = pq(i,j,k,3) + pq(i,j,k,1) - 2.0*pq(i,j,k,2)
           pq(i,j,k,2) = pq(i,j,k,2) + 0.5*alpha*wk*filter
           pq(i,j,k,3) = pq(i,j,k,3) - 0.5*alpha*(1.0 - wk)*filter            
         enddo
         filter = pi(i,j,3) + pi(i,j,1) - 2.0*pi(i,j,2)
         pi(i,j,2) = pi(i,j,2) + 0.5*alpha*wk*filter
         pi(i,j,3) = pi(i,j,3) - 0.5*alpha*(1.0 - wk)*filter            
       enddo
     enddo  
   endif
!
!  Swapp variables for next model time step
!
   if (nstep > 1) then
     pu(:,:,:,1) = pu(:,:,:,2)
     pu(:,:,:,2) = pu(:,:,:,3)
     pv(:,:,:,1) = pv(:,:,:,2) 
     pv(:,:,:,2) = pv(:,:,:,3)
     h (:,:,:,1) = h (:,:,:,2) 
     h (:,:,:,2) = h (:,:,:,3)
     pq(:,:,:,1) = pq(:,:,:,2) 
     pq(:,:,:,2) = pq(:,:,:,3)
     pi(:,:,1)   = pi(:,:,2) 
     pi(:,:,2)   = pi(:,:,3)
   endif
!
!  Define upper boundary conditions 
!
   pu(:,:,1,:) = pu(:,:,2,:)
   pv(:,:,1,:) = pv(:,:,2,:)
!  h(:,:,1,:)  = h(:,:,2,:)
   pq(:,:,1,:) = pq(:,:,2,:)
!
!  Define lateral boundary conditions (thermodynamical variables)
!
   h(nx,:,:,:)  = h(nx-1,:,:,:)
   pq(nx,:,:,:) = pq(nx-1,:,:,:)
   h(1,:,:,:)   = h(2,:,:,:)
   pq(1,:,:,:)  = pq(2,:,:,:)
   h(:,ny,:,:)  = h(:,ny-1,:,:)
   pq(:,ny,:,:) = pq(:,ny-1,:,:)
   h(:,1,:,:)   = h(:,2,:,:)
   pq(:,1,:,:)  = pq(:,2,:,:)
   pi(:,1,:)    = pi(:,2,:)
   pi(1,:,:)    = pi(2,:,:)
   pi(nx,:,:)   = pi(nx-1,:,:)
   pi(:,ny,:)   = pi(:,ny-1,:)
!
!  Corner values
!
   h(nx,1,:,:)  = h(nx-1,2,:,:)
   h(1,ny,:,:)  = h(2,ny-1,:,:)
   h(1,1,:,:)   = h(2,2,:,:)
   h(nx,ny,:,:) = h(nx-1,ny-1,:,:)
!
   pq(nx,1,:,:)  = pq(nx-1,2,:,:)
   pq(1,ny,:,:)  = pq(2,ny-1,:,:)
   pq(1,1,:,:)   = pq(2,2,:,:)
   pq(nx,ny,:,:) = pq(nx-1,ny-1,:,:)
!
   pi(nx,1,:)   = pi(nx-1,2,:)
   pi(1,ny,:)   = pi(2,ny-1,:)
   pi(1,1,:)    = pi(2,2,:)
   pi(nx,ny,:)  = pi(nx-1,ny-1,:)
!
!  Lower horizontal axis (wind components)
!
   do k=1,nz
     do i=1,nx+1
       if (pv(i,2,k,2) <= 0.0) then ! outflow
         pu(i,1,k,2) = pu(i,2,k,2)
         pv(i,1,k,2) = pv(i,2,k,2)
       endif
!
!  Upper horizontal axis (wind components)
!
       if (pv(i,ny,k,2) >= 0.0) then ! ouflow
         pu(i,ny+1,k,2) = pu(i,ny,k,2)
         pv(i,ny+1,k,2) = pv(i,ny,k,2)
       endif
     enddo
!
!  Left vertical axis (wind components)
!
     do j=1,ny+1
       if (pu(2,j,k,2) <= 0.0) then ! outflow
         pu(1,j,k,2) = pu(2,j,k,2)
         pv(1,j,k,2) = pv(2,j,k,2)
       endif
!
!  Right vertical axis (wind components)
!
       if (pu(nx,j,k,2) >= 0.0) then ! outflow
         pu(nx+1,j,k,2) = pu(nx,j,k,2)
         pv(nx+1,j,k,2) = pv(nx,j,k,2)
       endif
     enddo
   enddo
!
!  Diagnostics with fields at final time step (middle point at the next time step)
!  -------------------------------------------------------------------------------
!  
!  1) Temperature and liquid water from entropy
!
   call compute_t_qv_qc (h(:,:,:,2),pq(:,:,:,2),pi(:,:,2),qv,ql,t)
!
!  2) Geopotential height
!
   phi(:,:,nz+1) = rg*z_surf(:,:)
   do i=1,nx
     do j=1,ny
       p = sigma(nz)*pi(i,j,2) + ptop
       phatm = (p/p00)**rscp
       phatp = (pi(i,j,2)/p00)**rscp
       theta = 0.5*(t(i,j,nz)/phatm + ts(i,j)/phatp)
       qvm = 0.5*(qv(i,j,nz) + qs(i,j))
       phi(i,j,nz) = phi(i,j,nz+1) + cp*(theta*(1.+0.61*qvm)*(phatp-phatm))
       do k=nz-1,1,-1
         pm = sigma(k)*pi(i,j,2) + ptop
         pp = sigma(k+1)*pi(i,j,2) + ptop
         phatm = (pm/p00)**rscp
         phatp = (pp/p00)**rscp
         theta = 0.5*(t(i,j,k)/phatm + t(i,j,k+1)/phatp)
         qvm = 0.5*(qv(i,j,k) + qv(i,j,k+1))
         phi(i,j,k) = phi(i,j,k+1) + cp*(theta*(1.+0.61*qvm)*(phatp-phatm))
       enddo
     enddo
   enddo 
!
!  3) Surface fluxes (at mass points)
!
   do i=1,nx
     do j=1,ny
       pu_init(i,j,:) = 0.25*(pu(i+1,j+1,:,1) + pu(i,j+1,:,1) + pu(i+1,j,:,1) + pu(i,j,:,1))
       pv_init(i,j,:) = 0.25*(pv(i+1,j+1,:,1) + pv(i,j+1,:,1) + pv(i+1,j,:,1) + pv(i,j,:,1))
     enddo
   enddo
! 
   call surface_fluxes (t,qv,pu_init,pv_init,phi,pi(:,:,1),z0,ts,qs,ustar,tstar,qstar,ff,fg)
!
   fluxu(:,:) = pi(:,:,2)*ustar(:,:)*ustar(:,:)*pu_init(:,:,nz)/(sqrt(pu_init(:,:,nz)**2 + pv_init(:,:,nz)**2))
   fluxv(:,:) = pi(:,:,2)*ustar(:,:)*ustar(:,:)*pv_init(:,:,nz)/(sqrt(pu_init(:,:,nz)**2 + pv_init(:,:,nz)**2))
!
! Interpolate momentum fluxes at wind points 
!
   do i=1,nx+1
     do j=1,ny+1
       im = max(1,i-1)
       jm = max(1,j-1)
       ip = min(i,nx)
       jp = min(j,ny) 
       fluxum(i,j) = 0.25*(fluxu(ip,jp) + fluxu(ip,jm) + fluxu(im,jp) + fluxu(im,jm))
       fluxvm(i,j) = 0.25*(fluxv(ip,jp) + fluxv(ip,jm) + fluxv(im,jp) + fluxv(im,jm))
     enddo
   enddo
!
   fluxh(:,:) = pi(:,:,1)*(ustar(:,:)*tstar(:,:)*((pi(:,:,1)*sigma(nz) + ptop)/p00)**rscp/t(:,:,nz) + &
              &         lv*ustar(:,:)*qstar(:,:)/(cp*t(:,:,nz)))
   fluxq(:,:) = pi(:,:,1)*ustar(:,:)*qstar(:,:)
!
!  4) Exchange coefficients (at mass points)
!
   call eddy_coefficients (phi,ustar,ff,fg,kdifh,kdifm)
   print *,'exchange coefficients ok',lv*ustar(10,18)*qstar(10,18),cp*ustar(10,18)*tstar(10,18),pi(19,18,2)
   print *,'one time step ok'
!   
 enddo ! End of temporal loop
!
 return 
end subroutine time_evolution
