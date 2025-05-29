subroutine time_evolution_2Dversion
 use consts
 use vgrid , only  : sigma, nu, nz
 use hgrid , only  : nx, dx
! 
 implicit none
! 
 real, dimension (nx,nz)       :: pu_init, h_init, t
 real, dimension (nx)          :: pi_init
 real, dimension (nx,nz+1)     :: phi, nudotsigpi, nudot
 real, dimension (nx+1,nz,3)   :: pu
 real, dimension (nx,nz,3)     :: h
 real, dimension (nx,3)        :: pi
 real, dimension (nx)          :: z_surf, z0, ts, qs, dpidt
 real, dimension (nx)          :: phik, press
 real, dimension (nx+1)        :: pim
 real, dimension (nx+1)        :: zpu, pgradx
 real, dimension (nx+1)        :: advhu, diffhu
 real, dimension (nx)          :: advhh, diffhh
 real, dimension (nx+1,nz)     :: advvu
 real, dimension (nx,nz)       :: advvh
 real, dimension (nx,nz+1)     :: divu
 real, dimension (0:nz+1)      :: sigp
 real, dimension (nz)          :: kh, dnu
 integer                       :: i, k, nstep_max, im, ip, nstep, nt, nti, ntf, ntm, kk 
 real                          :: alpha, dt, wk, xpi, tmean, umean, psmean
 real                          :: kscale, filter, t1, t2, h0, b, x, zzz
 logical                       :: l_lbc_cst
!
! Define output files
!
 open (unit=88,file='../data/output1.dat')
!
 call cpu_time(time=t1) 
!
! Set-up various parameters
!
 xpi = acos(-1.0)
 alpha = 0.5
 wk = 0.0 ! pure Asselin filter
 kscale = (dx*dx)/(2.0*dt0)
 l_lbc_cst=.false.
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
!  Define depths of model layers in nu-coordinate system
!
 do k=1,nz   
   if (k == 1) then
     dnu(k) = 0.5*(nu(k) + nu(k+1)) - nu(k-1)
   elseif (k == nz) then
     dnu(k) = nu(k+1) - 0.5*(nu(k)+nu(k-1))
   else
     dnu(k) = 0.5*(nu(k+1) - nu(k-1))
   endif  
 enddo  
 
!
! Number of time steps
!
 nstep_max = 60
!
! vertical coordinate parameter : dsigma/dnu 
!
 sigp(:) = 4.0*(1. - nu(:)**3)/3.0
!
! Define initial conditions
!
 tmean  = 273.15
 umean  = 20.0
 psmean = 101315.0
!
!  Bell shaped mountain
!
 h0 = 5.0
 b  = 25.0E3 
! 
 do i=1,nx
   x = (float(i)-13.5)*dx
   z_surf(i) =  h0/(1.0 + (x/b)**2)
   pi_init(i) = psmean*exp(-rgsrd/tmean*z_surf(i))
 enddo
 
 ts(:) = tmean
 
 do i=1,nx
   do k=1,nz
     press(k) = sigma(k)*pi_init(i) + ptop
     h_init(i,k) = tmean*(p00/press(k))**rscp
     pu_init(i,k) = umean
   enddo  
 enddo 

 do i=1,nx
   do nt=1,3
     do k=1,nz
       h(i,k,nt)  = h_init(i,k)*pi_init(i)
       pu(i,k,nt) = pu_init(i,k)*pi_init(i)
     enddo
     pi(i,nt) = pi_init(i)
   enddo
 enddo
 
! Compute geopotential at initial time

 call compute_geopotential (z_surf,pi(:,1),ts,h(:,:,1),phi)
!
! Set-up boundary values for winds (add one line and one column)
!
 pu(nx+1,:,:) = pu(nx,:,:)
!
 print *,'Initial fields have been defined'
!
! Set-up vertical velocity to zero
!
 nudotsigpi(:,:) = 0.0
 dpidt(:) = 0.0
 nudot(:,:) = 0.0
!
! Start temporal loop
!
 do nstep=1,nstep_max
! 
! First time step = Euler forward - other time steps = Leapfrog
!
   if (nstep == 1) then
     nti = 1
     ntm = 1
     ntf = 2
     dt  = dt0
   else
     nti = 1
     ntm = 2
     ntf = 3
     dt  = 2.*dt0
   endif
!
! 1) Horizontal transfers (advection - pressure gradient - horizontal diffusion)   
!
   do k=1,nz ! Start vertical loop 
!
! a) Horizontal advection: leapfrog scheme
!
     do i=1,nx+1
       im = max(1,i-1)
       ip = min(i,nx)
       pim(i) = 0.5*(pi(ip,ntm) + pi(im,ntm))
     enddo 
     zpu(:) = pu(:,k,ntm)/pim(:)

     call advech_u (zpu,pu(:,k,ntm),advhu)
     call advech_h (zpu,h(:,k,ntm),advhh)    
!
! b) Horizontal diffusion: Euler forward scheme 
!    
     call diff_horiz_u (kh(k),pu(:,k,nti),diffhu)
     call diff_horiz_h (kh(k),h(:,k,nti),diffhh)
! 
     phik(:) = phi(:,k)
     press(:) = sigma(k)*pi(:,ntm) + ptop
!
     call pressure_gradx (k,pi(:,ntm),phik,t,press,pgradx)
!
     pu(:,k,ntf) = pu(:,k,nti) + dt*(advhu(:) + pgradx(:) + diffhu(:))
     h(:,k,ntf)  = h(:,k,nti)  + dt*(advhh(:) + diffhh(:))
     
     !i = 18
     !print *,'horizontal advection u',k,advhu(i)/pi(i,ntm)*86400.0
     !print *,'pressure gradient     ',k,pgradx(i)/pi(i,ntm)*86400.0
     !print *,'horizontal diffusion u',k,diffhu(i)/pi(i,ntm)*86400.0
     !print *,'horizontal advection h',k,advhh(i)/pi(i,ntm)*86400.0
     !print *,'horizontal diffusion h',k,diffhh(i)/pi(i,ntm)*86400.0
     !print *,'------'
!     
   enddo ! End vertical loop
!   
!  3- Compute surface pressure tendency
!   
   dpidt (:) = 0.0
   divu(:,:) = 0.0
   do i=1,nx
     do k=1,nz
       divu(i,k) = (pu(i+1,k,ntm) - pu(i,k,ntm))/dx 
       dpidt(i) =  dpidt(i) - dnu(k)*sigp(k)*divu(i,k)
     enddo
     pi(i,ntf) = pi(i,nti) + dt*dpidt(i)
   enddo 
!
!  4 - Compute vertical velocity
!
   nudotsigpi(:,:) = 0.0
   do i=1,nx
     do k=1,nz+1
       do kk=1,k
         nudotsigpi(i,k) = nudotsigpi(i,k) - dnu(kk)*sigp(kk)*(dpidt(i) + divu(i,kk))  
       enddo
     enddo
   enddo
!
   nudot(:,:) = 0.0
   do i=1,nx
     do k=2,nz
       nudot(i,k) = nudotsigpi(i,k)/(pi(i,ntm)*sigp(k))      
     enddo
   enddo
!   
!  2- Vertical transfers (vertical advection)
! 
   call advecv_u (dnu,nudot,sigp,pu(:,:,ntm),advvu)
   call advecv_h (dnu,nudot,sigp,h(:,:,ntm),advvh)
!
   pu(:,:,ntf) =  pu(:,:,ntf) + dt*advvu(:,:)
   h(:,:,ntf)  =  h(:,:,ntf)  + dt*advvh(:,:)   
   
!   do k=1,nz
!     i = 18
!     print *,'vertical advection',k,advvu(i,k)/pi(i,2)*86400.0,advvh(i,k)/pi(i,2)*86400.0
!   enddo
!   
!  5 - Apply an Asselin filter on prognostic variables (modified by Williams, 2009)
!
   if (nstep > 1) then
     do i=1,nx+1
       do k=1,nz
         filter = pu(i,k,3) + pu(i,k,1) - 2.0*pu(i,k,2)
         pu(i,k,2) = pu(i,k,2) + 0.5*alpha*wk*filter
         pu(i,k,3) = pu(i,k,3) - 0.5*alpha*(1.0 - wk)*filter    
       enddo
     enddo  
     do i=1,nx
       do k=1,nz
         filter = h(i,k,3) + h(i,k,1) - 2.0*h(i,k,2)
         h(i,k,2) = h(i,k,2) + 0.5*alpha*wk*filter
         h(i,k,3) = h(i,k,3) - 0.5*alpha*(1.0 - wk)*filter   
       enddo
       filter = pi(i,3) + pi(i,1) - 2.0*pi(i,2)
       pi(i,2) = pi(i,2) + 0.5*alpha*wk*filter
       pi(i,3) = pi(i,3) - 0.5*alpha*(1.0 - wk)*filter            
     enddo
   endif
!
!  Swapp variables for next model time step
!
   if (nstep > 1) then
     pu(:,:,1) = pu(:,:,2)
     pu(:,:,2) = pu(:,:,3)
     h (:,:,1) = h (:,:,2) 
     h (:,:,2) = h (:,:,3)
     pi(:,1)   = pi(:,2) 
     pi(:,2)   = pi(:,3)
   endif
!
!  Define upper boundary conditions (kept to the initial conditions)
!
   do nt=1,3
     pu(1:nx,1,nt) = umean*pi_init(:)
     h(:,1,nt)  = h_init(:,1)*pi_init(:)
   enddo
   pu(nx+1,1,:) = pu(nx,1,:)
!
!  Define lateral boundary conditions (thermodynamical variables)
!
   if (.not.l_lbc_cst) then 
!   
     h(nx,:,:)  = h(nx-1,:,:)
     h(1,:,:)   = h(2,:,:)
     pi(1,:)    = pi(2,:)
     pi(nx,:)   = pi(nx-1,:)

   else ! Keep boundary conditions fixed in time
!   
     do nt=1,3
       h(nx,:,nt)  = h_init(nx,:)
       h(1,:,nt)   = h_init(1,:)
       pi(1,nt)    = pi_init(1)
       pi(nx,nt)   = pi_init(nx)
     enddo
!
   endif 
        
   do k=1,nz
!
!  Left vertical axis (wind components)
!
     if (pu(2,k,2) <= 0.0) then ! outflow
       pu(1,k,2) = pu(2,k,2)
     endif
!
!  Right vertical axis (wind components)
!
     if (pu(nx,k,2) >= 0.0) then ! outflow
       pu(nx+1,k,2) = pu(nx,k,2)
     endif
     
   enddo 
! 
!  Diagnostics with fields at final time step (middle point at the next time step)
!  -------------------------------------------------------------------------------
!
!  Geopotential height
!
   call compute_geopotential (z_surf,pi(:,2),ts,h(:,:,2),phi)
!
!  Write fields in file
!
   if (nstep == nstep_max) then
     do i = 1,nx
       do k=1,nz  
         press(k) = sigma(k)*pi(i,ntm) + ptop
         zzz = (press(k)/p00)**rscp
         if (i == 12) print *,'temperature',k,h(i,k,ntm)/pi(i,ntm)*zzz - tmean         
         write (88,*) (i-1)*dx,phi(i,k)/rg,pu(i,k,ntm)/pi(i,ntm)-umean,h(i,k,ntm)/pi(i,ntm)*zzz - tmean
       enddo
     enddo
     close (unit=88)
   endif 
!   
 enddo ! End of temporal loop
!
 call cpu_time(time=t2)
 print *,'Total execution time for MESMOD 2D Version =',t2-t1,' sec for ',nstep_max*dt0/3600.0,' hours'
 return 
end subroutine time_evolution_2DVersion
