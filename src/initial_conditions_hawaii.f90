subroutine initial_conditions
use consts
use vgrid, only  : nu, sigma
integer, parameter           :: nx=26, ny=26, nz=15
real, dimension (nx,ny,nz)   :: u,  v,  t, qv, ql 
real, dimension (nx,ny,nz)   :: pu, pv, h, pq
real, dimension (nx,ny,nz+1) :: phi, kdifh, kdifm
real, dimension (nx,ny)      :: z_surf, z0, pi, ts, qs
real, dimension (nx,ny)      :: ustar, tstar, qstar, ff, fg
integer, dimension (nx,ny)   :: soil_type, veg_type
real, dimension (0:40)       :: tref, zref, pref, uref, qref, href
real, dimension (nz)         :: zalt
integer                      :: im, jm, i, j, k
!
z_surf(:,:)     = 0.0
soil_type (:,:) = 1
veg_type(:,:)   = 1
z0(:,:)         = 0.0001
!
!  Define orography (hawaii island) from Nickerson and Magaziner (1976)
!
z_surf(12,8)  = 1.0    ; z_surf(11,9)  = 60.0   ; z_surf(12,9)  = 340.0  ; z_surf(10,10) = 240.0 
z_surf(11,10) = 700.0  ; z_surf(12,10) = 840.0  ; z_surf(13,10) = 430.0  ; z_surf(10,11) = 550.0 
z_surf(11,11) = 1650.0 ; z_surf(12,11) = 1580.0 ; z_surf(13,11) = 820.0  ; z_surf(14,11) = 120.0 
z_surf(10,12) = 670.0  ; z_surf(11,12) = 2130.0 ; z_surf(12,12) = 2290.0 ; z_surf(13,12) = 1860.0 
z_surf(14,12) = 700.0  ; z_surf(15,12) = 590.0  ; z_surf(16,12) = 210.0  ; z_surf(17,12) = 30.0 
z_surf(10,13) = 580.0  ; z_surf(11,13) = 2100.0 ; z_surf(12,13) = 3020.0 ; z_surf(13,13) = 2760.0 
z_surf(14,13) = 1650.0 ; z_surf(15,13) = 930.0  ; z_surf(16,13) = 1000.0 ; z_surf(17,13) = 930.0 
z_surf(18,13) = 640.0  ; z_surf(19,13) = 90.0   ; z_surf(10,14) = 850.0  ; z_surf(11,14) = 1900.0 
z_surf(12,14) = 3200.0 ; z_surf(13,14) = 3960.0 ; z_surf(14,14) = 2650.0 ; z_surf(15,14) = 1650.0 
z_surf(16,14) = 1280.0 ; z_surf(17,14) = 910.0  ; z_surf(18,14) = 620.0  ; z_surf(19,14) = 410.0 
z_surf(20,14) = 170.0  ; z_surf(9,15)  = 30.0   ; z_surf(10,15) = 1160.0 ; z_surf(11,15) = 1740.0 
z_surf(12,15) = 2600.0 ; z_surf(13,15) = 3230.0 ; z_surf(14,15) = 3020.0 ; z_surf(15,15) = 2270.0 
z_surf(16,15) = 1490.0 ; z_surf(17,15) = 880.0  ; z_surf(18,15) = 430.0  ; z_surf(19,15) = 170.0
z_surf(20,15) = 80.0   ; z_surf(9,16)  = 370.0  ; z_surf(10,16) = 1510.0 ; z_surf(11,16) = 1580.0 
z_surf(12,16) = 1940.0 ; z_surf(13,16) = 2150.0 ; z_surf(14,16) = 2150.0 ; z_surf(15,16) = 1830.0 
z_surf(16,16) = 1330.0 ; z_surf(17,16) = 720.0  ; z_surf(18,16) = 260.0  ; z_surf(19,16) = 30.0 
z_surf(8,17)  = 20.0   ; z_surf(9,17)  = 880.0  ; z_surf(10,17) = 1300.0 ; z_surf(11,17) = 1340.0 
z_surf(12,17) = 1600.0 ; z_surf(13,17) = 1810.0 ; z_surf(14,17) = 2130.0 ; z_surf(15,17) = 1910.0 
z_surf(16,17) = 1250.0 ; z_surf(17,17) = 670.0  ; z_surf(18,17) = 60.0   ; z_surf(19,17) = 20.0 
z_surf(9,18)  = 150.0  ; z_surf(10,18) = 500.0  ; z_surf(11,18) = 900.0  ; z_surf(12,18) = 1310.0 
z_surf(13,18) = 2470.0 ; z_surf(14,18) = 3990.0 ; z_surf(15,18) = 2680.0 ; z_surf(16,18) = 1580.0 
z_surf(17,18) = 820.0  ; z_surf(18,18) = 90.0   ; z_surf(10,19) = 90.0   ; z_surf(11,19) = 430.0 
z_surf(12,19) = 1010.0 ; z_surf(13,19) = 1430.0 ; z_surf(14,19) = 2190.0 ; z_surf(15,19) = 2130.0 
z_surf(16,19) = 1370.0 ; z_surf(17,19) = 520.0  ; z_surf(11,20) = 400.0  ; z_surf(12,20) = 820.0 
z_surf(13,20) = 1010.0 ; z_surf(14,20) = 1110.0 ; z_surf(15,20) = 850.0  ; z_surf(16,20) = 380.0 
z_surf(10,21) = 110.0  ; z_surf(11,21) = 1070.0 ; z_surf(12,21) = 1040.0 ; z_surf(13,21) = 580.0 
z_surf(14,21) = 270.0  ; z_surf(10,22) = 340.0  ; z_surf(11,22) = 590.0  ; z_surf(12,22) = 1.0 
z_surf(10,23) = 3.0 
!
!
z_surf(:,:) = 0.0
!
! Define soil type (from Nickerson and Magaziner, 1976)
!
soil_type(12,8)  = 6  ; soil_type(11,9)  = 5  ; soil_type(12,9)  = 5  ; soil_type(10,10) = 5     
soil_type(11,10) = 3  ; soil_type(12,10) = 3  ; soil_type(13,10) = 3  ; soil_type(10,11) = 4     
soil_type(11,11) = 4  ; soil_type(12,11) = 4  ; soil_type(13,11) = 4  ; soil_type(14,11) = 4     
soil_type(10,12) = 4  ; soil_type(11,12) = 5  ; soil_type(12,12) = 5  ; soil_type(13,12) = 4      
soil_type(14,12) = 4  ; soil_type(15,12) = 5  ; soil_type(16,12) = 2  ; soil_type(17,12) = 2    
soil_type(10,13) = 4  ; soil_type(11,13) = 5  ; soil_type(12,13) = 5  ; soil_type(13,13) = 5      
soil_type(14,13) = 4  ; soil_type(15,13) = 4  ; soil_type(16,13) = 3  ; soil_type(17,13) = 4     
soil_type(18,13) = 4  ; soil_type(19,13) = 3  ; soil_type(10,14) = 4  ; soil_type(11,14) = 4      
soil_type(12,14) = 5  ; soil_type(13,14) = 5  ; soil_type(14,14) = 5  ; soil_type(15,14) = 4      
soil_type(16,14) = 4  ; soil_type(17,14) = 4  ; soil_type(18,14) = 4  ; soil_type(19,14) = 4     
soil_type(20,14) = 4  ; soil_type(9,15)  = 2  ; soil_type(10,15) = 4  ; soil_type(11,15) = 5      
soil_type(12,15) = 5  ; soil_type(13,15) = 5  ; soil_type(14,15) = 5  ; soil_type(15,15) = 5      
soil_type(16,15) = 4  ; soil_type(17,15) = 4  ; soil_type(18,15) = 4  ; soil_type(19,15) = 4    
soil_type(20,15) = 4  ; soil_type(9,16)  = 4  ; soil_type(10,16) = 4  ; soil_type(11,16) = 3      
soil_type(12,16) = 5  ; soil_type(13,16) = 5  ; soil_type(14,16) = 5  ; soil_type(15,16) = 5      
soil_type(16,16) = 4  ; soil_type(17,16) = 4  ; soil_type(18,16) = 4  ; soil_type(19,16) = 4    
soil_type(8,17)  = 3  ; soil_type(9,17)  = 4  ; soil_type(10,17) = 3  ; soil_type(11,17) = 5      
soil_type(12,17) = 5  ; soil_type(13,17) = 5  ; soil_type(14,17) = 2  ; soil_type(15,17) = 4      
soil_type(16,17) = 4  ; soil_type(17,17) = 4  ; soil_type(18,17) = 4  ; soil_type(19,17) = 4    
soil_type(9,18)  = 5  ; soil_type(10,18) = 5  ; soil_type(11,18) = 5  ; soil_type(12,18) = 2      
soil_type(13,18) = 2  ; soil_type(14,18) = 5  ; soil_type(15,18) = 3  ; soil_type(16,18) = 4      
soil_type(17,18) = 4  ; soil_type(18,18) = 4  ; soil_type(10,19) = 5  ; soil_type(11,19) = 2     
soil_type(12,19) = 2  ; soil_type(13,19) = 2  ; soil_type(14,19) = 2  ; soil_type(15,19) = 4      
soil_type(16,19) = 4  ; soil_type(17,19) = 4  ; soil_type(11,20) = 2  ; soil_type(12,20) = 2     
soil_type(13,20) = 3  ; soil_type(14,20) = 4  ; soil_type(15,20) = 4  ; soil_type(16,20) = 4     
soil_type(10,21) = 2  ; soil_type(11,21) = 3  ; soil_type(12,21) = 4  ; soil_type(13,21) = 4     
soil_type(14,21) = 4  ; soil_type(10,22) = 3  ; soil_type(11,22) = 4  ; soil_type(12,22) = 4   
soil_type(10,23) = 4
!
!  Define vegetation type (from Nickerson and Magaziner, 1976)
!
!  1: none - 2: short grass - 3: tall grass - 4: shrub - 5: forest
!
veg_type(12,8)  = 1  ; veg_type(11,9)  = 1  ; veg_type(12,9)  = 1  ; veg_type(10,10) = 1     
veg_type(11,10) = 4  ; veg_type(12,10) = 5  ; veg_type(13,10) = 3  ; veg_type(10,11) = 5     
veg_type(11,11) = 4  ; veg_type(12,11) = 5  ; veg_type(13,11) = 5  ; veg_type(14,11) = 3     
veg_type(10,12) = 5  ; veg_type(11,12) = 1  ; veg_type(12,12) = 1  ; veg_type(13,12) = 5      
veg_type(14,12) = 3  ; veg_type(15,12) = 1  ; veg_type(16,12) = 2  ; veg_type(17,12) = 2    
veg_type(10,13) = 5  ; veg_type(11,13) = 1  ; veg_type(12,13) = 1  ; veg_type(13,13) = 1      
veg_type(14,13) = 5  ; veg_type(15,13) = 4  ; veg_type(16,13) = 4  ; veg_type(17,13) = 4     
veg_type(18,13) = 4  ; veg_type(19,13) = 4  ; veg_type(10,14) = 5  ; veg_type(11,14) = 4      
veg_type(12,14) = 1  ; veg_type(13,14) = 1  ; veg_type(14,14) = 1  ; veg_type(15,14) = 4      
veg_type(16,14) = 5  ; veg_type(17,14) = 5  ; veg_type(18,14) = 5  ; veg_type(19,14) = 5     
veg_type(20,14) = 5  ; veg_type(9,15)  = 2  ; veg_type(10,15) = 4  ; veg_type(11,15) = 1      
veg_type(12,15) = 1  ; veg_type(13,15) = 1  ; veg_type(14,15) = 1  ; veg_type(15,15) = 1      
veg_type(16,15) = 5  ; veg_type(17,15) = 5  ; veg_type(18,15) = 4  ; veg_type(19,15) = 4    
veg_type(20,15) = 4  ; veg_type(9,16)  = 4  ; veg_type(10,16) = 4  ; veg_type(11,16) = 3      
veg_type(12,16) = 1  ; veg_type(13,16) = 1  ; veg_type(14,16) = 1  ; veg_type(15,16) = 1      
veg_type(16,16) = 5  ; veg_type(17,16) = 5  ; veg_type(18,16) = 5  ; veg_type(19,16) = 4    
veg_type(8,17)  = 3  ; veg_type(9,17)  = 4  ; veg_type(10,17) = 3  ; veg_type(11,17) = 1      
veg_type(12,17) = 1  ; veg_type(13,17) = 1  ; veg_type(14,17) = 4  ; veg_type(15,17) = 4      
veg_type(16,17) = 5  ; veg_type(17,17) = 5  ; veg_type(18,17) = 5  ; veg_type(19,17) = 3    
veg_type(9,18)  = 1  ; veg_type(10,18) = 1  ; veg_type(11,18) = 1  ; veg_type(12,18) = 3      
veg_type(13,18) = 3  ; veg_type(14,18) = 1  ; veg_type(15,18) = 4  ; veg_type(16,18) = 5      
veg_type(17,18) = 5  ; veg_type(18,18) = 4  ; veg_type(10,19) = 1  ; veg_type(11,19) = 4     
veg_type(12,19) = 4  ; veg_type(13,19) = 4  ; veg_type(14,19) = 2  ; veg_type(15,19) = 4      
veg_type(16,19) = 5  ; veg_type(17,19) = 5  ; veg_type(11,20) = 2  ; veg_type(12,20) = 2     
veg_type(13,20) = 4  ; veg_type(14,20) = 4  ; veg_type(15,20) = 5  ; veg_type(16,20) = 5     
veg_type(10,21) = 2  ; veg_type(11,21) = 4  ; veg_type(12,21) = 5  ; veg_type(13,21) = 5     
veg_type(14,21) = 5  ; veg_type(10,22) = 3  ; veg_type(11,22) = 5  ; veg_type(12,22) = 4   
veg_type(10,23) = 3  
!
!  Define pi quantity = (ps - pt) and surface roughness length (m)
! 
do i=1,nx
  do j=1,ny
    pi(i,j)=psurf*(1.0 - gradstd*z_surf(i,j)/tsurf)**rgsgrd - ptop
!
    if ((soil_type(i,j) == 5 .or. soil_type(i,j) == 2) .and. (veg_type(i,j) < 4)) then
      z0(i,j) = 0.01
    endif
    if ((soil_type(i,j) == 2) .and. (veg_type(i,j) == 4)) then
      z0(i,j) = 0.2
    endif   
    if ((soil_type(i,j) == 3) .and. ((veg_type(i,j) == 2) .or. veg_type(i,j) == 3)) then
      z0(i,j) = 0.05
    endif   
    if ((soil_type(i,j) == 3) .and. (veg_type(i,j) == 4)) then
      z0(i,j) = 0.5
    endif  
    if (soil_type(i,j) == 4) then
      if (veg_type(i,j) == 2 .or. veg_type(i,j) == 3) z0(i,j) = 0.5
      if (veg_type(i,j) == 4) z0(i,j) = 1.0
      if (veg_type(i,j) == 5) z0(i,j) = 3.0
    endif
  enddo
enddo

z0(:,:) = 0.1

!
! Define a reference temperature profile at sea level
!
tref(1) = tsurf - gradstd*500.0
zref(1) = 500.0
!
tref(0) = tsurf
zref(0) = 0.0
!
pref(0) = psurf
!
do kk=2,26
  tref(kk) = tref(kk-1) - gradstd*500.0
  if (kk == 3) tref(kk) = tref(kk) + 4.0 ! inversion at top of pbl
  zref(kk) = 500.*real(kk)
enddo
zref(27) = 17000.0 ; tref(27) = 213.0
zref(28) = 21000.0 ; tref(28) = 218.0
zref(29) = 24000.0 ; tref(29) = 223.0
zref(30) = 31000.0 ; tref(30) = 237.0
zref(31) = 34000.0 ; tref(31) = 244.0
zref(32) = 36000.0 ; tref(32) = 251.0
zref(33) = 40000.0 ; tref(33) = 263.0
zref(34) = 49000.0 ; tref(34) = 276.0
zref(35) = 52000.0 ; tref(35) = 273.0
zref(36) = 54000.0 ; tref(36) = 269.0
zref(37) = 58000.0 ; tref(37) = 256.0
zref(38) = 66000.0 ; tref(38) = 225.0
zref(39) = 68000.0 ; tref(39) = 216.0
zref(40) = 70000.0 ; tref(40) = 208.0
!
! Define the reference temperature profile at sea level
!
pref(1) = psurf*exp(-rgsrd*zref(1)/(0.5*(tref(1)+tsurf)))
do kk=2,40
  pref(kk) = pref(kk-1)*exp(-rgsrd*(zref(kk) - zref(kk-1))/(0.5*(tref(kk)+tref(kk-1))))
enddo
!
! Define a reference profile for wind and relative humidity 
!
uref(0) = 2.5 ! set surface value to non zero since 1st value from reference profile is at 500 m

do kk=1,40
  if (zref(kk) < 5000.) then
  !  uref(kk) = -5.0
  endif
  if (zref(kk) > 7000.) then
  !  uref(kk) =  5.0
  endif
  if (zref(kk) >= 5000. .and. zref(kk) <= 7000.) then
  !  uref(kk) = -5.0 + 10./2000.*(zref(kk) - 5000.)
  endif

  uref(kk) = 5.0
!
  if (zref(kk) < 3000.) then
    href(kk) = 0.80
  endif
  if (zref(kk) > 6000.) then
    href(kk) = 0.20
  endif
  if (zref(kk) >= 3000. .and. zref(kk) <= 6000.) then
    href(kk) = 0.80 - 0.6/3000.*(zref(kk) - 3000.)
  endif
  if (zref(kk) > 12000.) href(kk) = 0.05 ! dry stratosphere
enddo
!
! Define 3d-fields (temperature, wind, specific humidity)
!
do i=1,nx
  do j=1,ny
!
! Define surface href (surface relative humidity) according to soil type
!
    if (soil_type(i,j) == 1) href(0) = 1.0  ! sea
    if (soil_type(i,j) == 2) href(0) = 0.2  ! dry
    if (soil_type(i,j) == 3) href(0) = 0.6  ! semi-moist
    if (soil_type(i,j) == 4) href(0) = 0.8  ! wet
    if (soil_type(i,j) == 5) href(0) = 0.1  ! lava
    if (soil_type(i,j) == 6) href(0) = 0.5  ! sand

href(0) = 0.2

!
! Averaged surface pressure at wind points
!
    im = max(1,i-1)
    jm = max(1,j-1)
    pim = 0.25*(pi(i,j) + pi(im,j) + pi(i,jm) + pi(im,jm))
!
    do k=nz,1,-1
      p = sigma(k)*pi(i,j) + ptop
      do kk=0,39        
        if (pref(kk) >= p .and. pref(kk+1) < p) then
          weight1 = (p - pref(kk))/(pref(kk) - pref(kk+1)) ! linear interpolation in pressure (p)
          weight2 = (log(p) - log(pref(kk)))/(log(pref(kk))-log(pref(kk+1))) ! linear interpolation in log(p)
          t(i,j,k)  = tref(kk) + (tref(kk) - tref(kk+1))*weight2
          rh        = href(kk) + (href(kk) - href(kk+1))*weight2 
          qv(i,j,k) = 0.0 ! rh*qsat(p,t(i,j,k))
          u(i,j,k)  = uref(kk) + (uref(kk) - uref(kk+1))*weight2 
          zalt(k)   = zref(kk) + (zref(kk) - zref(kk+1))*weight2 ! height approximate
        endif        
      enddo  
!
!  Initial variables to be stored 
!
    h(i,j,k)  = pi(i,j)*(log(t(i,j,k)*(p00/p)**rscp) + lv*qv(i,j,k)/(cp*t(i,j,k)))
    pq(i,j,k) = pi(i,j)*qv(i,j,k)
    pu(i,j,k) = pim*u(i,j,k)
    pv(i,j,k) = 0.0  
! 
    enddo  
!
!  Smooth wind profile near the ground (log profile)
!
    if (zalt(nz) < z_surf(i,j)) then
      print *, 'wrong estimate of model height',zalt(nz),z_surf(i,j)
      stop
    endif
    u(i,j,nz)  = u(i,j,nz-1)/log((zalt(nz-1)-z_surf(i,j))/(zalt(nz)-z_surf(i,j)))
    pu(i,j,nz) = pim*u(i,j,nz)   
!
!  Define surface temperature and surface specific humidity
!
    ps = pi(i,j) + ptop
    do kk=0,39
      if (pref(kk) >= ps .and. pref(kk+1) < ps) then
        ts(i,j) = tref(kk) + (tref(kk) - tref(kk+1))/(pref(kk) - pref(kk+1))*(ps - pref(kk))
        qs(i,j) = href(0)*qsat(ps,ts(i,j))
      endif
    enddo 
  enddo
enddo
!
!  Compute geopotential height
!
phi(:,:,nz+1) = rg*z_surf(:,:)
do i=1,nx
  do j=1,ny
    p = sigma(nz)*pi(i,j) + ptop
    phatp = (p/p00)**rscp
    phatm = (pi(i,j)/p00)**rscp
    theta = 0.5*(t(i,j,nz)/phatp + ts(i,j)/phatm)
    qvm = 0.5*(qv(i,j,nz) + qs(i,j))
    phi(i,j,nz) = phi(i,j,nz+1) + cp*(theta*(1.+0.61*qvm)*(phatm-phatp))
    do k=nz-1,1,-1
      pp = sigma(k)*pi(i,j) + ptop
      pm = sigma(k+1)*pi(i,j) + ptop
      phatp = (pp/p00)**rscp
      phatm = (pm/p00)**rscp
      theta = 0.5*(t(i,j,k)/phatp + t(i,j,k+1)/phatm)
      qvm = 0.5*(qv(i,j,k) + qv(i,j,k+1))
      phi(i,j,k) = phi(i,j,k+1) + cp*(theta*(1.+0.61*qvm)*(phatm-phatp))
    enddo
  enddo
enddo 
!
! Write fields for initial conditions
!
 open (unit=300,file='initial_conditions_hawaii4.dat')
!
 write (300,*) h
 write (300,*) pu
 write (300,*) pv
 write (300,*) pq
 write (300,*) t
 write (300,*) phi
 write (300,*) pi
 write (300,*) z0
 write (300,*) ts
 write (300,*) qs
 write (300,*) z_surf
!
 close (unit=300)
!
! Diagnostic of clouds at the end of the time step - t, h and qv are recomputed accordingly
!
call compute_t_qv_qc (h,pq,pi,qv,ql,t)
!
! Compute surface fluxes
!
call surface_fluxes (t,pq,pu,pv,phi,pi,z0,ts,qs,ustar,tstar,qstar,ff,fg)
!
! Compute exchange coefficients
!
call eddy_coefficients (phi,ustar,ff,fg,kdifh,kdifm)
!
! Write in file for plotting
!
do i=1,nx
  do j=1,ny
    write (110,*) i,j,z_surf(i,j)
    write (111,*) i,j,pi(i,j)
    write (112,*) i,j,ts(i,j)
    write (113,*) i,j,z0(i,j)
  enddo
  write (110,*)
  write (111,*)
  write (112,*)
  write (113,*)
enddo
!
end subroutine initial_conditions
