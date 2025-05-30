module consts
real, parameter ::  dt0 = 10.0,        &
                    psurf = 101325.0,  &
  &                 p00 = 100000.0,    &
  &                 ptop = 10000.0,    & 
  &                 tsurf = 299.0,     &
  &                 rg = 9.80665,      &
  &                 gradstd = 6.5E-3,  &
  &                 rd = 287.05,       & 
  &                 cp = 1005.0,       &
  &                 lv = 2.5E6,        &  
  &                 karman = 0.41,     &
  &                 rscp = rd/cp,      &    
  &                 rgsrd = rg/rd,     &
  &                 rgsgrd=rg/(gradstd*rd)
end module consts
