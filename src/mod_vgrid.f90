module vgrid
implicit none
integer, parameter :: nz=15
real, dimension (0:nz+1) :: nu, sigma
data nu    /0.0000, 0.0333, 0.1000, 0.1667, 0.2333, 0.3000, &
   &        0.3667, 0.4333, 0.5000, 0.5667, 0.6333, 0.7000, &
   &        0.7667, 0.8333, 0.9000, 0.9667, 1.0000/
data sigma /0.0000, 0.0444, 0.1333, 0.2220, 0.3101, 0.3973, &
   &        0.4829, 0.5660, 0.6458, 0.7212, 0.7908, 0.8533, &
   &        0.9071, 0.9504, 0.9813, 0.9978, 1.0000/
end module vgrid
