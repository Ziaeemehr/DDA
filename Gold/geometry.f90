program geometry
implicit none
real(8) :: d, d0, R, R2, u
real(8), allocatable :: x(:), y(:), z(:),rat(:,:)
integer :: t, N, nat

print *, "Enter diameter in nanometer :"
read *, d
R  =  d/2.d0
d0 = (4.d0 * pi/(3.d0 * N ))**(1./3.) * R;
n = floor(2.d0*R/d0)  

