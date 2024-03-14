module physical_parameters
  implicit none
  real(8), parameter :: m = 166 * 1.66E-27, kB = 1.38E-23, e0 = 8.854E-12
  real(8), parameter :: c = 299792458, g = 9.81, hbar = 1.05457E-34
  real(8), parameter :: pi = 3.141592, e = 2.71828
  real(8), parameter :: conversion = 0.1482E-24 * 1.113E-16
  real(8), parameter :: alpha_GS = 430 * conversion, lambd_trap = 488E-9
  real(8), parameter :: lambd583 = 583E-9, lambd631 = 631E-9, lambd841 = 841E-9
  real(8), parameter :: Gamma583 = 2*pi*180E3, Gamma631 = 2*pi*28E3
  real(8), parameter :: Gamma841 = 2*pi*8E3

end module physical_parameters
