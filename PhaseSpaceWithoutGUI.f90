program PhaseSpace
  real :: phi0 = 0
  real :: phiN
  real :: phiN1
  
  real :: v0 = 5
  real :: vN
  real :: vN1
  
  real :: a0
  real :: aN
  real :: aN1
  
  real :: dt = 0.3
  
  integer :: j
  integer(kind=8) :: cycles
  
  real :: gamma = 0.00
  real :: omega = 1
  
  cycles = 200000
  
  ! Euler method:
  j = 0
  phiN = phi0
  vN = v0
  aN = -2 * gamma * v0 - omega**2 * sin(phi0/180*3.1415926535) 
  do while (j < cycles)
    aN = -2 * gamma * vN - omega**2 * sin(phiN/180*3.1415926535) 

    phiN = phiN + vN * dt * 0.01
    vN = vN + aN * dt * 0.01

    
    print *, j * dt, phiN, vN, (aN**2 / 2 - cos(phiN/180*3.1415926535)) 
    j = j + 1
  end do 
  
  ! Euler-Cromer method:
  j = 0
  phiN = phi0
  vN = v0
  aN = -2 * gamma * v0 - omega**2 * sin(phi0/180*3.1415926535) 
  do while (j < cycles)
    aN = -2 * gamma * vN - omega**2 * sin(phiN/180*3.1415926535) 
    vN = vN + aN * dt * 0.01
    phiN = phiN + vN * dt * 0.01
    
    print *, j * dt, phiN, vN, (aN**2 / 2 - cos(phiN/180*3.1415926535)) 
    j = j + 1
  end do

  ! Predictor - Corrector Method
  j = 0
  phiN = phi0
  vN = v0
  aN = -2 * gamma * v0 - omega**2 * sin(phi0/180*3.1415926535) 
  do while (j < cycles)
    aN = -2 * gamma * vN - omega**2 * sin(phiN/180*3.1415926535)
    phiN1 = phiN + vN * dt
    vN1 = vN + aN * dt
    aN1 = -2 * gamma * vN1 - omega**2 * sin(phiN1/180*3.1415926535)
    phiN = phiN + vN * dt + aN * dt**2 / 2 
    vN = vN + (aN + aN1) * dt / 2
    
    print *, j * dt, phiN, vN, (aN**2 / 2 - cos(phiN/180*3.1415926535)) 
    j = j + 1
  end do 
  
   
end program PhaseSpace
