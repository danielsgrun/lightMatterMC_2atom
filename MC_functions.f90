! Suite of functions to be used by the LightMatterMC "simulationLoss.f90"
! subroutine for the simulation of a single atom inside an optical tweezer
! in the presence of a light field

module MC_functions
  use physical_parameters
  implicit none
  
contains
  
  function readSpontEmiFile(spontCase)
    implicit none
    integer :: i
    character(len=1) :: spontCase
    integer, parameter :: arrayLength = 10000
    real(8), dimension(arrayLength, 3) :: readSpontEmiFile

    if (spontCase == 's') then
       open(unit=10, file="spontEmission_Sigma.txt")
    else if (spontCase == 'i') then
       open(unit=10, file="spontEmission_isotropic.txt")
    endif
    do i=1,arrayLength
       read(10,*) readSpontEmiFile(i,:)
    enddo
    close(unit=10)
  end function readSpontEmiFile


  
  real(8) function Rscatt(Gamma, s, Delta)
    implicit none
    real(8) :: Gamma, s, Delta
    Rscatt = Gamma/2. * (s/(1+s+4*(Delta/Gamma)**2))
  end function Rscatt


  
  function getR12(rv1,rv2)
    implicit none
    real(8), dimension(6) :: rv1, rv2
    real(8) :: x1,x2, y1,y2, z1,z2
    real(8) :: getR12
    real(8), parameter :: epsilon=1E-10 ! parameter to avoid singularity

    x1=rv1(1); y1=rv1(2); z1=rv1(3);
    x2=rv2(1); y2=rv2(2); z2=rv2(3)

    getR12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2+epsilon**2)
  end function getR12
    
  
  real(8) function trapPotential(rv1, rv2, P, w0, alpha, lambd, C3)
    implicit none 
    real(8) :: P,w0,alpha,lambd,zR,U0 
    real(8), dimension(6) :: rv1, rv2
    real(8) :: C3, rr,LJ, x1,y1,z1, x2,y2,z2

    x1=rv1(1); y1=rv1(2); z1=rv1(3) 
    x2=rv2(1); y2=rv2(2); z2=rv2(3) 
    
    rr = getR12(rv1,rv2)

    LJ = C3/(rr**3)

!!$    LJ = 0.
    
    U0 = P*alpha / (pi*e0*c*w0**2)
    zR = pi * w0**2 / lambd

    trapPotential = -U0/(1+(z1/zR)**2)*exp(-2/(w0**2)*(x1**2+y1**2)/(1+(z1/zR)**2)) + LJ
    
  end function trapPotential


  
  function trapPotDerivs(rv1, rv2, U, w0, lambd, C3) result(derivs)
    implicit none
    real(8) :: U,w0,lambd,ax,ay,az,z_term,zR
    real(8) :: C3, rr, x1,y1,z1, x2,y2,z2
    real(8), dimension(6) :: rv1, rv2
    real(8), dimension(3) :: derivs
    real(8) :: LJ_diff_x, LJ_diff_y, LJ_diff_z
    
    x1=rv1(1); y1=rv1(2); z1=rv1(3)
    x2=rv2(1); y2=rv2(2); z2=rv2(3)
    
    zR = pi*w0**2/lambd
    z_term = sqrt(1+(z1/zR)**2)
    
    rr = getR12(rv1,rv2)
    
    LJ_diff_x = -(3*C3*(x1-x2))/rr**5
    LJ_diff_y = -(3*C3*(y1-y2))/rr**5
    LJ_diff_z = -(3*C3*(z1-z2))/rr**5

!!$    LJ_diff_x = 0.
!!$    LJ_diff_y = 0.
!!$    LJ_diff_z = 0.
    
    ax = 4/m *x1/(w0**2*z_term**2)*U -1/m* LJ_diff_x
    ay = 4/m *y1/(w0**2*z_term**2)*U -1/m* LJ_diff_y
    az = 2/m *z1/(zR**4*w0**2*z_term**4)*(zR**2*(w0**2-2*(x1**2+y1**2))+w0**2*z1**2)*U - 1/m * LJ_diff_z
    
    derivs(1)=ax; derivs(2)=ay; derivs(3)=az
    
  end function trapPotDerivs

  

  real(8) function kineticEnergy(v)
    implicit none
    real(8), dimension(3) :: v
    kineticEnergy = m/2 * (v(1)**2+v(2)**2+v(3)**2)
  end function kineticEnergy


  
  function checkLost(rvVector, P, w0, C3)
    implicit none
    real(8), dimension(6) :: rvVector
    real(8) :: P, w0, potEnergy, kinEnergy, totalEnergy, C3
    !real(8), external :: trapPotential, kineticEnergy
    real(8), dimension(3) :: r,v
    real(8) :: checkLost

    r = rvVector(1:3); v = rvVector(4:6)
	
    potEnergy = trapPotential(rvVector,rvVector,P,w0,alpha_GS,lambd_trap,C3)
    kinEnergy = kineticEnergy(v)
    totalEnergy = potEnergy + kinEnergy

    checkLost = 0.0
    
    if (totalEnergy >= 0) then
       checkLost = 1.0
    endif
	
  end function checkLost

  
    
  function trapEvol(rv1, rv2, P, w0, alpha, C3)
    implicit none
    real(8), dimension(6) :: rv1, rv2, trapEvol 
    real(8), dimension(3) :: r1, r2, v1, accels1 
    real(8) :: P, w0, ax1, ay1, az1, U1, alpha, C3
    integer :: i

    r1 = rv1(1:3); r2 = rv2(1:3); v1 = rv1(4:6)
    
    U1 = trapPotential(rv1,rv2,P,w0,alpha,lambd_trap,C3)
    accels1 = trapPotDerivs(rv1,rv2,U1,w0,lambd_trap,C3)

    ax1=accels1(1); ay1=accels1(2); az1=accels1(3)
    
    do i=1,3 
       trapEvol(i) = v1(i)
    enddo

    trapEvol(4)=ax1; trapEvol(5)=ay1; trapEvol(6)=az1-g
    
  end function trapEvol
  

  
  function odeRK4_solver(y01, y02, dt, P, w0, alpha, C3)
    implicit none
    real(8) :: dt, P, w0, alpha, C3
    real(8), dimension(6) :: y01, y02, dtArray, k1, k2, k3, k4
    real(8), dimension(6) :: odeRK4_solver

    dtArray = [dt,dt,dt,dt,dt,dt]
    
    k1 = dtArray * trapEvol(y01, y02, P, w0, alpha,C3)
    k2 = dtArray * trapEvol(y01+0.5*k1, y02, P, w0, alpha,C3)
    k3 = dtArray * trapEvol(y01+0.5*k2, y02, P, w0, alpha,C3)
    k4 = dtArray * trapEvol(y01+k3, y02, P, w0, alpha,C3)

    !print*, k1
    
    odeRK4_solver = y01 + (k1 + 2*k2 + 2*k3 + k4)/6.0

  end function odeRK4_solver

  
  
  function timeEvolution(rv1, rv2, dt, P, w0, alphas, C3, GSFlags)

    implicit none
    real(8) :: dt, P, w0, C3
    real(8), dimension(2) :: alphas, GSFlags
    real(8) :: C3zero
    real(8), dimension(6) :: rv1, rv2, timeEvolution, GSTE
    real(8), dimension(6) :: ex1bTE, ex2bTE, exTE, nullVec
    
    C3zero = 0.0
    nullVec = [0., 0., 0., 0., 0., 0.]

    if (GSFlags(1) > 0 .and. GSFlags(2) > 0) then   
       GSTE = odeRK4_solver(rv1, rv2, dt, P, w0, alphas(1), C3zero)
       ex1bTE = nullVec
       ex2bTE = nullVec
    elseif (GSFlags(1) == 0 .and. GSFlags(2) == 0) then
       GSTE = nullVec
       ex1bTE = odeRK4_solver(rv1, rv2, dt, P, w0, alphas(2), C3zero)
       ex2bTE = nullVec
    else
       GSTE = nullVec
       ex1bTE = nullVec
       ex2bTE = odeRK4_solver(rv1, rv2, dt, P, w0, alphas(2), C3)
    endif

!!$    GSTE = odeRK4_solver(rv1, rv2, dt, P, w0, alphas(1), C3zero)
!!$    ex1bTE = [0.,0.,0.,0.,0.,0.]
!!$    ex2bTE = [0.,0.,0.,0.,0.,0.]
       
    timeEvolution = GSTE + ex1bTE + ex2bTE
    ! timeEvolution =  odeRK4_solver(rv1, rv2, dt, P, w0, alphas(1), C3zero)

  end function timeEvolution

    

  real(8) function getAcStarkShift(rvVector, P, w0, alpha_E)
    implicit none
    real(8), dimension(6) :: rvVector
    real(8), parameter :: C3=0      ! ACS doesn't depend on DDI
    real(8) :: P, w0, alpha_E, U_g, U_e
    	
    U_g = trapPotential(rvVector,rvVector,P,w0,alpha_GS,lambd_trap, C3)
    U_e = trapPotential(rvVector,rvVector,P,w0,alpha_E,lambd_trap, C3)
  
    getAcStarkShift = -(U_e-U_g)/hbar
  end function getAcStarkShift


  function getOneBodyDelta(rv, P, w0, delta, alpha, lambd, absProj)
    implicit none
    real(8), dimension(6) :: rv
    real(8), dimension(3) :: absProj
    real(8) :: P, w0, alpha, lambd, delta
    real(8) :: acStarkShift, DopplerShift, getOneBodyDelta

    acStarkShift = getACStarkShift(rv, P, w0, alpha)
    DopplerShift = getDopplerShift(rv, absProj, lambd)

    getOneBodyDelta = 2*pi*delta + DopplerShift + acStarkShift

  end function getOneBodyDelta


  real(8) function getDopplerShift(rvVector, absProj, lambd)
    implicit none
    real(8), dimension(3) :: absProj, v
    real(8), dimension(6) :: rvVector
    real(8) :: absProjLength, k, lambd, DopplerShift
    integer :: i
    
    v = rvVector(4:6)
    absProjLength = 0.0
    do i=1,3
       absProjLength = absProjLength + absProj(i)**2
    enddo

    absProjLength = sqrt(absProjLength)
    do i=1,3
       absProj(i) = absProj(i)/absProjLength
    enddo

    k = 2*pi/lambd
    DopplerShift = dot_product(absProj, v)

    getDopplerShift = -k*DopplerShift
  end function getDopplerShift

  
  function recoilVel(absProj, lambd, case, externalSpontDir)
    implicit none
    real(8), dimension(3) :: absProj, recoilDir, spontDir, externalSpontDir
    real(8) :: lambd, k, absProjLength
    real(8), dimension(3) :: random3D
    character(len=3) :: case
    integer :: i
    real(8), dimension(3) :: recoilVel
    
    absProjLength = 0.0

    call random_number(random3D)
    
    !spontDir = 1.0-2*random3D
    !spontDirLength = 0.0
    
    do i=1,3
       absProjLength = absProjLength + absProj(i)**2
       !spontDirLength = spontDirLength + spontDir(i)**2
    enddo
    
    !spontDirLength = sqrt(spontDirLength)
    absProjLength = sqrt(absProjLength)
    
    do i=1,3
       absProj(i) = absProj(i)/absProjLength
       !spontDir(i) = spontDir(i)/spontDirLength
    enddo

    if (case=="abs") then
       recoilDir = absProj
       !print*, "absorption"
    else if (case=="emi") then
       spontDir = externalSpontDir
       recoilDir = spontDir
       !print*, "emission"
    endif

    k = 2*pi/lambd

    recoilVel = hbar*k/m * recoilDir

    !print*, recoilDir

  end function recoilVel


  function getDetuningFromC3(C3, rv1, rv2)
    implicit none
    real(8), dimension(6) :: rv1, rv2
    real(8) :: C3, LJ, rr, getDetuningFromC3
    real(8), parameter :: epsilon = 1e-10 ! parameter to avoid singularities
    
    rr = getR12(rv1,rv2)

    LJ = C3/(rr**3)

    getDetuningFromC3 = -LJ/hbar

  end function getDetuningFromC3
       
end module MC_functions

    
  
  
