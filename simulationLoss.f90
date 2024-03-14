! Monte-Carlo simulation of survival of a single atom inside the Tweezer during
! imaging with near-resonant light.
! Includes photon absorption and spontaneous emission in the presence of
! a light field

subroutine create_simulation_loss(C3Vals, firstInitialCond, P,T,w0,titf,s0,nBeams,lambd,&
     Gammas,absProj,delta,alpha,spontCase,solution)

  use physical_parameters
  use MC_functions
  
  implicit none

  integer :: i, j, k, index, atom1, atom2, availablePos
  integer :: kickIndex, indexReversed, case, lostFlagSum
  
  integer, parameter :: arrayLength = 1E4
  real(8), parameter :: arraySize = 1E4, C3zero = 0.0
  real(8), dimension(arrayLength, 3) :: spontEmissionArray
  real(8), dimension(3) :: recoil, externalSpontDir
  integer, intent(in) :: nBeams
  integer, parameter :: nAtoms = 2
  real(8), intent(in) ::  P,T,w0
  real(8), dimension(3), intent(in) :: titf
  real(8), dimension(nBeams, 3), intent(in) :: absProj
  real(8), dimension(nAtoms,6), intent(in) :: firstInitialCond
  real(8), dimension(nBeams), intent(in) :: s0, lambd, Gammas, delta, alpha
  real(8), dimension(nAtoms+1), intent(out) :: solution
  real(8), dimension(nAtoms,6) :: sols, initialCond
  real(8) :: ti, tf, dt, decayCond, auxNum, phScatt
    
  character(len=1), intent(in) :: spontCase

  real(8), dimension(nBeams, nAtoms) :: DopplerShift, acStarkShift, deltaOneBody
  real(8), dimension(nBeams, 2) :: deltaTotal
  real(8), dimension(nAtoms) :: passedTime, waitTime, lost, C3List
  real(8), dimension(nAtoms) :: oneBodyExcit, kickFlag, excitFlag
  real(8), dimension(nAtoms) :: twoBodyExcit, GSFlag, randomNumbers, auxRatios
  
  real(8), dimension(3), intent(in) :: C3Vals
  real(8), dimension(3) :: DeltaC3
  
  integer, dimension(nAtoms) :: partnerAtom
  
  real(8), dimension(2) :: auxRatio12, scattProb, alphas, phScattEach, auxValC3
  integer, dimension(2) :: atoms
  real(8) :: currentTime, auxRatio
  
  call random_seed
  
  spontEmissionArray = readSpontEmiFile(spontCase)

  open(unit=11, file="dynamics.dat")
  
  ti=titf(1) 
  tf=titf(2)
  dt=titf(3)
  
    
  do i=1,nAtoms
     
     excitFlag(i) = 0; lost(i) = 0.0; C3List(i) = 0.0
     passedTime(i) = 0.0; waitTime(i) = 0.0; GSFlag(i) = 1.0
     phScattEach(i) = 0.0
     
  enddo
  
  do i=1,3
     DeltaC3(i) = 0.0
  enddo
  
  partnerAtom(1) = 2
  partnerAtom(2) = 1
  
  initialCond = firstInitialCond
  phScatt = 0.0
  currentTime = 0.0
  
  alphas = [alpha_GS, alpha(1)] ! for now, a single beam
  
  do while (currentTime < tf)
     
     sols = initialCond
     
     do i=1,nAtoms

        C3List(i) = C3List(i) * (1-GSFlag(i))*GSFlag(partnerAtom(i)) &
             + C3List(i) * (1-GSFLag(partnerAtom(i)))*GSFlag(i)
        
        lost(i) = checkLost(initialCond(i,:), P, w0, C3zero) ! check if atom "i" is lost
        
        ! print*, i, partnerAtom(i)
        
        sols(i,:) = timeEvolution(sols(i,:), sols(partnerAtom(i),:), dt, P, w0, &
             alphas, C3List(i), [GSFlag(i), GSFlag(partnerAtom(i))])
        
        passedTime(i) = passedTime(i) + dt*excitFlag(i)
        
        ! print*, PassedTime(1), GSFlag(1), excitFlag(1), waitTime(1), phScatt
        
     enddo
     
     currentTime = currentTime + dt
     
     initialCond = sols
     
     ! write(11,*) initialCond(1,:)
     ! write(11,*) GSFlag(2), excitFlag(2), passedTime(2)-waitTime(2)
     
     
     ! Checking the probabilities and availabilities for excitation within all cases (1-body and 2-body) 
     
     do i=1,nBeams
        
        call random_number(auxNum)
        index = nint(1+auxNum)
        
        do j=1,nAtoms
           deltaOneBody(i,j) = getOneBodyDelta(sols(j,:), P, w0, delta(i), &
                alpha(i), lambd(i), absProj(i,:))
           !write(11,*) deltaOneBody(i,1)
        enddo
        
        do k=1,2
           call random_number(auxNum)
           DeltaC3(k) = getDetuningFromC3(C3Vals(k), sols(index,:), sols(partnerAtom(index),:))
           deltaTotal(i,index) = deltaOneBody(i,index)+DeltaC3(k)
           ! auxValC3(k) = dabs(auxNum * deltaTotal(i,index))
           auxValC3(k) = Rscatt(Gammas(i), s0(i), deltaTotal(i,index))*dt / auxNum
        enddo
        
        !index = findloc(auxValC3, minval(auxValC3), 1)
        index = findloc(auxValC3, maxval(auxValC3), 1)
        
        do j=1,nAtoms
           deltaTotal(i,j) = deltaOneBody(i,j) + DeltaC3(index)
           scattProb(j) = Rscatt(Gammas(i), s0(i), deltaTotal(i,j))*dt
           
           call random_number(auxNum)
           auxRatios(j) = scattProb(j) / auxNum * GSFlag(j)
           
           if (auxRatios(j) >= 1) then

              C3List(1) = C3Vals(index)
              C3List(2) = C3Vals(index)
              
              recoil = recoilVel(absProj(i,:), lambd(i), "abs", absProj(i,:)) ! parsing absProj(1,:) on last arg. only for type consistency
              
              sols(j, 4:6) = sols(j, 4:6) + recoil ! just one of the atoms receives the kick
              
              initialCond(j, 4:6) = sols(j, 4:6) ! updates the "global" solution array with the velocity kick on one of the atoms
              
              phScatt = phScatt + 1.0 ! updates the total scattered # of photons by 1
              
              excitFlag(j) = 1.0
              
              ! For how long will the excitation last? 
              
              call random_number(auxNum)
              
              waitTime(j) = -log(1-auxNum)/Gammas(i)
              passedTime(j) = 0.0
              
           endif
        enddo
     enddo

     ! write(11,*) C3List(1), C3List(2)
     
     ! Now, we check if the excited atom(s) will decay, with a subsequent change on its/their flags and random kick.
     
     do i=1,nAtoms

        ! see if the atom(s) stayed the required time in |e>
        decayCond = nint(passedTime(i) - waitTime(i) + 0.5 - dt)
        
        call random_number(auxNum)
        
        externalSpontDir = spontEmissionArray(floor(arraySize*auxNum),:)
        recoil = recoilVel(absProj(1,:), lambd(1), "emi", externalSpontDir)

        ! Quick regard for the spontaneous emission:
        ! 1. If the atoms were not "sharing" the excitation: normal recoil
        !                                                    on one atom
        ! 2. If the atoms were "sharing" the excitation: recoil acts on
        !                                                both atoms without
        !                                                changing the CoM vel.
        
        sols(i,4:6) = sols(i,4:6) + recoil * decayCond * excitFlag(i) * excitFlag(partnerAtom(i))
        sols(i,4:6) = sols(i,4:6) + recoll/2 * decayCond * excitFlag(i) * GSFlag(partnerAtom(i))
        sols(partnerAtom(i),4:6) = sols(partnerAtom(i),4:6) &
             - recoil/2 * decayCond * excitFlag(i) * GSFlag(partnerAtom(i))
                
        passedTime(i) = passedTime(i) * (1-decayCond) ! zero the "passedTime" flag if atom(s) should decay
        waitTime(i) = waitTime(i) * (1-decayCond) ! zero the "waitTime" flag if atom(s) should decay
        excitFlag(i) = excitFlag(i) * (1-decayCond) ! excitationFlag should become zero if atom "i" decays
        kickFlag(i) = kickFlag(i) * (1-decayCond) ! changing the kickFlag to 0 upon de-excitation
        
        GSFlag(i) = 1 - excitFlag(i) ! correctly set the GSFlag depending on if the atom "i" is excited
        
     enddo
     
     initialCond = sols
     
     lostFlagSum = int(sum(lost))
     
     if (lostFlagSum == nAtoms) then
        currentTime = tf+1
     endif
     
  end do
  
  close(11)
  
  do i=1,nAtoms
     solution(i) = 1-lost(i) ! in the solution, 1 means the atom survived, 0 means it was lost
  enddo
  
  solution(nAtoms+1) = phScatt ! last entry of "solution" array is the number of scattered photons
  
end subroutine create_simulation_loss

   
   



    



    

    
  
  
    
  
