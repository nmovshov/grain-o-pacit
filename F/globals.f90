MODULE globals
!----------------------------------------------------------------------
! This module holds variables and routines that need to be visible 
! throughout the program. Also, it contains the oSetup subroutine that
! is responsible for reading from the configuration file, initializing
! all the program parameters, reading the input files, and performing
! all one-time calculations that are needed before the evolution code
! can run.
!----------------------------------------------------------------------

implicit none
real(8) :: params_coreMass, params_envMass, params_gasMass, &
           params_sigmaGas
real(8) :: params_silicateDensity, params_accretionRate, &
           params_dust2gasRatio, params_monomerSize
real(8) :: params_binSpacingParameter, params_maxTimeStep, &
           params_stickCoef
integer :: params_nBins
real(8) :: atm_D(200), atm_Dd(200), atm_dR(200), atm_dV(200), &
           atm_g(200), atm_meanFreePath(200), atm_nDensity(200), &
           atm_R(200), atm_T(200), atm_Tt(200), atm_viscosity(200), &
           atm_vThermal(200), atm_Z(200), atm_opacity(200),&
           atm_opticalDepth(200)
integer :: atm_nZones
real(8) :: dust_index(100,3), dust_kernel(200,100,100), &
           dust_massBin(100), dust_nDensity(200,100), &
           dust_redist(100,100,100), dust_sizeBin(100), &
           dust_vSed(200,100),dust_source(200)
integer :: dust_aBins
real(8) :: debug_massInAtm(200),debug_massInCore,debug_ndCore(100)
logical :: cfg_flxFlg,cfg_cogFlg,cfg_srcFlg,cfg_newFlg, cfg_tesFlg,&
           cfg_msgFlg
character(30) :: cfg_atmFile,cfg_indFile,cfg_distFile,cfg_outFile, &
                 cfg_tesFile,cfg_srcSize
real(8) :: tTimes(10)

CONTAINS

subroutine oSetup(fileName)
!**********************************************************
! Read configuration and free parameters from a file and perform all
! one-time calculations to create model variables.
!**********************************************************
implicit none
character(*) :: fileName ! The configuration file
real(8) :: kBoltz, univG, avogad, earthM, pi, volume, kn, psi, lB, &
           D, del, v, r, P1, P2, x, y(200,2), bigEm, diEm(200)
real(8) :: difuse(200,100), vDust(200,100), delta(200,100),masBin(100)
integer :: i,j,k,nLines,nZones,nBins
character(60) :: cLine,cProp,cVal

!**********************************************************  
! Define some physical constants.
!**********************************************************
  kBoltz = 1.3806d-16 ! Boltzmann's constant [erg/K]
  univG  = 6.673d-8   ! Universal gravity constant [cm^3 g^-1 s^-2]
  avogad = 6.0221d23  ! Avogadro's number [mol^-1]
  earthM = 5.974d27   ! Earth mass [g]
  pi     = 3.14159265 ! Pi

!**********************************************************  
! Initialize all arrays and all flags.
!**********************************************************
  atm_R=0
  atm_Dd=0
  atm_Tt=0
  atm_Z=0
  atm_dR=0
  atm_dV=0
  atm_D=0
  atm_T=0
  atm_g=0
  atm_vThermal=0
  atm_meanFreePath=0
  atm_nDensity=0
  atm_viscosity=0
  atm_opacity=0
  atm_opticalDepth=0
  dust_vSed=0
  dust_sizeBin=0
  dust_massBin=0
  dust_nDensity=0
  dust_kernel=0
  dust_redist=0
  dust_index=0
  dust_source=0
  debug_massInAtm=0
  debug_ndCore=0
  difuse=0
  diEm=0
  vDust=0
  delta=0
  tTimes=0
  cfg_flxFlg=.true.
  cfg_cogFlg=.true.
  cfg_srcFlg=.true.
  cfg_tesFlg=.true.
  cfg_newFlg=.true.
  cfg_msgFlg=.true.


!**********************************************************  
! Read parameters from the configuration file.
!**********************************************************
  open(unit=11,file=fileName,action='read',status='old')
  cLine=''
  cProp=''
  cVal=''
  do ! Skip the header
    read(11,'(a)')cLine
    if (trim(adjustl(cLine))=='--- BEGIN CONFIG. ---') exit
  end do
  do ! Read from file line-by-line
    read(11,'(a)')cLine
    if (trim(adjustl(cLine))=='--- END CONFIG. ---') exit
    if (trim(adjustl(cLine))=='') cycle ! Skip empty lines
    k=index(cLine,'=')
    if (k==0) cycle ! Skip description lines
    cProp=trim(adjustl(cLine(1:k-1)))
    select case (cProp) ! Assign values to the appropriate variable
    
    case ('Core mass') ! Core mass
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_coreMass=x*earthM
            
    case ('Envelope mass') ! Envelope mass
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_envMass=x*earthM      
      
    case ('Maximum time step') ! Maximum allowed time step
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_maxTimeStep=x
      
    case ('Target time(s)') ! Target time(s)
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*,end=91)tTimes
91    close(12)
            
    case ('Grain matter density') ! Density of dust matter
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_silicateDensity=x
      
    case ('Sticking coefficient') ! Sticking coefficient
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_stickCoef=x
      
    case ('Number of size bins') ! Number of size bins
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_nBins=int(x)
      
    case ('Initial distribution file') ! Filename or 'new'
      cVal=trim(adjustl(cLine(k+1:)))
      if (cVal=='new') then
        cfg_newFlg=.true.
        cfg_distFile=''
      else
        cfg_newFlg=.false.
        cfg_distFile=cVal
      end if

    case ('Planetesimals source file') ! Filename or 'none'
      cVal=trim(adjustl(cLine(k+1:)))
      if (cVal=='none') then
        cfg_tesFlg=.false.
        cfg_tesFile=''
      else
        cfg_tesFlg=.true.
        cfg_tesFile=cVal
      end if

    case ('Source size') ! 'small' or 'equal'
      cVal=trim(adjustl(cLine(k+1:)))
      cfg_srcSize=cVal

    case ('Atmosphere file') ! Filename
      cVal=trim(adjustl(cLine(k+1:)))
      cfg_atmFile=cVal
           
    case ('Refractive index table') ! Filename
      cVal=trim(adjustl(cLine(k+1:)))
      cfg_indFile=cVal

    case ('Monomer size (smallest bin)') ! Size of smallest bin/grain
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_monomerSize=x
      
    case ('Solids accretion rate') !  Solids accretion rate
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_accretionRate=x
      
    case ('Bin spacing parameter') ! Bin spacing parameter
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_binSpacingParameter=x
    
    case ('Molecular weight') ! Molecular weight of atmospheric gas
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_gasMass=x/avogad
    
    case ('Collision cross-section') ! Collision cross-section of a gas
      cVal=trim(adjustl(cLine(k+1:)))! molecule
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_sigmaGas=x
    
    case ('Dust-to-gas ratio') ! Solids-to-gas ratio
      cVal=trim(adjustl(cLine(k+1:)))
      open(12,status='scratch')
      write(12,*)cVal
      rewind(12)
      read(12,*)x
      close(12)
      params_dust2gasRatio=x
      
    case ('Coagulation') ! Coagulation process flag
      cVal=trim(adjustl(cLine(k+1:)))
      if (cVal=='off') then
        cfg_cogFlg=.false.
      else if (cVal=='on') then
        cfg_cogFlg=.true.
      else
        write(*,*)'Error! Wrong value for property: coagulation'
        stop
      end if
      
    case ('Sedimentation') ! Sedimentation process flag
      cVal=trim(adjustl(cLine(k+1:)))
      if (cVal=='off') then
        cfg_flxFlg=.false.
      else if (cVal=='on') then
        cfg_flxFlg=.true.
      else
        write(*,*)'Error! Wrong value for property: sedimentation'
        stop
      end if
      
    case ('Source term') ! Source process flag
      cVal=trim(adjustl(cLine(k+1:)))
      if (cVal=='off') then
        cfg_srcFlg=.false.
      else if (cVal=='on') then
        cfg_srcFlg=.true.
      else
        write(*,*)'Error! Wrong value for property: source'
        stop
      end if
          
    case ('Progress messages') ! Progress messages flag
      cVal=trim(adjustl(cLine(k+1:)))
      if (cVal=='off') then
        cfg_msgFlg=.false.
      else if (cVal=='on') then
        cfg_msgFlg=.true.
      else
        write(*,*)'Error! Wrong value for property: progress messages'
        stop
      end if

    case default
      write(*,*)'Error! Unknown configuration parameter: ',cProp
      stop
      
    end select
    
  end do
  close(11)

!**********************************************************
! Read the model atmosphere.
!**********************************************************
  open(unit=11,file=cfg_atmFile,action='read',status='old')
  do ! Skip the header
    read(11,'(a)')cLine
    if (trim(adjustl(cLine))=='--- BEGIN DATA ---') exit
  end do
  do k=1,200
    read(11,*,end=97) atm_R(k),atm_Dd(k),x,atm_Tt(k),x,x
  end do
97 continue
  close(11)
  do k=1,200 ! Determine actual number of lines in the model.
    if (atm_R(k).eq.0) exit
  end do
  nLines=k-1

!**********************************************************
! Calculate derived properties.
!**********************************************************
  atm_nZones=nLines-1 !*
!* There are length(R)-1 'zones'. Zone(i) is the volume
! between R(i) and R(i+1). Any property of zone(i), e.g. temperature,
! will be in the ith place of the appropriate array, and averaged
! from the values on the boundaries.
  nZones=atm_nZones
  do k=1,nZones
    atm_Z(k)=(atm_R(k)+atm_R(k+1))/2
    atm_dR(k)=atm_R(k)-atm_R(k+1)
    atm_dV(k)=4*pi/3*(atm_R(k)**3-atm_R(k+1)**3)
    atm_D(k)=(atm_Dd(k)+atm_Dd(k+1))/2
    atm_T(k)=(atm_Tt(k)+atm_Tt(k+1))/2
    diEm=atm_D*atm_dV
    bigEm=params_coreMass+params_envMass-sum(diEm)
    atm_g(k)=univG*(bigEm+sum(diEm(k+1:)))/atm_Z(k)**2
    atm_vThermal(k)=sqrt((8*kBoltz*atm_T(k))/(pi*params_gasMass))
    atm_nDensity(k)=atm_D(k)/params_gasMass
    atm_meanFreePath(k)=1/(sqrt(2.d0)*params_sigmaGas*atm_nDensity(k))
    atm_viscosity(k)=8.6d-6*sqrt(atm_T(k))
  end do

!**********************************************************
! Create the size bins and corresponding mass bins.
!**********************************************************
  nBins=params_nBins
  do k=1,nBins+1
! * Note. The +1 is there to facilitate coding of the mass
!   redistribution coefficients. The nBins+1 elements are never
!   used.
    x=params_binSpacingParameter
    dust_sizeBin(k)=params_monomerSize*2**(k/x)
    volume=4*pi/3*dust_sizeBin(k)**3
    dust_massBin(k)=volume*params_silicateDensity
  end do
      
!**********************************************************
! Calculate the largest allowed active bin. *
! * We use the last few bins as "recycle bins", to make sure no mass
!   is lost due to coagulation into too large particles. We need the
!   largest "active" bin to be half the mass of the largest bin.
!   (With a spacing parameter of 3 aBins=nBins-1.)
!**********************************************************
  do k=1,nBins
    if (2*dust_massBin(k)-dust_massBin(nBins).gt.&
        1e-3*dust_massBin(nBins)) exit
  end do
  dust_aBins=k-1

!**********************************************************
! Calculate sedimentation speed for each size bin in each zone.
!**********************************************************
  do j=1,nZones
    do k=1,nBins
      kn=atm_meanFreePath(j)/dust_sizeBin(k)
      psi=1+kn*(1.249+0.42*exp(-0.87/kn))
      dust_vSed(j,k)=(dust_massBin(k)*atm_g(j)*psi)/&
                     (6*pi*dust_sizeBin(k)*atm_viscosity(j))+&
                     params_accretionRate/params_dust2gasRatio/atm_D(j)
    end do
  end do

!**********************************************************
! Initialize the nd array. *
! * This is the array that holds the number density of dust particles
!   in each size bin and each zone in the atmosphere.
!*********************************************************
  if (cfg_newFlg) then
    do k=1,nZones
      dust_nDensity(k,1)=params_dust2gasRatio*atm_D(k)/dust_massBin(1)
    end do
  else
    open(unit=11,file=cfg_distFile,action='read',status='old')
    do ! Skip the header
      read(11,'(a)')cLine
      if (trim(adjustl(cLine))=='--- BEGIN DATA ---') exit
    end do
    do j=1,nZones
      read(11,*)(dust_nDensity(j,k),k=1,nBins)
    end do
    close(11)
  end if

!**********************************************************
! Create the collision kernel. *
! * The 'kernel' array holds the probabilities of collision, as a 
!   result of Brownian motion and 'overtaking'. kernel(i,j,k) is the
!   probability of collision between a j-bin particle and a k-bin
!   particle in the ith zone. This is unnormalized probability, i.e.,
!   kernel*density*density is the number of collisions per unit volume
!   per unit time.
!   (Note to myself: this is a semi-black-box calculation for me.
!   This is section 2.2 in podolak03.pdf.)
!**********************************************************
  do j=1,nZones
    do k=1,nBins
      difuse(j,k)=(3*kBoltz*atm_T(j))/&
                  (4*pi*dust_sizeBin(k)**2*atm_nDensity(j)*&
                  params_gasMass*atm_vThermal(j))
      vDust(j,k)=atm_vThermal(j)*sqrt(params_gasMass/dust_massBin(k))
      lB=(8*difuse(j,k))/(pi*vDust(j,k))
      delta(j,k)=(sqrt(2.d0)/(6*dust_sizeBin(k)*lB))*&
                 ((2*dust_sizeBin(k)+lB)**3-&
                 (4*dust_sizeBin(k)**2+lB**2)**(3./2.))-&
                 2*dust_sizeBin(k)
    end do
  end do
  do i=1,nZones
    do j=1,nBins
      do k=1,nBins
        D=(difuse(i,j)+difuse(i,k))/2
        del=sqrt(delta(i,j)**2+delta(i,k)**2)
        v=sqrt(vDust(i,j)**2+vDust(i,k)**2)
        r=(dust_sizeBin(j)+dust_sizeBin(k))/2
        P1=8*pi*r*D/(r/(r+del/2)+(4*D)/(r*v))
        P2=4*pi*r**2*abs(dust_vSed(i,j)-dust_vSed(i,k))
        dust_kernel(i,j,k)=params_stickCoef*(P1+P2)
      end do
    end do
  end do
      
!**********************************************************
! Create coefficients of mass redistribution. *
! * The 'redist' array holds the probability of mass redistribution
!   options after collision. redist(i,j,k) is the probability that
!   after a collision of an i-bin particle with a j-bin particle, a
!   k-bin particle will be created.
!**********************************************************
  masBin=dust_massBin
  if (2*masBin(1).le.masBin(2)) then
    dust_redist(1,1,1)=(masBin(2)-2*masBin(1))/(masBin(2)-masBin(1))
  end if
  do i=1,nBins
    do j=1,nBins
      do k=2,nBins
        if ((masBin(i)+masBin(j).gt.masBin(k-1)).and.&
              (masBin(i)+masBin(j).le.masBin(k))) then
          dust_redist(i,j,k)=(masBin(i)+masBin(j)-masBin(k-1))/&
                          (masBin(k)-masBin(k-1))
        end if
        if ((masBin(i)+masBin(j).gt.masBin(k)).and.&
             (masBin(i)+masBin(j).le.masBin(k+1))) then
          dust_redist(i,j,k)=(masBin(k+1)-masBin(i)-masBin(j))/&
                          (masBin(k+1)-masBin(k))
        end if
      end do
    end do
  end do

!**********************************************************
! Read in a table of refractive indices for different wavelengths
!**********************************************************
  open(unit=11,file=cfg_indFile,action='read',status='old')
  do ! Skip the header
    read(11,'(a)')cLine
    if (trim(adjustl(cLine))=='--- BEGIN DATA ---') exit
  end do
  do j=1,size(dust_index,1)
    read(11,*,end=89)dust_index(j,1),x,x,dust_index(j,2),&
                     dust_index(j,3)
  end do
89 continue
  close(11)
  dust_index(:,1)=dust_index(:,1)*1d-4 ! Table values are in microns

!**********************************************************
! Read in a source of solids from planetesimal ablation
!**********************************************************
  if (cfg_tesFlg) then
    open(unit=11,file=cfg_tesFile,action='read',status='old')
  do ! Skip the header
    read(11,'(a)')cLine
    if (trim(adjustl(cLine))=='--- BEGIN DATA ---') exit
  end do
    y=0
    do k=1,200
      read(11,*,end=88)y(k,:)
    end do
88  i=k-1
    do k=1,atm_nZones
      dust_source(k)=interp(y(1:i,1),y(1:i,2),atm_Z(k))
    end do
    close(11)
    dust_source(1)=dust_source(1)+params_accretionRate/atm_dR(1)
  else
    dust_source(1)=params_accretionRate/atm_dR(1)
  end if

  return
end subroutine oSetup

function interp(xx,yy,x) result(y)
! Linear interpolation, XX monotonic. (From Numerical Recipes.)
implicit none
real(8), intent(in) :: xx(:),yy(:),x
real(8)             :: y
real(8)             :: A,B
integer             :: j,jl,jm,ju,n

  jl=0
  n=size(xx)
  ju=n+1
  do while (ju-jl.gt.1)
    jm=(ju+jl)/2
    if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
      jl=jm
    else
      ju=jm
    endif
  end do
  if(x.eq.xx(1)) then
    j=1
  else if(x.eq.xx(n)) then
    j=n-1
  else
    j=jl
  end if
  
  if (j==0) then
    y=yy(1)
  else if (j==n) then
    y=yy(n)
  else
    A=(xx(j+1)-x)/(xx(j+1)-xx(j))
    B=1-A
    y=A*yy(j)+B*yy(j+1)
  end if

return
end function interp

subroutine debugMass()
! Calculate how much mass is in each layer, and how much is in the core
implicit none
real(8) :: x
integer :: k,j
  
  do j=1,atm_nZones
    x=0
    do k=1,params_nBins
      x=x+dust_massBin(k)*dust_nDensity(j,k)*atm_dV(j)
    end do
    debug_massInAtm(j)=x
  end do
  x=0
  do k=1,params_nBins
    x=x+debug_ndCore(k)*dust_massBin(k)
  end do
  debug_massInCore=x

return
end subroutine debugMass



END MODULE globals