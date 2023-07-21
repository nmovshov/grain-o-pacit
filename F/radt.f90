module radt
!----------------------------------------------------------------------
! This module contains the subroutine that calculates the dust's
! contribution to the opacity of the atmosphere. The Rosseland mean
! opacity and the optical depth are calculated for each layer. This
! subroutine is based on Morris's RADT.FOR and uses Cuzzie's MIE
! subroutine for the computation of extinction efficiencies. Values for
! the refractive index of the grain material are taken from a table.
!----------------------------------------------------------------------

use globals
implicit none

private mie, dBdT, interp

CONTAINS

subroutine radt1()
! ------------  Dust contribution to opacity -----------
implicit none
real(8) :: h, kB, c, T, nupeak, numin, numax, nr, ni, invkap, normal
integer :: z,j,k
integer, parameter :: MNU=21
real(8) :: nu(MNU)=0, lambda(MNU)=0
real(8) :: crosec(MNU,params_nBins), Kay(MNU,params_nBins), &
           kappa(MNU), weight(MNU)
real(8), allocatable :: table(:,:)

!**********************************************************
! Some physical constants.
!**********************************************************
  h=6.626d-27  ! Planck's constant [erg s]
  kB=1.381d-16 ! Boltzmann's constant [erg/K]
  c=2.998d10   ! Speed of light [cm/s]
  
!**********************************************************
! Prepare the refractive index table.
!**********************************************************
  do k=1,100
    if (dust_index(k,1)==0) exit
  end do
  allocate (table(k-1,3))
  table=dust_index(1:k-1,:)

!**********************************************************
! For each zone in the atmosphere...
!**********************************************************
  do z=1,atm_nZones
  
!**********************************************************
! Define the range of relevant frequencies/wavelengths.
!**********************************************************
    T=atm_T(z)
    nupeak=5.879d10*T ! Location of peak in Planck's function
    numin=nupeak/10
    numax=nupeak*3.3
    nu=(/(numin+k*(numax-numin)/(MNU-1),k=0,MNU-1)/)
    lambda=c/nu
    
!**********************************************************
! Calculate the specific cross-section for all wave lengths and size
! bins.
!**********************************************************
    do j=1,MNU
      nr=interp(table(:,1),table(:,2),lambda(j))
      ni=interp(table(:,1),table(:,3),lambda(j))
      !nr=1.5; ni=0.01
      do k=1,params_nBins
        crosec(j,k)=mie(nu(j),dust_sizeBin(k),nr,ni)
      end do
    end do

!**********************************************************
! Define the extinction coefficient and opacity.
!**********************************************************
    do k=1,MNU
      Kay(k,:)=crosec(k,:)*dust_nDensity(z,1:params_nBins)
      kappa(k)=sum(Kay(k,:))/atm_D(z)
    end do
  
!**********************************************************
! Calculate Rosseland mean opacity.
!**********************************************************
    weight=dBdT(nu,T)
    
! Integrate inverse kappa weighted by dBdT.
    invkap=kappa(1)**(-1)*weight(1)            !
    do k=2,MNU-1,2                             !
      invkap=invkap+4*kappa(k)**(-1)*weight(k) !
    end do                                     !
    do k=3,MNU-2,2                             ! Simpson-rule
      invkap=invkap+2*kappa(k)**(-1)*weight(k) !
    end do                                     !
    invkap=invkap+kappa(MNU)**(-1)*weight(MNU) !
    invkap=invkap*(nu(2)-nu(1))/3.0            !
    
! Integrate dBdT as the normalization factor.
    normal=dBdT(nu(1),T)                     !
    do k=2,MNU-1,2                           !
      normal=normal+4*dBdT(nu(k),T)          !
    end do                                   !
    do k=3,MNU-2,2                           ! Simpson-rule
      normal=normal+2*dBdT(nu(k),T)          !
    end do                                   !
    normal=normal+dBdT(nu(MNU),T)            !
    normal=normal*(nu(2)-nu(1))/3.0          !
    
    invkap=invkap/normal
    
!**********************************************************
! Save opacity and optical depth of layer.
!**********************************************************
    atm_opacity(z)=invkap**(-1)
    atm_opticalDepth(z)=atm_opacity(z)*atm_D(z)*atm_dR(z)

!**********************************************************
! End loop over layers.
!**********************************************************
  end do

return
end subroutine radt1


function mie(nu,a,nr,ni) result(sigma)
!----------------------------------------------------------------------
! This function approximates the scattering and absorption efficiencies
! and the anisotropy parameter, for a given frequency and particle
! size, and returns the effective cross-section of the particle. This 
! is a rough estimate that circumvents that lengthy Mie theory
! calculations. (This is the subroutine from Cuzzie.)
!----------------------------------------------------------------------
implicit none
real(8), intent(in)  :: nu,a,nr,ni
real(8)              :: sigma
real(8)              :: c, lambda, x, pi, cosbar, Qs, Qa, Qe, Q

  c=2.998e10      ! Speed of light [cm/s]
  lambda=c/nu     ! Wavelength [cm]
  pi=3.141592
  x=2*pi*a/lambda ! Size parameter
  
  if (ni.lt.3) then
    cosbar=0.2
    if (x.gt.2.5) cosbar=0.8
  else
    cosbar=-0.2
    if (x.gt.2.5) cosbar=0.5
  end if
  
  if (x<1.3) then
    Qs=8.0/3.0*x**4*((nr**2-1)/(nr**2+2))**2
  else
    Qs=2*x**2*((nr-1)**2+ni**2)
  end if
  if (Qs>1) Qs=1
  
  Qa=12*x*(2*nr*ni)/((nr**2+ni**2+2)**2+(4*nr**2*ni**2))
  if (Qa>1) Qa=1
  
  Qe=Qs+Qa
  
  Q=Qe*(1-Qs/Qe*cosbar)
  
  sigma=Q*pi*a**2
  
return
end function mie

elemental function dBdT(nu,T)
! The Planck's function's temperature derivative.
real(8) :: dBdT
real(8), intent(in) :: nu,T
real(8) :: h, kB, c

  h=6.626d-27  ! Planck's constant [erg s]
  kB=1.381d-16 ! Boltzmann's constant [erg/K]
  c=2.998d10   ! Speed of light [cm/s]

  dBdT=(2*h**2*c**(-2)*kB**(-1)*T**(-2))*(nu**4)*(exp(h*nu/kB/T))*&
       ((exp(h*nu/kB/T)-1)**(-2))
     
return
end function dBdT

end module radt