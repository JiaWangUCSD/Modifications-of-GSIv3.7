module intreflmod

use m_obsNode,  only: obsNode
use m_reflNode, only: reflNode
use m_reflNode, only: reflNode_typecast
use m_reflNode, only: reflNode_nextcast

implicit none

PRIVATE
PUBLIC intrefl

interface intrefl; module procedure &
          intrefl_
end interface

contains

subroutine intrefl_(reflhead,rval,sval)
! abstract: apply observation operator for radar winds
!             with nonlinear qc operator
! usage: call intw(ru,rv,su,sv)
!   input argument list:
!     reflhead   - obs type pointer to obs structure     
!     su       - current u solution increment 
!     sv       - current v solution increment 
!     ru
!     rv
!
!   output argument list:
!     ru       - u results from observation operator 
!     rv       - v results from observation operator 
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: half,one,tiny_r_kind,cg_term,r3600
  use obsmod, only: lsaveobsens,l_do_adjoint,luse_obsdiag
  use qcmod, only: nlnqc_iter,varqc_iter
  use qcmod, only: vqc
  use jfunc, only: jiter
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_4dvar, only: ladtest_obs
  implicit none

! Declare passed variables
  class(obsNode), pointer, intent(in   ) :: reflhead
  type(gsi_bundle),        intent(in   ) :: sval
  type(gsi_bundle),        intent(inout) :: rval

! Declare local varibles
  logical include_w
  integer(i_kind) j1,j2,j3,j4,j5,j6,j7,j8,ier,istatus
! real(r_kind) penalty
  real(r_kind) val,valrefl,w1,w2,w3,w4,w5,w6,w7,w8
  real(r_kind) cg_refl,p0,grad,wnotgross,wgross,pg_refl
  real(r_kind),pointer,dimension(:) :: srefl
  real(r_kind),pointer,dimension(:) :: rrefl
  type(reflNode), pointer :: reflptr

!  If no refl data return
  if(.not. associated(reflhead))return

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'refl',srefl,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'refl',rrefl,istatus);ier=istatus+ier
  if(ier/=0)return


  reflptr => reflNode_typecast(reflhead)
  do while (associated(reflptr))
     j1=reflptr%ij(1)
     j2=reflptr%ij(2)
     j3=reflptr%ij(3)
     j4=reflptr%ij(4)
     j5=reflptr%ij(5)
     j6=reflptr%ij(6)
     j7=reflptr%ij(7)
     j8=reflptr%ij(8)
     w1=reflptr%wij(1)
     w2=reflptr%wij(2)
     w3=reflptr%wij(3)
     w4=reflptr%wij(4)
     w5=reflptr%wij(5)
     w6=reflptr%wij(6)
     w7=reflptr%wij(7)
     w8=reflptr%wij(8)

!    Forward model (Tangent Linear; TL)
!    TLVr  =  TLu*costilt*cosazm  +  TLv*costilt*sinazm  +  TLw*sintilt
     val=(w1*srefl(j1)+ w2*srefl(j2)+ w3*srefl(j3)+ w4*srefl(j4)+ w5*srefl(j5)+    &
          w6*srefl(j6)+ w7*srefl(j7)+ w8*srefl(j8))


     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*reflptr%raterr2*reflptr%err2
           reflptr%diags%obssen(jiter) = grad
        else
           if (reflptr%luse) reflptr%diags%tldepart(jiter)=val
        endif
     endif

     ! l_do_adjoint=.true.: apply H^T when in int routines
     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           ! ladtest_obs: Run adjoint test for obervation
           ! val (RHS): HX
           ! val (LHS): HX-y, y=O-B
           if( .not. ladtest_obs ) val=val-reflptr%res

!          gradient of nonlinear operator
           ! vqc: when true, use ECMWF's non linear QC
           ! nlnqc_iter: logical flag (T=nonlinear qc on, F=nonlinear qc off) for iteration
           ! for reflectivity: reflptr%pg=0.0 in "global_convinfo.txt"
           if (vqc .and. nlnqc_iter .and. reflptr%pg > tiny_r_kind .and. &
                                reflptr%b  > tiny_r_kind) then
              pg_refl=reflptr%pg*varqc_iter
              cg_refl=cg_term/reflptr%b
              wnotgross= one-pg_refl
              wgross = pg_refl*cg_refl/wnotgross
              p0   = wgross/(wgross+exp(-half*reflptr%err2*val**2))
              val = val*(one-p0)
           endif

           if( ladtest_obs)  then
              grad = val
           else
              ! grad = (HX-y)/E^2
              grad = val*reflptr%raterr2*reflptr%err2
           end if

        endif

!       Adjoint (AD)
        valrefl=1.0_r_kind*grad  ! ADVr_u = costilt*cosazm*ADVr
        rrefl(j1)=rrefl(j1)+w1*valrefl                 ! ADu = ADu + ADVr_u
        rrefl(j2)=rrefl(j2)+w2*valrefl
        rrefl(j3)=rrefl(j3)+w3*valrefl
        rrefl(j4)=rrefl(j4)+w4*valrefl
        rrefl(j5)=rrefl(j5)+w5*valrefl
        rrefl(j6)=rrefl(j6)+w6*valrefl
        rrefl(j7)=rrefl(j7)+w7*valrefl
        rrefl(j8)=rrefl(j8)+w8*valrefl

     endif

     reflptr => reflNode_nextcast(reflptr)
  end do
  return
end subroutine intrefl_

end module intreflmod
