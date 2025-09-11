module stpreflmod

! module:   stpreflmod    module for stprefl and its tangent linear stprefl_tl

implicit none

PRIVATE
PUBLIC stprefl

contains

subroutine stprefl(reflhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stprefl       calculate penalty and contribution to
!                            stepsize with nonlinear qc added.
!
!   input argument list:
!     reflhead
!     ru       - search direction for u
!     rv       - search direction for v
!     rw       - search direction for w
!     su       - analysis increment for u
!     sv       - analysis increment for v
!     sw       - analysis increment for w
!     sges     - step size estimates (nstep)
!     nstep    - number of step sizes (== 0 means use outer iteration value)
!
!   output argument list     - output for step size calculation
!     out(1:nstep)   - penalty from radar winds sges(1:nstep)
!
!$$$
  use kinds, only: r_kind,i_kind,r_quad
  use qcmod, only: nlnqc_iter,varqc_iter
  use qcmod, only: vqc
  use constants, only: half,one,two,tiny_r_kind,cg_term,zero_quad,r3600
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use m_obsNode,  only: obsNode
  use m_reflNode, only: reflNode
  use m_reflNode, only: reflNode_typecast
  use m_reflNode, only: reflNode_nextcast

  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: reflhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables
  logical include_w
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4,j5,j6,j7,j8,kk
  real(r_kind) valrefl,facrefl,w1,w2,w3,w4,w5,w6,w7,w8
  real(r_kind) cg_refl,r,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep))::pen
  real(r_kind) pg_refl
  real(r_kind),pointer,dimension(:) :: srefl
  real(r_kind),pointer,dimension(:) :: rrefl
  type(reflNode), pointer :: reflptr

  out=zero_quad

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
     if(reflptr%luse)then
        if(nstep > 0)then
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

!          Gradient
           valrefl=(w1*rrefl(j1)+ w2*rrefl(j2)+ w3*rrefl(j3)+ &
                    w4*rrefl(j4)+ w5*rrefl(j5)+ w6*rrefl(j6)+ &
                    w7*rrefl(j7)+ w8*rrefl(j8))

!          Gradient - residual
           facrefl=(w1* srefl(j1)+w2* srefl(j2)+w3* srefl(j3)+w4* srefl(j4)+w5* srefl(j5)+  &
                   w6* srefl(j6)+w7* srefl(j7)+w8* srefl(j8))
           facrefl=facrefl-reflptr%res

           do kk=1,nstep
              r=facrefl+sges(kk)*valrefl
              pen(kk)=r*r*reflptr%err2
           end do
        else
           pen(1)=reflptr%res*reflptr%res*reflptr%err2
        end if

!  Modify penalty term if nonlinear QC
        if (vqc .and. nlnqc_iter .and. reflptr%pg > tiny_r_kind .and.  &
                             reflptr%b  > tiny_r_kind) then
           pg_refl=reflptr%pg*varqc_iter
           cg_refl=cg_term/reflptr%b
           wnotgross= one-pg_refl
           wgross = pg_refl*cg_refl/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif

        out(1) = out(1)+pen(1)*reflptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*reflptr%raterr2
        end do
     end if

     reflptr => reflNode_nextcast(reflptr)

  end do
  return
end subroutine stprefl

end module stpreflmod
