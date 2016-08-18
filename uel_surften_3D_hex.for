!************************************************************************
! User element (UEL) for large-deformation viscoelastic deformation 
!  in three-dimensions with surface tension
!
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-3).
!
! Material behavior is Neo-Hookean or Gent rubber elasticity 
!  with a finite-deformation Maxwell element in parallel.
!
!                 |
!           -------------
!           |           |  Maxwell element:
! Neo-Hooke |           |
!  or Gent  |           \
!  spring:  \           /    Hencky spring:Gneq
!   G0      /           |    
!   Kbulk   \           |
!  (Im)     /           |    Linearly viscous
!           |         | - |   Mises dashpot:eta
!           |         |___|
!           |           |
!           -------------
!                 |
! 
! This subroutine is for a three-dimensional 8-node isoparametric
!  brick element as shown below with 8pt (full) integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
! Surface tension contributions are included in this element.  
!  Based on our convention, the face on which the surface tension
!  acts is the "label", i.e.
!  - U1,U2,U3,U4,U5,U6 refer to surface tension acting
!    on faces 1,2,3,4,5,6, respectively,
!  Mechanical, traction- and pressure-type boundary conditions 
!  may be applied to the dummy mesh using the Abaqus built-in 
!  commands *Dload or *Dsload.
!
!
!  8-node     8-----------7
!  brick     /|          /|       xi_3
!           / |         / |       
!          5-----------6  |       |     xi_2
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi_1
!          1-----------2        origin at cube center
!
!
! David L. Henann, May 2016
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=8,Type=U1,Iproperties=1,Properties=5,Coordinates=3,Variables=72,Unsymm
!  1,2,3
!
!     State Variables
!       jj = 0
!       do k = 1,nIntPt
!          svars(1+jj) = Fv(1,1) ---- Fv(1,1) at integ pt k
!          svars(2+jj) = Fv(1,2) ---- Fv(1,2) at integ pt k
!          svars(3+jj) = Fv(1,3) ---- Fv(1,3) at integ pt k
!          svars(4+jj) = Fv(2,1) ---- Fv(2,1) at integ pt k
!          svars(5+jj) = Fv(2,2) ---- Fv(2,2) at integ pt k
!          svars(6+jj) = Fv(2,3) ---- Fv(2,3) at integ pt k
!          svars(7+jj) = Fv(3,1) ---- Fv(3,1) at integ pt k
!          svars(8+jj) = Fv(3,2) ---- Fv(3,2) at integ pt k
!          svars(9+jj) = Fv(3,3) ---- Fv(3,3) at integ pt k
!          jj = jj + nlSdv
!       end loop over k
!
! In the input file, set the parameter matflag=1 for Neo-Hookean material 
!  behavior or matflag=2 for Gent.
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     G0      = props(1)  ! Ground-state shear modulus
!     Kbulk   = props(2)  ! Bulk modulus
!     Im      = props(3)  ! Limited chain extensibility parameter (Gent only)
!     Gneq    = props(4)  ! Nonequilibrium (Maxwell element) stiffness
!     eta     = props(5)  ! Maxwell element viscosity
!     matflag = jprops(1) ! Neo-Hookean=1, Gent=2
!
!************************************************************************

      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     +     props,nprops,coords,mcrd,nnode,uall,duall,vel,accn,jtype,
     +     time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     +     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,
     +     njprop,period)

      implicit none
      !
      ! variables defined in uel, passed back to Abaqus
      !
      real*8 rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
      !
      ! variables passed into UEL
      !
      integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      real*8 props(*),coords(mcrd,nnode),uall(ndofel),duall(mlvarx,*),
     +  vel(ndofel),accn(ndofel),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      !
      ! variables defined and used in the UEL
      !
      integer i,j,k,jj,A11,B11,A12,B12,nDim,nSdv,nInt,nIntSurf,nIntPt,
     +  intpt,stat,matflag,face,nNodeSurf,nIntPtSurf,intptSurf
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=3)      ! number of spatial dimensions, do not change
      parameter(nSdv=9)      ! number of state variables per integ pt, do not change
      parameter(nNodeSurf=4) ! number of nodes on a face, do not change
      parameter(nInt=8)      ! number of integration points, do not change
      parameter(nIntSurf=4)  ! number of surface integration points, do not change
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 u(nNode,3),du(nNode,ndofel),uOld(nNode,ndofel),
     +  coordsC(mcrd,nNode),Ru(3*nNode,1),Kuu(3*nNode,3*nNode),
     +  Iden(3,3),xi(nInt,3),w(nInt),sh0(nNode),sh(nNode),dsh0(nNode,3),
     +  dshC0(nNode,3),dsh(nNode,3),dshC(nNode,3),dshxi(nNode,3),
     +  detMapJ0,detMapJ0C,detMapJ,detMapJC,Fc_tau(3,3),detFc_tau,
     +  F_tau(3,3),detF_tau,T_tau(3,3),Fv_t(3,3),Fv_tau(3,3),
     +  SpTanMod(3,3,3,3),Smat(6,1),Bmat(6,3*nNode),Gmat(9,3*nNode),
     +  G0mat(9,3*nNode),Amat(9,9),Qmat(9,9),body(3),
     +  BodyForceRes(3*nNode,1),gamma,Rsurf(3*nNodeSurf,1),
     +  Ksurf(3*nNodeSurf,3*nNodeSurf),coordsCSurf(3*nNodeSurf,1),
     +  xiSurf(nIntSurf,2),wSurf(nIntSurf),shSurf(nNodeSurf),
     +  dshxiSurf(nNodeSurf,2),dshxiMat1(3,3*nNodeSurf),
     +  dshxiMat2(3,3*nNodeSurf),E,F,G,H2,H,dE(1,3*nNodeSurf),
     +  dG(1,3*nNodeSurf),dF(1,3*nNodeSurf),dH(1,3*nNodeSurf)
      !
      integer nListSurf(nNodeSurf)
      !
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Check the procedure type; this should be a 
      !  *Static step, which is either 1 or 2
      !
      if((lflags(1).eq.1).or.(lflags(1).eq.2)) then
         !
         ! correct procedure specified
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear perturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear perturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear perturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return


      ! Get flag for material behavior
      !
      matflag = jprops(1)


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Kuu = zero
      Energy = zero


      ! Obtain nodal displacements
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         end do
      end do


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33, 3277-3296.
      !
      ! Obtain shape functions and their local gradients at the element
      !  centroid, that means xi_1=xi_2=x_3=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centroid
      !  at the end of the increment for use in the `F-bar' method
      !  The subscript tau denotes the time at the end of the increment.
      !
      Fc_tau = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
            enddo
         enddo
      enddo
      call mdet(Fc_tau,detFc_tau)
      !
      ! With the deformation gradient known at the element centroid
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.8) then
         if(nInt.eq.8) then
            call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
         else
            write(*,*) 'Incorrect number of int points: nInt.ne.8'
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! This is the first increment of the first step.
            !  Give initial conditions.
            !
            Fv_t     = Iden
            !
         else
            !
            ! This is not the first increment; read old values.
            !
            Fv_t(1,1) = svars(1+jj)
            Fv_t(1,2) = svars(2+jj)
            Fv_t(1,3) = svars(3+jj)
            Fv_t(2,1) = svars(4+jj)
            Fv_t(2,2) = svars(5+jj)
            Fv_t(2,3) = svars(6+jj)
            Fv_t(3,1) = svars(7+jj)
            Fv_t(3,2) = svars(8+jj)
            Fv_t(3,3) = svars(9+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain the deformation gradient at this integration point.
         !  The subscript tau denotes the time at the end of the increment.
         !
         F_tau = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         
         
         ! Modify the deformation gradient for the `F-bar' method.
         !
         call mdet(F_tau,detF_tau)
         F_tau = ((detFc_tau/detF_tau)**third)*F_tau


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive update at this integ. point
         !
         if (matflag.eq.1) then
            call NeoHookean(props,nprops,dtime,F_tau,Fv_t,
     +           T_tau,Fv_tau,SpTanMod,stat)
         elseif (matflag.eq.2) then
            call Gent(props,nprops,dtime,F_tau,Fv_t,
     +           T_tau,Fv_tau,SpTanMod,stat)
         else
            write(*,*) 'Invalid matflag: matflag.ne.1 or 2'
            call xit
         endif
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment.
         !
         svars(1+jj)  = Fv_tau(1,1)
         svars(2+jj)  = Fv_tau(1,2)
         svars(3+jj)  = Fv_tau(1,3)
         svars(4+jj)  = Fv_tau(2,1)
         svars(5+jj)  = Fv_tau(2,2)
         svars(6+jj)  = Fv_tau(2,3)
         svars(7+jj)  = Fv_tau(3,1)
         svars(8+jj)  = Fv_tau(3,2)
         svars(9+jj)  = Fv_tau(3,3)
         jj = jj + nSdv ! setup for the next intPt


         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         !
         Bmat = zero
         do k=1,nNode
            Bmat(1,1+nDim*(k-1)) = dshC(k,1)
            Bmat(2,2+nDim*(k-1)) = dshC(k,2)
            Bmat(3,3+nDim*(k-1)) = dshC(k,3)
            Bmat(4,1+nDim*(k-1)) = dshC(k,2)
            Bmat(4,2+nDim*(k-1)) = dshC(k,1)
            Bmat(5,2+nDim*(k-1)) = dshC(k,3)
            Bmat(5,3+nDim*(k-1)) = dshC(k,2)
            Bmat(6,1+nDim*(k-1)) = dshC(k,3)
            Bmat(6,3+nDim*(k-1)) = dshC(k,1)
         enddo
         !
         body = zero ! The body force vector may be specified here
         !
         BodyForceRes = zero
         do k=1,nNode
            BodyForceRes(1+nDim*(k-1),1) = sh(k)*body(1)
            BodyForceRes(2+nDim*(k-1),1) = sh(k)*body(2)
            BodyForceRes(3+nDim*(k-1),1) = sh(k)*body(3)
         enddo
         !
         Ru = Ru - matmul(transpose(Bmat),Smat)*detmapJC*w(intpt)
     +        + BodyForceRes*detmapJC*w(intpt)


         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do k=1,nNode
            Gmat(1,1+nDim*(k-1)) = dshC(k,1)
            Gmat(2,2+nDim*(k-1)) = dshC(k,1)
            Gmat(3,3+nDim*(k-1)) = dshC(k,1)
            Gmat(4,1+nDim*(k-1)) = dshC(k,2)
            Gmat(5,2+nDim*(k-1)) = dshC(k,2)
            Gmat(6,3+nDim*(k-1)) = dshC(k,2)
            Gmat(7,1+nDim*(k-1)) = dshC(k,3)
            Gmat(8,2+nDim*(k-1)) = dshC(k,3)
            Gmat(9,3+nDim*(k-1)) = dshC(k,3)
         enddo
         !
         G0mat = zero
         do k=1,nNode
            G0mat(1,1+nDim*(k-1)) = dshC0(k,1)
            G0mat(2,2+nDim*(k-1)) = dshC0(k,1)
            G0mat(3,3+nDim*(k-1)) = dshC0(k,1)
            G0mat(4,1+nDim*(k-1)) = dshC0(k,2)
            G0mat(5,2+nDim*(k-1)) = dshC0(k,2)
            G0mat(6,3+nDim*(k-1)) = dshC0(k,2)
            G0mat(7,1+nDim*(k-1)) = dshC0(k,3)
            G0mat(8,2+nDim*(k-1)) = dshC0(k,3)
            G0mat(9,3+nDim*(k-1)) = dshC0(k,3)
         enddo
         !
         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,3,1)
         Amat(1,4) = SpTanMod(1,1,1,2)
         Amat(1,5) = SpTanMod(1,1,2,2)
         Amat(1,6) = SpTanMod(1,1,3,2)
         Amat(1,7) = SpTanMod(1,1,1,3)
         Amat(1,8) = SpTanMod(1,1,2,3)
         Amat(1,9) = SpTanMod(1,1,3,3)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,3,1)
         Amat(2,4) = SpTanMod(2,1,1,2)
         Amat(2,5) = SpTanMod(2,1,2,2)
         Amat(2,6) = SpTanMod(2,1,3,2)
         Amat(2,7) = SpTanMod(2,1,1,3)
         Amat(2,8) = SpTanMod(2,1,2,3)
         Amat(2,9) = SpTanMod(2,1,3,3)
         Amat(3,1) = SpTanMod(3,1,1,1)
         Amat(3,2) = SpTanMod(3,1,2,1)
         Amat(3,3) = SpTanMod(3,1,3,1)
         Amat(3,4) = SpTanMod(3,1,1,2)
         Amat(3,5) = SpTanMod(3,1,2,2)
         Amat(3,6) = SpTanMod(3,1,3,2)
         Amat(3,7) = SpTanMod(3,1,1,3)
         Amat(3,8) = SpTanMod(3,1,2,3)
         Amat(3,9) = SpTanMod(3,1,3,3)
         Amat(4,1) = SpTanMod(1,2,1,1)
         Amat(4,2) = SpTanMod(1,2,2,1)
         Amat(4,3) = SpTanMod(1,2,3,1)
         Amat(4,4) = SpTanMod(1,2,1,2)
         Amat(4,5) = SpTanMod(1,2,2,2)
         Amat(4,6) = SpTanMod(1,2,3,2)
         Amat(4,7) = SpTanMod(1,2,1,3)
         Amat(4,8) = SpTanMod(1,2,2,3)
         Amat(4,9) = SpTanMod(1,2,3,3)
         Amat(5,1) = SpTanMod(2,2,1,1)
         Amat(5,2) = SpTanMod(2,2,2,1)
         Amat(5,3) = SpTanMod(2,2,3,1)
         Amat(5,4) = SpTanMod(2,2,1,2)
         Amat(5,5) = SpTanMod(2,2,2,2)
         Amat(5,6) = SpTanMod(2,2,3,2)
         Amat(5,7) = SpTanMod(2,2,1,3)
         Amat(5,8) = SpTanMod(2,2,2,3)
         Amat(5,9) = SpTanMod(2,2,3,3)
         Amat(6,1) = SpTanMod(3,2,1,1)
         Amat(6,2) = SpTanMod(3,2,2,1)
         Amat(6,3) = SpTanMod(3,2,3,1)
         Amat(6,4) = SpTanMod(3,2,1,2)
         Amat(6,5) = SpTanMod(3,2,2,2)
         Amat(6,6) = SpTanMod(3,2,3,2)
         Amat(6,7) = SpTanMod(3,2,1,3)
         Amat(6,8) = SpTanMod(3,2,2,3)
         Amat(6,9) = SpTanMod(3,2,3,3)
         Amat(7,1) = SpTanMod(1,3,1,1)
         Amat(7,2) = SpTanMod(1,3,2,1)
         Amat(7,3) = SpTanMod(1,3,3,1)
         Amat(7,4) = SpTanMod(1,3,1,2)
         Amat(7,5) = SpTanMod(1,3,2,2)
         Amat(7,6) = SpTanMod(1,3,3,2)
         Amat(7,7) = SpTanMod(1,3,1,3)
         Amat(7,8) = SpTanMod(1,3,2,3)
         Amat(7,9) = SpTanMod(1,3,3,3)
         Amat(8,1) = SpTanMod(2,3,1,1)
         Amat(8,2) = SpTanMod(2,3,2,1)
         Amat(8,3) = SpTanMod(2,3,3,1)
         Amat(8,4) = SpTanMod(2,3,1,2)
         Amat(8,5) = SpTanMod(2,3,2,2)
         Amat(8,6) = SpTanMod(2,3,3,2)
         Amat(8,7) = SpTanMod(2,3,1,3)
         Amat(8,8) = SpTanMod(2,3,2,3)
         Amat(8,9) = SpTanMod(2,3,3,3)
         Amat(9,1) = SpTanMod(3,3,1,1)
         Amat(9,2) = SpTanMod(3,3,2,1)
         Amat(9,3) = SpTanMod(3,3,3,1)
         Amat(9,4) = SpTanMod(3,3,1,2)
         Amat(9,5) = SpTanMod(3,3,2,2)
         Amat(9,6) = SpTanMod(3,3,3,2)
         Amat(9,7) = SpTanMod(3,3,1,3)
         Amat(9,8) = SpTanMod(3,3,2,3)
         Amat(9,9) = SpTanMod(3,3,3,3)
         !
         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         !
         ! Finally, add the contribution of this integration point
         !  to the tangent matrix.
         !
         Kuu = Kuu + detMapJC*w(intpt)*
     +             (
     +             matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +             )

      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface tension terms
      !
      ! Based on our convention, the face on which the surface tension
      !  acts is the "label", i.e.
      !  - U1,U2,U3,U4,U5,U6 refer to surface tension acting
      !    on faces 1,2,3,4,5,6, respectively,
      ! Mechanical, traction- and pressure-type boundary conditions 
      !  may be applied to the dummy mesh using the Abaqus built-in 
      !  commands *Dload or *Dsload.
      !
      if(ndload.gt.0) then
        !
        ! Loop over faces and make proper modifications to
        !  residuals and tangents as needed.
        !
        do i=1,ndload
          !
          face = jdltyp(i,1)  ! face label
          gamma = adlmag(i,1) ! surface tension
          !
          if(face.eq.1) then
            !
            ! surface tension on face 1 of the element
            !
            ! Initialize
            !
            Rsurf = zero
            Ksurf = zero
            nListSurf(1) = 1
            nListSurf(2) = 2
            nListSurf(3) = 3
            nListSurf(4) = 4
            
            
            ! Obtain current nodal coordinates on the surface
            !
            do k=1,nNodeSurf
               coordsCSurf(1+nDim*(k-1),1) = coordsC(1,nListSurf(k))
               coordsCSurf(2+nDim*(k-1),1) = coordsC(2,nListSurf(k))
               coordsCSurf(3+nDim*(k-1),1) = coordsC(3,nListSurf(k))
            end do


            ! Obtain integration point local coordinates and weights
	    !
	    if(nIntSurf.eq.4) then
              call xint2D4pt(xiSurf,wSurf,nIntPtSurf) ! 4-pt integration, nIntPtSurf=4 above
            else
              write(*,*) 'Incorrect number of surface int points'
              call xit
            endif
             
	     
	    ! Loop over surface integration points
	    !
            do intptSurf=1,nIntPtSurf
            
            
              ! Obtain shape functions and their local gradients
              !
              call calcShape2DLinear(nIntPtSurf,xiSurf,intptSurf,
     +                                 shSurf,dshxiSurf)
                
                
              ! Assemble matrices of local shape function gradients
              !
              dshxiMat1 = zero
	      dshxiMat2 = zero
	      do k=1,nNodeSurf
		dshxiMat1(1,1+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(2,2+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(3,3+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat2(1,1+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(2,2+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(3,3+nDim*(k-1)) = dshxiSurf(k,2)
              end do
                
                
              ! Calculate the covariant metric components
              !
              E = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat1,coordsCSurf)))
              G = sum((matmul(dshxiMat2,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              F = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              H = dsqrt(E*G - F*F)
                
                
              ! Calculate the derivatives of the metric components wrt
              !  the current coordinates (or equivalently the dofs)
              !
              dE = transpose(two*matmul(transpose(dshxiMat1),
     +                        matmul(dshxiMat1,coordsCSurf)))
              dG = transpose(two*matmul(transpose(dshxiMat2),
     +                        matmul(dshxiMat2,coordsCSurf)))
              dF = transpose(matmul(transpose(dshxiMat1),
     +                    matmul(dshxiMat2,coordsCSurf)) + 
     +             matmul(transpose(dshxiMat2),
     +                    matmul(dshxiMat1,coordsCSurf)))
              dH = (half/H)*(G*dE + E*dG - two*F*dF)
                
                
              ! Compute the surface residual vector
              !
              Rsurf = Rsurf - (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      )
                
                
              ! Compute the surface tangent matrix
              !
              Ksurf = Ksurf + (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(transpose(dshxiMat1),dshxiMat1) - 
     +                      F*matmul(transpose(dshxiMat2),dshxiMat1) + 
     +                      E*matmul(transpose(dshxiMat2),dshxiMat2) - 
     +                      F*matmul(transpose(dshxiMat1),dshxiMat2) +
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat1),coordsCSurf),dG) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat1),coordsCSurf),dF) + 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat2),coordsCSurf),dE) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat2),coordsCSurf),dF) - 
     +                      (one/H)*matmul(
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      ),dH)
     +                      )
     
     
            end do
            
            
            ! Modify the residual vector, loop over nodes
            !
            do k = 1,nNodeSurf
              Ru(1+nDim*(nListSurf(k)-1),1) = 
     +        Ru(1+nDim*(nListSurf(k)-1),1) + Rsurf(1+nDim*(k-1),1)
              Ru(2+nDim*(nListSurf(k)-1),1) = 
     +        Ru(2+nDim*(nListSurf(k)-1),1) + Rsurf(2+nDim*(k-1),1)
              Ru(3+nDim*(nListSurf(k)-1),1) = 
     +        Ru(3+nDim*(nListSurf(k)-1),1) + Rsurf(3+nDim*(k-1),1)
            end do
            
            
            ! Modify the tangent matrix
            !
            do j = 1,nNodeSurf
              do k = 1,nNodeSurf
                Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),1+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),2+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),3+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),1+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),2+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),3+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),1+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),2+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),3+nDim*(k-1))
              end do
            end do
            !
          elseif(face.eq.2) then
            !
            ! surface tension on face 2 of the element
            !
            ! Initialize
            !
            Rsurf = zero
            Ksurf = zero
            nListSurf(1) = 5
            nListSurf(2) = 6
            nListSurf(3) = 7
            nListSurf(4) = 8
            
            
            ! Obtain current nodal coordinates on the surface
            !
            do k=1,nNodeSurf
               coordsCSurf(1+nDim*(k-1),1) = coordsC(1,nListSurf(k))
               coordsCSurf(2+nDim*(k-1),1) = coordsC(2,nListSurf(k))
               coordsCSurf(3+nDim*(k-1),1) = coordsC(3,nListSurf(k))
            end do


            ! Obtain integration point local coordinates and weights
	    !
	    if(nIntSurf.eq.4) then
              call xint2D4pt(xiSurf,wSurf,nIntPtSurf) ! 4-pt integration, nIntPtSurf=4 above
            else
              write(*,*) 'Incorrect number of surface int points'
              call xit
            endif
             
	     
	    ! Loop over surface integration points
	    !
            do intptSurf=1,nIntPtSurf
            
            
              ! Obtain shape functions and their local gradients
              !
              call calcShape2DLinear(nIntPtSurf,xiSurf,intptSurf,
     +                                 shSurf,dshxiSurf)
                
                
              ! Assemble matrices of local shape function gradients
              !
              dshxiMat1 = zero
	      dshxiMat2 = zero
	      do k=1,nNodeSurf
		dshxiMat1(1,1+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(2,2+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(3,3+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat2(1,1+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(2,2+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(3,3+nDim*(k-1)) = dshxiSurf(k,2)
              end do
                
                
              ! Calculate the covariant metric components
              !
              E = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat1,coordsCSurf)))
              G = sum((matmul(dshxiMat2,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              F = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              H = dsqrt(E*G - F*F)
                
                
              ! Calculate the derivatives of the metric components wrt
              !  the current coordinates (or equivalently the dofs)
              !
              dE = transpose(two*matmul(transpose(dshxiMat1),
     +                        matmul(dshxiMat1,coordsCSurf)))
              dG = transpose(two*matmul(transpose(dshxiMat2),
     +                        matmul(dshxiMat2,coordsCSurf)))
              dF = transpose(matmul(transpose(dshxiMat1),
     +                    matmul(dshxiMat2,coordsCSurf)) + 
     +             matmul(transpose(dshxiMat2),
     +                    matmul(dshxiMat1,coordsCSurf)))
              dH = (half/H)*(G*dE + E*dG - two*F*dF)
                
                
              ! Compute the surface residual vector
              !
              Rsurf = Rsurf - (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      )
                
                
              ! Compute the surface tangent matrix
              !
              Ksurf = Ksurf + (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(transpose(dshxiMat1),dshxiMat1) - 
     +                      F*matmul(transpose(dshxiMat2),dshxiMat1) + 
     +                      E*matmul(transpose(dshxiMat2),dshxiMat2) - 
     +                      F*matmul(transpose(dshxiMat1),dshxiMat2) +
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat1),coordsCSurf),dG) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat1),coordsCSurf),dF) + 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat2),coordsCSurf),dE) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat2),coordsCSurf),dF) - 
     +                      (one/H)*matmul(
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      ),dH)
     +                      )
     
     
            end do
            
            
            ! Modify the residual vector, loop over nodes
            !
            do k = 1,nNodeSurf
              Ru(1+nDim*(nListSurf(k)-1),1) = 
     +        Ru(1+nDim*(nListSurf(k)-1),1) + Rsurf(1+nDim*(k-1),1)
              Ru(2+nDim*(nListSurf(k)-1),1) = 
     +        Ru(2+nDim*(nListSurf(k)-1),1) + Rsurf(2+nDim*(k-1),1)
              Ru(3+nDim*(nListSurf(k)-1),1) = 
     +        Ru(3+nDim*(nListSurf(k)-1),1) + Rsurf(3+nDim*(k-1),1)
            end do
            
            
            ! Modify the tangent matrix
            !
            do j = 1,nNodeSurf
              do k = 1,nNodeSurf
                Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),1+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),2+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),3+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),1+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),2+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),3+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),1+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),2+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),3+nDim*(k-1))
              end do
            end do
            !
          elseif(face.eq.3) then
            !
            ! surface tension on face 3 of the element
            !
            ! Initialize
            !
            Rsurf = zero
            Ksurf = zero
            nListSurf(1) = 1
            nListSurf(2) = 5
            nListSurf(3) = 6
            nListSurf(4) = 2
            
            
            ! Obtain current nodal coordinates on the surface
            !
            do k=1,nNodeSurf
               coordsCSurf(1+nDim*(k-1),1) = coordsC(1,nListSurf(k))
               coordsCSurf(2+nDim*(k-1),1) = coordsC(2,nListSurf(k))
               coordsCSurf(3+nDim*(k-1),1) = coordsC(3,nListSurf(k))
            end do


            ! Obtain integration point local coordinates and weights
	    !
	    if(nIntSurf.eq.4) then
              call xint2D4pt(xiSurf,wSurf,nIntPtSurf) ! 4-pt integration, nIntPtSurf=4 above
            else
              write(*,*) 'Incorrect number of surface int points'
              call xit
            endif
             
	     
	    ! Loop over surface integration points
	    !
            do intptSurf=1,nIntPtSurf
            
            
              ! Obtain shape functions and their local gradients
              !
              call calcShape2DLinear(nIntPtSurf,xiSurf,intptSurf,
     +                                 shSurf,dshxiSurf)
                
                
              ! Assemble matrices of local shape function gradients
              !
              dshxiMat1 = zero
	      dshxiMat2 = zero
	      do k=1,nNodeSurf
		dshxiMat1(1,1+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(2,2+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(3,3+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat2(1,1+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(2,2+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(3,3+nDim*(k-1)) = dshxiSurf(k,2)
              end do
                
                
              ! Calculate the covariant metric components
              !
              E = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat1,coordsCSurf)))
              G = sum((matmul(dshxiMat2,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              F = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              H = dsqrt(E*G - F*F)
                
                
              ! Calculate the derivatives of the metric components wrt
              !  the current coordinates (or equivalently the dofs)
              !
              dE = transpose(two*matmul(transpose(dshxiMat1),
     +                        matmul(dshxiMat1,coordsCSurf)))
              dG = transpose(two*matmul(transpose(dshxiMat2),
     +                        matmul(dshxiMat2,coordsCSurf)))
              dF = transpose(matmul(transpose(dshxiMat1),
     +                    matmul(dshxiMat2,coordsCSurf)) + 
     +             matmul(transpose(dshxiMat2),
     +                    matmul(dshxiMat1,coordsCSurf)))
              dH = (half/H)*(G*dE + E*dG - two*F*dF)
                
                
              ! Compute the surface residual vector
              !
              Rsurf = Rsurf - (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      )
                
                
              ! Compute the surface tangent matrix
              !
              Ksurf = Ksurf + (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(transpose(dshxiMat1),dshxiMat1) - 
     +                      F*matmul(transpose(dshxiMat2),dshxiMat1) + 
     +                      E*matmul(transpose(dshxiMat2),dshxiMat2) - 
     +                      F*matmul(transpose(dshxiMat1),dshxiMat2) +
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat1),coordsCSurf),dG) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat1),coordsCSurf),dF) + 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat2),coordsCSurf),dE) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat2),coordsCSurf),dF) - 
     +                      (one/H)*matmul(
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      ),dH)
     +                      )
     
     
            end do
            
            
            ! Modify the residual vector, loop over nodes
            !
            do k = 1,nNodeSurf
              Ru(1+nDim*(nListSurf(k)-1),1) = 
     +        Ru(1+nDim*(nListSurf(k)-1),1) + Rsurf(1+nDim*(k-1),1)
              Ru(2+nDim*(nListSurf(k)-1),1) = 
     +        Ru(2+nDim*(nListSurf(k)-1),1) + Rsurf(2+nDim*(k-1),1)
              Ru(3+nDim*(nListSurf(k)-1),1) = 
     +        Ru(3+nDim*(nListSurf(k)-1),1) + Rsurf(3+nDim*(k-1),1)
            end do
            
            
            ! Modify the tangent matrix
            !
            do j = 1,nNodeSurf
              do k = 1,nNodeSurf
                Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),1+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),2+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),3+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),1+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),2+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),3+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),1+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),2+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),3+nDim*(k-1))
              end do
            end do
            !
          elseif(face.eq.4) then
            !
            ! surface tension on face 4 of the element
            !
            ! Initialize
            !
            Rsurf = zero
            Ksurf = zero
            nListSurf(1) = 2
            nListSurf(2) = 3
            nListSurf(3) = 7
            nListSurf(4) = 6
            
            
            ! Obtain current nodal coordinates on the surface
            !
            do k=1,nNodeSurf
               coordsCSurf(1+nDim*(k-1),1) = coordsC(1,nListSurf(k))
               coordsCSurf(2+nDim*(k-1),1) = coordsC(2,nListSurf(k))
               coordsCSurf(3+nDim*(k-1),1) = coordsC(3,nListSurf(k))
            end do


            ! Obtain integration point local coordinates and weights
	    !
	    if(nIntSurf.eq.4) then
              call xint2D4pt(xiSurf,wSurf,nIntPtSurf) ! 4-pt integration, nIntPtSurf=4 above
            else
              write(*,*) 'Incorrect number of surface int points'
              call xit
            endif
             
	     
	    ! Loop over surface integration points
	    !
            do intptSurf=1,nIntPtSurf
            
            
              ! Obtain shape functions and their local gradients
              !
              call calcShape2DLinear(nIntPtSurf,xiSurf,intptSurf,
     +                                 shSurf,dshxiSurf)
                
                
              ! Assemble matrices of local shape function gradients
              !
              dshxiMat1 = zero
	      dshxiMat2 = zero
	      do k=1,nNodeSurf
		dshxiMat1(1,1+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(2,2+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(3,3+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat2(1,1+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(2,2+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(3,3+nDim*(k-1)) = dshxiSurf(k,2)
              end do
                
                
              ! Calculate the covariant metric components
              !
              E = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat1,coordsCSurf)))
              G = sum((matmul(dshxiMat2,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              F = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              H = dsqrt(E*G - F*F)
                
                
              ! Calculate the derivatives of the metric components wrt
              !  the current coordinates (or equivalently the dofs)
              !
              dE = transpose(two*matmul(transpose(dshxiMat1),
     +                        matmul(dshxiMat1,coordsCSurf)))
              dG = transpose(two*matmul(transpose(dshxiMat2),
     +                        matmul(dshxiMat2,coordsCSurf)))
              dF = transpose(matmul(transpose(dshxiMat1),
     +                    matmul(dshxiMat2,coordsCSurf)) + 
     +             matmul(transpose(dshxiMat2),
     +                    matmul(dshxiMat1,coordsCSurf)))
              dH = (half/H)*(G*dE + E*dG - two*F*dF)
                
                
              ! Compute the surface residual vector
              !
              Rsurf = Rsurf - (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      )
                
                
              ! Compute the surface tangent matrix
              !
              Ksurf = Ksurf + (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(transpose(dshxiMat1),dshxiMat1) - 
     +                      F*matmul(transpose(dshxiMat2),dshxiMat1) + 
     +                      E*matmul(transpose(dshxiMat2),dshxiMat2) - 
     +                      F*matmul(transpose(dshxiMat1),dshxiMat2) +
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat1),coordsCSurf),dG) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat1),coordsCSurf),dF) + 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat2),coordsCSurf),dE) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat2),coordsCSurf),dF) - 
     +                      (one/H)*matmul(
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      ),dH)
     +                      )
     
     
            end do
            
            
            ! Modify the residual vector, loop over nodes
            !
            do k = 1,nNodeSurf
              Ru(1+nDim*(nListSurf(k)-1),1) = 
     +        Ru(1+nDim*(nListSurf(k)-1),1) + Rsurf(1+nDim*(k-1),1)
              Ru(2+nDim*(nListSurf(k)-1),1) = 
     +        Ru(2+nDim*(nListSurf(k)-1),1) + Rsurf(2+nDim*(k-1),1)
              Ru(3+nDim*(nListSurf(k)-1),1) = 
     +        Ru(3+nDim*(nListSurf(k)-1),1) + Rsurf(3+nDim*(k-1),1)
            end do
            
            
            ! Modify the tangent matrix
            !
            do j = 1,nNodeSurf
              do k = 1,nNodeSurf
                Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),1+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),2+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),3+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),1+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),2+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),3+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),1+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),2+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),3+nDim*(k-1))
              end do
            end do
            !
          elseif(face.eq.5) then
            !
            ! surface tension on face 5 of the element
            !
            ! Initialize
            !
            Rsurf = zero
            Ksurf = zero
            nListSurf(1) = 4
            nListSurf(2) = 8
            nListSurf(3) = 7
            nListSurf(4) = 3
            
            
            ! Obtain current nodal coordinates on the surface
            !
            do k=1,nNodeSurf
               coordsCSurf(1+nDim*(k-1),1) = coordsC(1,nListSurf(k))
               coordsCSurf(2+nDim*(k-1),1) = coordsC(2,nListSurf(k))
               coordsCSurf(3+nDim*(k-1),1) = coordsC(3,nListSurf(k))
            end do


            ! Obtain integration point local coordinates and weights
	    !
	    if(nIntSurf.eq.4) then
              call xint2D4pt(xiSurf,wSurf,nIntPtSurf) ! 4-pt integration, nIntPtSurf=4 above
            else
              write(*,*) 'Incorrect number of surface int points'
              call xit
            endif
             
	     
	    ! Loop over surface integration points
	    !
            do intptSurf=1,nIntPtSurf
            
            
              ! Obtain shape functions and their local gradients
              !
              call calcShape2DLinear(nIntPtSurf,xiSurf,intptSurf,
     +                                 shSurf,dshxiSurf)
                
                
              ! Assemble matrices of local shape function gradients
              !
              dshxiMat1 = zero
	      dshxiMat2 = zero
	      do k=1,nNodeSurf
		dshxiMat1(1,1+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(2,2+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(3,3+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat2(1,1+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(2,2+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(3,3+nDim*(k-1)) = dshxiSurf(k,2)
              end do
                
                
              ! Calculate the covariant metric components
              !
              E = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat1,coordsCSurf)))
              G = sum((matmul(dshxiMat2,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              F = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              H = dsqrt(E*G - F*F)
                
                
              ! Calculate the derivatives of the metric components wrt
              !  the current coordinates (or equivalently the dofs)
              !
              dE = transpose(two*matmul(transpose(dshxiMat1),
     +                        matmul(dshxiMat1,coordsCSurf)))
              dG = transpose(two*matmul(transpose(dshxiMat2),
     +                        matmul(dshxiMat2,coordsCSurf)))
              dF = transpose(matmul(transpose(dshxiMat1),
     +                    matmul(dshxiMat2,coordsCSurf)) + 
     +             matmul(transpose(dshxiMat2),
     +                    matmul(dshxiMat1,coordsCSurf)))
              dH = (half/H)*(G*dE + E*dG - two*F*dF)
                
                
              ! Compute the surface residual vector
              !
              Rsurf = Rsurf - (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      )
                
                
              ! Compute the surface tangent matrix
              !
              Ksurf = Ksurf + (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(transpose(dshxiMat1),dshxiMat1) - 
     +                      F*matmul(transpose(dshxiMat2),dshxiMat1) + 
     +                      E*matmul(transpose(dshxiMat2),dshxiMat2) - 
     +                      F*matmul(transpose(dshxiMat1),dshxiMat2) +
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat1),coordsCSurf),dG) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat1),coordsCSurf),dF) + 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat2),coordsCSurf),dE) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat2),coordsCSurf),dF) - 
     +                      (one/H)*matmul(
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      ),dH)
     +                      )
     
     
            end do
            
            
            ! Modify the residual vector, loop over nodes
            !
            do k = 1,nNodeSurf
              Ru(1+nDim*(nListSurf(k)-1),1) = 
     +        Ru(1+nDim*(nListSurf(k)-1),1) + Rsurf(1+nDim*(k-1),1)
              Ru(2+nDim*(nListSurf(k)-1),1) = 
     +        Ru(2+nDim*(nListSurf(k)-1),1) + Rsurf(2+nDim*(k-1),1)
              Ru(3+nDim*(nListSurf(k)-1),1) = 
     +        Ru(3+nDim*(nListSurf(k)-1),1) + Rsurf(3+nDim*(k-1),1)
            end do
            
            
            ! Modify the tangent matrix
            !
            do j = 1,nNodeSurf
              do k = 1,nNodeSurf
                Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),1+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),2+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),3+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),1+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),2+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),3+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),1+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),2+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),3+nDim*(k-1))
              end do
            end do
            !
          elseif(face.eq.6) then
            !
            ! surface tension on face 6 of the element
            !
            ! Initialize
            !
            Rsurf = zero
            Ksurf = zero
            nListSurf(1) = 1
            nListSurf(2) = 4
            nListSurf(3) = 8
            nListSurf(4) = 5
            
            
            ! Obtain current nodal coordinates on the surface
            !
            do k=1,nNodeSurf
               coordsCSurf(1+nDim*(k-1),1) = coordsC(1,nListSurf(k))
               coordsCSurf(2+nDim*(k-1),1) = coordsC(2,nListSurf(k))
               coordsCSurf(3+nDim*(k-1),1) = coordsC(3,nListSurf(k))
            end do


            ! Obtain integration point local coordinates and weights
	    !
	    if(nIntSurf.eq.4) then
              call xint2D4pt(xiSurf,wSurf,nIntPtSurf) ! 4-pt integration, nIntPtSurf=4 above
            else
              write(*,*) 'Incorrect number of surface int points'
              call xit
            endif
             
	     
	    ! Loop over surface integration points
	    !
            do intptSurf=1,nIntPtSurf
            
            
              ! Obtain shape functions and their local gradients
              !
              call calcShape2DLinear(nIntPtSurf,xiSurf,intptSurf,
     +                                 shSurf,dshxiSurf)
                
                
              ! Assemble matrices of local shape function gradients
              !
              dshxiMat1 = zero
	      dshxiMat2 = zero
	      do k=1,nNodeSurf
		dshxiMat1(1,1+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(2,2+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat1(3,3+nDim*(k-1)) = dshxiSurf(k,1)
		dshxiMat2(1,1+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(2,2+nDim*(k-1)) = dshxiSurf(k,2)
		dshxiMat2(3,3+nDim*(k-1)) = dshxiSurf(k,2)
              end do
                
                
              ! Calculate the covariant metric components
              !
              E = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat1,coordsCSurf)))
              G = sum((matmul(dshxiMat2,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              F = sum((matmul(dshxiMat1,coordsCSurf))*
     +              (matmul(dshxiMat2,coordsCSurf)))
              H = dsqrt(E*G - F*F)
                
                
              ! Calculate the derivatives of the metric components wrt
              !  the current coordinates (or equivalently the dofs)
              !
              dE = transpose(two*matmul(transpose(dshxiMat1),
     +                        matmul(dshxiMat1,coordsCSurf)))
              dG = transpose(two*matmul(transpose(dshxiMat2),
     +                        matmul(dshxiMat2,coordsCSurf)))
              dF = transpose(matmul(transpose(dshxiMat1),
     +                    matmul(dshxiMat2,coordsCSurf)) + 
     +             matmul(transpose(dshxiMat2),
     +                    matmul(dshxiMat1,coordsCSurf)))
              dH = (half/H)*(G*dE + E*dG - two*F*dF)
                
                
              ! Compute the surface residual vector
              !
              Rsurf = Rsurf - (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      )
                
                
              ! Compute the surface tangent matrix
              !
              Ksurf = Ksurf + (gamma/H)*wSurf(intptSurf)*
     +                      (
     +                      G*matmul(transpose(dshxiMat1),dshxiMat1) - 
     +                      F*matmul(transpose(dshxiMat2),dshxiMat1) + 
     +                      E*matmul(transpose(dshxiMat2),dshxiMat2) - 
     +                      F*matmul(transpose(dshxiMat1),dshxiMat2) +
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat1),coordsCSurf),dG) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat1),coordsCSurf),dF) + 
     +                      matmul(matmul(matmul(transpose(dshxiMat2),
     +                            dshxiMat2),coordsCSurf),dE) - 
     +                      matmul(matmul(matmul(transpose(dshxiMat1),
     +                            dshxiMat2),coordsCSurf),dF) - 
     +                      (one/H)*matmul(
     +                      (
     +                      G*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat1),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat1),coordsCSurf) + 
     +                      E*matmul(matmul(transpose(dshxiMat2),
     +                                      dshxiMat2),coordsCSurf) - 
     +                      F*matmul(matmul(transpose(dshxiMat1),
     +                                      dshxiMat2),coordsCSurf)
     +                      ),dH)
     +                      )
     
     
            end do
            
            
            ! Modify the residual vector, loop over nodes
            !
            do k = 1,nNodeSurf
              Ru(1+nDim*(nListSurf(k)-1),1) = 
     +        Ru(1+nDim*(nListSurf(k)-1),1) + Rsurf(1+nDim*(k-1),1)
              Ru(2+nDim*(nListSurf(k)-1),1) = 
     +        Ru(2+nDim*(nListSurf(k)-1),1) + Rsurf(2+nDim*(k-1),1)
              Ru(3+nDim*(nListSurf(k)-1),1) = 
     +        Ru(3+nDim*(nListSurf(k)-1),1) + Rsurf(3+nDim*(k-1),1)
            end do
            
            
            ! Modify the tangent matrix
            !
            do j = 1,nNodeSurf
              do k = 1,nNodeSurf
                Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),1+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),2+nDim*(k-1))
                Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(1+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(1+nDim*(j-1),3+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),1+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),2+nDim*(k-1))
                Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(2+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(2+nDim*(j-1),3+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),1+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),1+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),2+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),2+nDim*(k-1))
                Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) = 
     +          Kuu(3+nDim*(nListSurf(j)-1),3+nDim*(nListSurf(k)-1)) + 
     +            Ksurf(3+nDim*(j-1),3+nDim*(k-1))
              end do
            end do
            !
          else
            write(*,*) 'Incorrect dload type',face
            call xit
          end if
          !
        end do
        !
      end if
      !
      ! End loop over surface tension terms
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Return to Abaqus the RHS vector and the Stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return to Abaqus the right hand side vector.
      !
      do i=1,nNode
         !
         rhs(nDim*(i-1)+1,1) = Ru(nDim*(i-1)+1,1)
         rhs(nDim*(i-1)+2,1) = Ru(nDim*(i-1)+2,1)
         rhs(nDim*(i-1)+3,1) = Ru(nDim*(i-1)+3,1)
         !
      enddo
      !
      ! Return to Abaqus the tangent matrix.
      !
      amatrx = Kuu
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine uel

!************************************************************************
!     Material subroutines
!************************************************************************

      subroutine NeoHookean(props,nprops,dtime,F_tau,Fv_t,
     +        T_tau,Fv_tau,SpTanMod,stat)

      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fv_t(3,3),T_tau(3,3),
     +  Fv_tau(3,3),SpTanMod(3,3,3,3),Iden(3,3),G0,Kbulk,eta,
     +  Gneq,detF,Finv(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),
     +  I1_tau,I1bar,TR_tau(3,3),dTRdF(3,3,3,3),Fv_t_inv(3,3),
     +  det_Fv_t,Fe_tr(3,3),Re_tr(3,3),Ue_tr(3,3),Ee_tr(3,3),
     +  trEe_tr,Ee0_tr(3,3),Me_tr(3,3),tauBar_tr,Nv_tau(3,3),
     +  nuv_tau,tauBar_tau,Dv_tau(3,3),Dv_eig(3),Dv_vec(3,3),
     +  expdtDv(3,3),tmp,Me_tau(3,3),Gshear_tilde,Lambda_tilde,
     +  fac,c3
      !
      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0)
 

      ! Identity matrix
      !
      call onem(Iden)


      ! Obtain material properties
      !
      G0     = props(1) ! Ground-state shear modulus
      Kbulk  = props(2) ! Bulk modulus
      Gneq   = props(4) ! Nonequilibrium (Maxwell element) stiffness
      eta    = props(5) ! Maxwell element viscosity
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Neo-Hookean spring
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT = transpose(F_tau)
      FinvT = transpose(Finv)


      ! Compute the right Cauchy-Green tensor
      !  and the first stretch invariant
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      I1_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 

      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*I1_tau
 

      ! Compute the 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*G0*(F_tau-third*I1_tau*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
 

      ! Compute the Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))


      ! Calculate the material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*G0*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*I1_tau*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*I1_tau*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the spatial tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) +
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Maxwell element
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the trial elastic deformation gradient
      !
      call matInv3D(Fv_t,Fv_t_inv,det_Fv_t,stat)
      Fe_tr = matmul(F_tau,Fv_t_inv)


      ! Compute the trial kinematics
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr,stat)
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3)
      Ee0_tr = Ee_tr - (one/three)*trEe_tr*Iden


      ! Compute the trial Mandel stress, which is deviatoric
      !
      Me_tr = two*Gneq*Ee0_tr


      ! Compute the trial equiv. tensile stress
      !
      tauBar_tr = dsqrt(one/two)*dsqrt(sum(Me_tr*Me_tr))


      ! Compute the direction of viscous flow
      !
      if(tauBar_tr.le.zero) then
         Nv_tau = zero
      else
         Nv_tau = dsqrt(one/two)*(Me_tr/tauBar_tr)
      endif


      ! Compute the equivalent shear viscous strain rate
      !
      nuv_tau = tauBar_tr/(eta + Gneq*dtime)


      ! Compute the equivalent shear stress
      !
      tauBar_tau = tauBar_tr - Gneq*dtime*nuv_tau


      ! Compute the viscous stretching
      !
      Dv_tau = dsqrt(one/two)*nuv_tau*Nv_tau


      ! Compute the viscous deformation gradient at the
      !  end of the increment using the exponential map
      !
      if(nuv_tau.le.zero) then
         Fv_tau = Fv_t
      else
         call spectral(dtime*Dv_tau,Dv_eig,Dv_vec,stat)
         expdtDv = zero
         expdtDv(1,1) = dexp(Dv_eig(1))
         expdtDv(2,2) = dexp(Dv_eig(2))
         expdtDv(3,3) = dexp(Dv_eig(3))
         expdtDv = matmul(matmul(Dv_vec,expdtDv),transpose(Dv_vec))
         Fv_tau = matmul(expdtDv,Fv_t)
      endif


      ! Check to make sure that det(Fv_tau)>0
      !
      call mdet(Fv_tau,tmp)
      if(tmp.le.zero) then
         write(*,*) 'det(Fv_tau).le.zero in INTEG'
        stat = 0
        return
      endif


      ! Compute the Mandel stres at the end of the increment
      !
      Me_tau = Me_tr - two*Gneq*dtime*Dv_tau


      ! Compute the contribution to the Cauchy stress due to
      !  viscoelasticity
      !
      T_tau = T_tau + 
     +        matmul(Re_tr,matmul(Me_tau,transpose(Re_tr)))/detF


      ! Compute contribution to the spatial tangent due to
      !  viscoelasticity
      !
      Nv_tau = matmul(Re_tr,matmul(Nv_tau,transpose(Re_tr)))
      if (tauBar_tr.gt.zero) then
	   Gshear_tilde = (tauBar_tau/tauBar_tr)*Gneq
      else
	   Gshear_tilde = Gneq
      end if
      fac = (one + Gneq*dtime/eta)**(-one)
      Lambda_tilde =  - Gshear_tilde*two/three
      if (tauBar_tr.gt.zero) then
         c3 = -two*Gneq*((tauBar_tau/tauBar_tr) - fac)
      else
         c3 = -two*Gneq*(one - fac)
      end if
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) 
     +          +Gshear_tilde*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k))
     +          +Lambda_tilde*Iden(i,j)*Iden(k,l)
     +          +c3*Nv_tau(i,j)*Nv_tau(k,l)
               enddo
            enddo
         enddo
      enddo


      return
      end subroutine NeoHookean 
      
!****************************************************************************

      subroutine Gent(props,nprops,dtime,F_tau,Fv_t,
     +        T_tau,Fv_tau,SpTanMod,stat)
      !
      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fv_t(3,3),T_tau(3,3),
     +  Fv_tau(3,3),SpTanMod(3,3,3,3),Iden(3,3),G0,Kbulk,eta,Im,
     +  Gneq,detF,Finv(3,3),FT(3,3),FinvT(3,3),C_tau(3,3),Cinv(3,3),
     +  detC,trC,I1bar,fac,GShearGent,dGdF(3,3),TR_tau(3,3),
     +  dTRdF(3,3,3,3),Fv_t_inv(3,3),det_Fv_t,Fe_tr(3,3),Re_tr(3,3),
     +  Ue_tr(3,3),Ee_tr(3,3),trEe_tr,Ee0_tr(3,3),Me_tr(3,3),
     +  tauBar_tr,Nv_tau(3,3),nuv_tau,tauBar_tau,Dv_tau(3,3),Dv_eig(3),
     +  Dv_vec(3,3),expdtDv(3,3),tmp,Me_tau(3,3),Gshear_tilde,
     +  Lambda_tilde,fac2,c3
      !
      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
      !
      G0     = props(1) ! Ground-state shear modulus
      Kbulk  = props(2) ! Bulk modulus
      Im     = props(3) ! Limiting chain extensibility parameter
      Gneq   = props(4) ! Nonequilibrium (Maxwell element) stiffness
      eta    = props(5) ! Maxwell element viscosity
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Gent spring
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the determinant, the inverse, the transpose, 
      !  and the inverse transpose of the deformation gradient
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FT = transpose(F_tau)
      FinvT = transpose(Finv)


      ! Compute the right Cauchy-Green tensor, its inverse, 
      !  its determinant, and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      call matInv3D(C_tau,Cinv,detC,stat)
      trC = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
 

      ! Compute the trace of the distortional right Cauchy-Green tensor
      !
      I1bar = (detF**(-two/three))*trC
 

      ! Compute the ``I-3'' factor appearing in the Gent model
      !
      fac = (I1bar - three)/Im
      if(fac.gt.0.95d0) fac = 0.95d0
      fac = one/(one - fac)
 

      ! Compute the ``shear'' modulus. Note: this is not really the shear
      !  modulus, but it will help when computing the material tangent later
      !
      GShearGent = G0*fac
 

      ! Compute the derivative of the ``shear modulus'' with respect
      !  to the deformation gradient for use in the material tangent
      !
      dGdF = two*(G0/Im)*(detF**(-two/three))*
     +     fac*fac*(F_tau - third*trC*FinvT)
 

      ! Compute the 1st Piola stress
      !
      TR_tau = (detF**(-two/three))*GShearGent*(F_tau-third*trC*FinvT)
     +     + Kbulk*detF*(detF - one)*FinvT
 

      ! Compute the Cauchy stress
      !
      T_tau = (one/detF)*matmul(TR_tau,transpose(F_tau))


      ! Compute the material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3                  
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + (detF**(-two/three))*dGdF(k,l)*
     +                 (
     +                 F_tau(i,j) - third*trC*Finv(j,i)
     +                 )
     +                 + (detF**(-two/three))*GshearGent*
     +                 (
     +                 (-two/three)*F_tau(i,j)*Finv(l,k)
     +                 + (two/9.d0)*trC*Finv(j,i)*Finv(l,k)
     +                 + Iden(i,k)*Iden(j,l)
     +                 + third*trC*Finv(l,i)*Finv(j,k)
     +                 - (two/three)*Finv(j,i)*F_tau(k,l)
     +                 )
     +                 + detF*Kbulk*
     +                 (
     +                 (detF-one)*Finv(j,i)*Finv(l,k)
     +                 + detF*Finv(j,i)*Finv(l,k)
     +                 - (detF-one)*Finv(l,i)*Finv(j,k)
     +                 )
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the spatial tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute the contribution due to the Maxwell element
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Compute the trial elastic deformation gradient
      !
      call matInv3D(Fv_t,Fv_t_inv,det_Fv_t,stat)
      Fe_tr = matmul(F_tau,Fv_t_inv)


      ! Compute the trial kinematics
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr,stat)
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3)
      Ee0_tr = Ee_tr - (one/three)*trEe_tr*Iden


      ! Compute the trial Mandel stress, which is deviatoric
      !
      Me_tr = two*Gneq*Ee0_tr


      ! Compute the trial equiv. tensile stress
      !
      tauBar_tr = dsqrt(one/two)*dsqrt(sum(Me_tr*Me_tr))


      ! Compute the direction of viscous flow
      !
      if(tauBar_tr.le.zero) then
         Nv_tau = zero
      else
         Nv_tau = dsqrt(one/two)*(Me_tr/tauBar_tr)
      endif


      ! Compute the equivalent shear viscous strain rate
      !
      nuv_tau = tauBar_tr/(eta + Gneq*dtime)


      ! Compute the equivalent shear stress
      !
      tauBar_tau = tauBar_tr - Gneq*dtime*nuv_tau


      ! Compute the viscous stretching
      !
      Dv_tau = dsqrt(one/two)*nuv_tau*Nv_tau


      ! Compute the viscous deformation gradient at the
      !  end of the increment using the exponential map
      !
      if(nuv_tau.le.zero) then
         Fv_tau = Fv_t
      else
         call spectral(dtime*Dv_tau,Dv_eig,Dv_vec,stat)
         expdtDv = zero
         expdtDv(1,1) = dexp(Dv_eig(1))
         expdtDv(2,2) = dexp(Dv_eig(2))
         expdtDv(3,3) = dexp(Dv_eig(3))
         expdtDv = matmul(matmul(Dv_vec,expdtDv),transpose(Dv_vec))
         Fv_tau = matmul(expdtDv,Fv_t)
      endif


      ! Check to make sure that det(Fv_tau)>0
      !
      call mdet(Fv_tau,tmp)
      if(tmp.le.zero) then
         write(*,*) 'det(Fv_tau).le.zero in INTEG'
        stat = 0
        return
      endif


      ! Compute the Mandel stres at the end of the increment
      !
      Me_tau = Me_tr - two*Gneq*dtime*Dv_tau


      ! Compute the contribution to the Cauchy stress due to
      !  viscoelasticity
      !
      T_tau = T_tau + 
     +        matmul(Re_tr,matmul(Me_tau,transpose(Re_tr)))/detF


      ! Compute contribution to the spatial tangent due to
      !  viscoelasticity
      !
      Nv_tau = matmul(Re_tr,matmul(Nv_tau,transpose(Re_tr)))
      if (tauBar_tr.gt.zero) then
	   Gshear_tilde = (tauBar_tau/tauBar_tr)*Gneq
      else
	   Gshear_tilde = Gneq
      end if
      fac2 = (one + Gneq*dtime/eta)**(-one)
      Lambda_tilde =  - Gshear_tilde*two/three
      if (tauBar_tr.gt.zero) then
         c3 = -two*Gneq*((tauBar_tau/tauBar_tr) - fac2)
      else
         c3 = -two*Gneq*(one - fac2)
      end if
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) 
     +          +Gshear_tilde*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k))
     +          +Lambda_tilde*Iden(i,j)*Iden(k,l)
     +          +c3*Nv_tau(i,j)*Nv_tau(k,l)
               enddo
            enddo
         enddo
      enddo
      

      return
      end subroutine Gent

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi_1,xi_2,xi_3 coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(8,3),w(8)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

!************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i):      shape function of node i at the intpt.
      ! dshxi(i,j): derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt,i,j
      !
      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3),xi,eta,zeta
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      
      
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      
      
      return
      end subroutine calcShape3DLinear

!*************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi_1-xi_2-xi_3 domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode),mapJ(3,3),
     +  mapJ_inv(3,3),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      return
      end subroutine mapShape3D
      
!****************************************************************************
!     Surface element subroutines
!****************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi_1,xi_2 coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2),w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt
      
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i):      shape function of node i at the intpt.
      ! dshxi(i,j): derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear
      
!****************************************************************************
!     The next subroutine calculates various kinematical quantities
!      associated with the deformation gradient
!****************************************************************************

      subroutine skinem(F,R,U,E,istat)
      !
      ! This subroutine performs the right polar decomposition
      !  F = RU of the deformation gradient F into a rotation
      !  R and the right stretch tensor U.  The logarithmic 
      !  strain E = ln(U) is also returned.
      !
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	U(3,3):       the right stretch tensor; output
      !	Uinv(3,3):    the inverse of U
      !	C(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Ueigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of U
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
      implicit none
      !
      integer istat
      !
      real*8 F(3,3),C(3,3),omega(3),Ueigval(3),eigvec(3,3),
     +  U(3,3),E(3,3),Uinv(3,3),R(3,3),detF
     

      !	Store the identity matrix in R, U, and Uinv
      !
      call onem(R)
      call onem(U)
      call onem(Uinv)
      

      ! Store the zero matrix in E
      !
      E = 0.d0
      

      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the 
      !  time increment.
      !
      call mdet(F,detF)
      if (detF.le.0.d0) then
        write(*,'(/5X,A/)') '--problem in kinematics-- the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
      

      ! Calculate the right Cauchy-Green tensor C
      !
      C = matmul(transpose(F),F)
      
 
      ! Calculate the eigenvalues and eigenvectors of C
      !
      call spectral(C,omega,eigvec,istat)
      

      ! Calculate the principal values of U and E
      !
      Ueigval(1) = dsqrt(omega(1))
      Ueigval(2) = dsqrt(omega(2))
      Ueigval(3) = dsqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
      !
      E(1,1) = dlog(Ueigval(1))
      E(2,2) = dlog(Ueigval(2))
      E(3,3) = dlog(Ueigval(3))
      

      ! Calculate the complete tensors U and E
      !
      U = matmul(matmul(eigvec,U),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
      

      ! Calculate Uinv
      !
      call matInv3D(U,Uinv,detF,istat)
      

      ! calculate R
      !
      R = matmul(F,Uinv)
      

      return
      end subroutine skinem

!****************************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!****************************************************************************

      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)


      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      call eigsrt(D,V,3,np)
	

      return
      end subroutine spectral
	
!****************************************************************************

      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
     
      
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
	B(ip) = A(ip,ip)
	D(ip) = B(ip)
	Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the 
              !  off-diagonal element is small.
              !
	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +            .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
	          T =A(ip,iq)/H
	        else
	          theta = 0.5d0*H/A(ip,iq)
	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
	          if (theta.lt.0.d0) T = -T
	        end if
	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
	        S = T*C
	        tau = S/(1.d0 + C)
	        H = T*A(ip,iq)
	        Z(ip) = Z(ip) - H
	        Z(iq) = Z(iq) + H
	        D(ip) = D(ip) - H
	        D(iq) = D(iq) + H
	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
		!		
	        do j=1,ip-1
	          G = A(j,ip)
	          H = A(j,iq)
	          A(j,ip) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations P < J < Q
                !
	        do j=ip+1,iq-1
	          G = A(ip,j)
	          H = A(j,iq)
	          A(ip,j) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations Q < J <= N
                !
	        do j=iq+1,n
                  G = A(ip,j)
	          H = A(iq,j)
	          A(ip,j) = G - S*(H + G*tau)
	          A(iq,j) = H + S*(G - H*tau)
	        end do
	        do j = 1,n
	          G = V(j,ip)
	          H = V(j,iq)
	          V(j,ip) = G - S*(H + G*tau)
	          V(j,iq) = H + S*(G - H*tau)
	        end do
	        nrot = nrot + 1
              end if
	    end do
	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
	  do ip=1,n
	    B(ip) = B(ip) + Z(ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  end do
	end do


      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
      

      return
      end subroutine jacobi
	
!****************************************************************************

      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  jacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      

      do i=1,n-1
	k = i
	P = D(i)
	do j=i+1,n
	  if (D(j).ge.P) then
	    k = j
	    P = D(j)
	  end if
	end do
	if (k.ne.i) then
	  D(k) = D(i)
	  D(i) = P
	  do j=1,n
	    P = V(j,i)
	    V(j,i) = V(j,k)
	    V(j,k) = P
	  end do
  	end if
      end do
      

      return
      end subroutine eigsrt

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: SUBROUTINE matInv3:'
        write(*,*) 'WARNING: DET of MAT=',DET_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************