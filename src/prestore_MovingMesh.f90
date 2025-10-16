MODULE PRESTORE_MOVINGMESH
USE DECLARATION
USE DERIVATIVES
USE LIBRARY
USE BASIS
USE TRANSFORM
USE TRANSFORM_MovingMesh
USE LOCAL
USE LAPCK
USE DG_FUNCTIONS
USE PRESTORE

IMPLICIT NONE


 CONTAINS
 
 
SUBROUTINE RE_PRESTORE_1(N, node_position_index)
	!> @brief
	!> This subroutine calls other subroutines to store again the pseudoinverse reconstruction least square matrices after the mesh movement
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N, node_position_index
	INTEGER::KMAXE,i,m,K
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
	REAL,ALLOCATABLE,DIMENSION(:,:)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
	REAL,ALLOCATABLE,DIMENSION(:,:)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
	REAL,ALLOCATABLE,DIMENSION(:,:)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
	REAL,ALLOCATABLE,DIMENSION(:,:)::ILOX_VOLUME    !CELL VOLUME
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::ILOX_PERIODICFLAG
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:)::ILON_NODCOUNT  !NUMBER OF NODES
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::ILON_X           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::ILON_Y           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::ILON_Z           !COORDINATES OF EACH NODE IN X

	character*3 N_string
	character*50 filename_invmat_stencilt,filename_lscqmat

	KMAXE=XMPIELRANK(N)
	M=TYPESTEN

	ALLOCATE (ILOX_IHEXG(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_IHEXL(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_IHEXB(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_IHEXN(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_PERIODICFLAG(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_ISHAPE(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_XXC(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_VOLUME(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOX_YYC(M,NUMNEIGHBOURS*iextend))
	IF (DIMENSIONA.EQ.3)THEN
		ALLOCATE (ILOX_ZZC(M,NUMNEIGHBOURS*iextend))
		ILOX_ZZC=0.d0
	END IF
	ILOX_IHEXG=0
	ILOX_IHEXL=0
	ILOX_IHEXB=0
	ILOX_ISHAPE=0
	ILOX_IHEXN=0
	ILOX_VOLUME=0.d0
	ILOX_XXC=0.d0
	ILOX_YYC=0.d0

	IF (DIMENSIONA.EQ.3)THEN
		K=8
	ELSE
		K=4
	END IF

	ALLOCATE (ILON_NODCOUNT(M,NUMNEIGHBOURS*iextend,K))
	ALLOCATE (ILON_X(M,NUMNEIGHBOURS*iextend,K))
	ALLOCATE (ILON_Y(M,NUMNEIGHBOURS*iextend,K))
	IF (DIMENSIONA.EQ.3)THEN
		ALLOCATE (ILON_Z(M,NUMNEIGHBOURS*iextend,K))
		ILON_Z=0.D0
	END IF
	ILON_NODCOUNT=0
	ILON_X=0.D0
	ILON_Y=0.D0

	KMAXE=XMPIELRANK(N)

	if (dimensiona.eq.3)then
		! !$OMP DO
		! DO I=1,KMAXE
		! 	CALL LOCALISE_STENCIL_movement(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
		! 	CALL LOCALISE_STEN2_movement(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
		! 	! call CHECKGRADS(N,I)
		! 	CALL FIND_ROT_ANGLES_movement(N,I)

		! 	CALL PRESTORE_RECONSTRUCTION3_movement(N,i,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)

		! 	! if ((dg.eq.1).OR.(ADDA_type.EQ.2))then
		! 	! 	CALL PRESTORE_DG1(i)
		! 	! end if
		! END DO
		! !$OMP END DO
	else
		!$OMP DO
		DO I=1,KMAXE
			CALL LOCALISE_STENCIL_MovingMesh_2D(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z, node_position_index)
			CALL LOCALISE_STENCIL_STEP2_MovingMesh_2D(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z, node_position_index)
			! call CHECKGRADS2d(N,I)
			CALL FIND_ROT_ANGLES_MovingMesh_2D(N,I,node_position_index)
			CALL PRESTORE_RECONSTRUCTION_MovingMesh_2D(N,i,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z,node_position_index)
			! if ((dg.eq.1).OR.(ADDA_type.EQ.2))then
			! 	CALL PRESTORE_DG1(i)
			! end if
		END DO
		!$OMP END DO
	end if

	deallocate(ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y)

	IF (DIMENSIONA.EQ.3)THEN
		deallocate(ILOX_ZZC,ILON_Z)
	end if

END SUBROUTINE RE_PRESTORE_1






SUBROUTINE PRESTORE_RECONSTRUCTION_MovingMesh_3D(N,iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	!> @brief
	!> This subroutine prestores the pseudoinverse reconstruction least square matrices
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,iconsi
	INTEGER::I,J,K,llco,ll,ii,igf,IGF2,IFD2,idum,idum2,iq,jq,lq,IHGT,IHGJ,iqp,iqp2,NND,k0,g0,lcou,lcc,iqqq,ICOND1,ICOND2,N_NODE
	integer::ideg,ideg2,imax,imax2,ivgt,jxx,ixx,lxx1,kxx,icompwrt,number_of_dog,ELTYPE,inumo,inumo2,inum,facex,iconsidered
	real::ssss,gggg,UPTEMP,LOTEMP,X_STENCIL,Y_STENCIL,Z_STENCIL,DIST_STEN,DIST_STEN2
	real::ax,ay,az,ANGLE1,angle2,nx,ny,nz,nnx,nnz,nny,x1,y1,z1
	REAL,DIMENSION(1:8,1:dimensiona)::NODES_LIST
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	real,DIMENSION(1:dims,1:dims)::AINVJT
	real,dimension(1:dimensiona)::cords
	real,allocatable,dimension(:)::INTBS,BASEFACEVAL,BASEFACGVAL,PERMUTATION,PERMUTATIONG,xder,yder,zder
	real,allocatable,dimension(:,:)::stencil,invmat,WLSQR,LSQM,LSCQM,QFF,RFF,QTFF,INVRFF,VELLSQMAT
	real,allocatable,dimension(:,:,:)::stencilS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X


	! allocate(LSCQM(1:idegfree,1:idegfree),QFF(1:idegfree,1:idegfree),RFF(1:idegfree,1:idegfree),QTFF(1:idegfree,1:idegfree),INVRFF(1:idegfree,1:idegfree))
	! allocate(VELLSQMAT(1:idegfree-1,1:idegfree-1))
	! allocate(LSQM(1:IMAXDEGFREE,1:IDEGFREE-1))	
	! allocate(INTBS(1:idegfree),BASEFACEVAL(1:idegfree),BASEFACGVAL(1:idegfree),PERMUTATION(1:idegfree),PERMUTATIONG(1:idegfree),xder(1:idegfree),yder(1:idegfree),zder(1:idegfree))
	! allocate(invmat(1:idegfree,1:idegfree))
	! allocate(WLSQR(1:20,1:numneighbours-1))
	! allocate(stencil(1:numneighbours-1,1:idegfree))
	! allocate(stencilS(1:20,1:numneighbours-1,1:idegfree))


	! i=iconsi

	! IDUM=0;
	! if (ielem(n,i)%interior.eq.1)then
	! 	DO j=1,IELEM(N,I)%IFCA
	! 		if (ielem(n,i)%ibounds(J).gt.0)then
	! 			if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
	! 				IDUM=1
	! 			end if
	! 		END IF
	! 	END DO
	! end if      

	! iconsidered=i

	! !firstly compute the integrals of the basis functions
	! INTBS=zero;JXX=1;IXX=i;LXX1=1;number_of_dog=ielem(n,i)%idegfree;kxx=ielem(n,i)%iorder;ELTYPE=ielem(n,i)%ishape
	! icompwrt=0

	! INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)

	! INTEG_BASIS(I)%VALUE(1:ielem(n,i)%IDEGFREE)=INTBS(1:ielem(n,i)%IDEGFREE)
			
	! !the indicator matrix for the smaller polynomial
	! IF (IWENO.EQ.1)THEN
	! 	CALL INDICATORMATRIX(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! END IF
		
	! !the indicator matrix for the lower order polynomial
	! if (ees.eq.5)then
	! 	INTBS=zero;JXX=1;IXX=i;LXX1=1;number_of_dog=idegfree2;kxx=IORDER2;ELTYPE=ielem(n,i)%ishape
	! 	icompwrt=1
	! 	INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 	INTEG_BASIS(I)%VALUEc(1:number_of_dog)=INTBS(1:number_of_dog)
	! 	IF (IWENO.EQ.1)THEN
	! 		CALL INDICATORMATRIX2(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 		icompwrt=0
	! 	END IF
	! end if

	! LLCO=IELEM(N,I)%ADMIS;IMAX=IELEM(N,I)%inumneighbours-1;INUM=IELEM(N,I)%inumneighbours;IDEG=IELEM(N,I)%iDEGFREE
	! INUMO=ielem(n,i)%iorder;imax2=numneighbours2-1;inum2=numneighbours2;ideg2=iDEGFREE2;inumo2=iorder2

	! !this section computes the ratio of volumes inside each stencil
	! DIST_STEN=ZERO

	! DIST_STEN=-tolbig; DIST_STEN2=tolbig
	! DO LL=1,ielem(n,i)%iNUMNEIGHBOURS
	! 	if (ILOX_VOLUME(1,LL).lt.DIST_STEN2)then
	! 		DIST_STEN2=ILOX_VOLUME(1,LL)
	! 	end if
	! 	if (ILOX_VOLUME(1,LL).gt.DIST_STEN)then
	! 		DIST_STEN=ILOX_VOLUME(1,LL)
	! 	end if
	! end do
			
	! IELEM(N,I)%STENCIL_DIST=MAX((DIST_STEN/DIST_STEN2),(DIST_STEN2/DIST_STEN))
				
	! !now we start with looping all the admissible stencils
	! DO LL=1,LLCO	!for all stencils
	! 	if((ees.ne.5).or.(ll.eq.1))then
	! 		IMAX=IELEM(N,I)%inumneighbours-1;INUM=IELEM(N,I)%inumneighbours;IDEG=IELEM(N,I)%iDEGFREE;INUMO=ielem(n,i)%iorder
	! 		icompwrt=0;number_of_dog=IDEG
	! 	else
	! 		imax=numneighbours2-1;inum=numneighbours2;ideg=iDEGFREE2;inumo=IORDER2;number_of_dog=IDEG
	! 		icompwrt=1
	! 	end if

	! 	DO K=1,imax	!for all neighbours
	! 		ixx=i
	! 		kxx=INUMO

	! 		X_STENCIL=(ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1))**2
	! 		Y_STENCIL=(ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1))**2
	! 		Z_STENCIL=(ILOX_ZZC(ll,k+1)-ILOX_ZZC(ll,1))**2

	! 		DIST_STEN2=SQRT(X_STENCIL+Y_STENCIL+Z_STENCIL)

	! 		IF (WEIGHT_LSQR.EQ.1)THEN
	! 			WLSQR(ll,K)=1.0D0/((DIST_STEN2))
	! 		ELSE
	! 			WLSQR(ll,K)=1.0D0
	! 		END IF

	! 		if (fastest.eq.1)then	!this is when transformation is not active (Rarely used)
	! 			x1 = ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1)
	! 			y1 = ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1)
	! 			z1 = ILOX_ZZC(ll,k+1)-ILOX_ZZC(ll,1)

	! 			IF (GREENGO.EQ.0)THEN	!for the least squares green gauss gradients
	! 				if (idum.eq.1)then
	! 					if((ees.ne.5).or.(ll.eq.1))then
	! 						icompwrt=0
	! 						ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,ielem(n,i)%iorder,I,ielem(n,i)%idegfree,icompwrt)
	! 						ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
	! 					else
	! 						icompwrt=1
	! 						ILOCAL_RECON3(I)%STENCILSc(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,I,IDEG,icompwrt)
	! 						ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
	! 					end if
	! 				else
	! 					if((ees.ne.5).or.(ll.eq.1))then
	! 						icompwrt=0
	! 						STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)
	! 					else
	! 						icompwrt=1
	! 						STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)
	! 					end if
	! 				end if
	! 			ELSE
	! 				if((ees.ne.5).or.(ll.eq.1))then
	! 					icompwrt=0
	! 					STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)
	! 				else
	! 					icompwrt=1
	! 					STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)

	! 				end if
	! 			END IF

	! 		else		!with coordinate transformation

	! 			IXX=i
	! 			jxx=k+1
	! 			lxx1=ll

	! 			ELTYPE=ILOX_ishape(ll,k+1)

	! 			IF (GREENGO.EQ.0)THEN
	! 				if (idum.eq.1)then	!smaller memory footprint for non boundary elements
	! 					if((ees.ne.5).or.(ll.eq.1))then
	! 						icompwrt=0

	! 						ILOCAL_RECON3(I)%STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 						ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
	! 					else
	! 						icompwrt=1
	! 						ILOCAL_RECON3(I)%STENCILSc(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 						ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
	! 					end if
	! 				else
	! 					if((ees.ne.5).or.(ll.eq.1))then
	! 						icompwrt=0
	! 						STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 					else
	! 						icompwrt=1
	! 						STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 					end if
	! 				end if
	! 			ELSE
	! 				if((ees.ne.5).or.(ll.eq.1))then
	! 					icompwrt=0
	! 					STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 				else
	! 					icompwrt=1
	! 					STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	! 				end if
	! 			END IF
	! 		end if
	! 	END DO	!close loop for all neighbours

	! 	IF ((GREENGO.EQ.0))THEN
	! 		if (idum.eq.1)then
	! 			invmat=zero;LSCQM=zero
	! 			DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
	! 				if((ees.ne.5).or.(ll.eq.1))then
	! 					LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
	! 						+((ILOCAL_RECON3(I)%STENCILS(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)))
	! 				else
	! 					LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
	! 						+((ILOCAL_RECON3(I)%STENCILSc(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILSc(LL,LQ,IQ)))
	! 				end if
	! 			END DO;END DO;END DO
	! 		else
	! 			invmat=zero;LSCQM=zero
	! 			DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
	! 				LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
	! 					+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
	! 			END DO;END DO;END DO
	! 		end if
	! 	ELSE
	! 		invmat=zero;LSCQM=zero
	! 		DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
	! 			LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
	! 				+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
	! 		END DO;END DO;END DO
	! 	END IF

	! 	! svd_diagonal = ZERO
	! 	! matrix_copy_u = ZERO
	! 	! matrix_copy_u(1:ideg,1:ideg) = LSCQM(1:ideg,1:ideg)
	! 	! ! WRITE(*,*) N, "copied"
	! 	! CALL svdcmp(matrix_copy_u, ideg, ideg, idegfree, 100, svd_diagonal, matrix_copy_v)
	! 	! ! WRITE(*,*) N, "SVD done"
	! 	! max_eigv = 0.0
	! 	! min_eigv = 1000000000000.0
	! 	! do iter = 1,imax
	! 	! 	if (abs(svd_diagonal(iter)) >= 0.000000001) then
	! 	! 		if (abs(svd_diagonal(iter)) > max_eigv) then
	! 	! 			max_eigv = abs(svd_diagonal(iter))
	! 	! 		end if
	! 	! 		if (abs(svd_diagonal(iter)) < min_eigv) then
	! 	! 			min_eigv = abs(svd_diagonal(iter))
	! 	! 		end if
	! 	! 	end if
	! 	! end do
	! 	! condition_number = max_eigv / min_eigv

	! 	! write(8000+N, '(G0, X, G0, X, G0, X, G0, X, G0, X, G0)') IELEM(N,i )%XXC, IELEM(N,i)%YYC, IELEM(N,i)%ZZC, condition_number, max_eigv, min_eigv

	! 	QFF(:,:)=zero; RFF(:,:)=zero; QTFF(:,:)=zero; RFF(:,:)=zero;  INVRFF(:,:)=zero
	! 	CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IDEG)
	! 	CALL TRANSPOSEMATRIX(QFF,QTFF,IDEG)
	! 	IVGT=IDEG+1
	! 	CALL INVERT(RFF,INVRFF,IVGT)
	! 	invmat(1:IDEg,1:ideg)=MATMUL(INVRFF(1:ideg,1:IDEG),QTFF(1:IDEG,1:IDEG))

	! 	IF (GREENGO.EQ.0)THEN
	! 		if (idum.eq.1)then
	! 			if((ees.ne.5).or.(ll.eq.1))then
	! 				stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:ideg)
	! 			else
	! 				stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILSc(LL,1:imax,1:ideg)
	! 			end if
	! 		else
	! 			stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)
	! 		end if
	! 	ELSE
	! 		stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)
	! 	END IF

	! 	if((ees.ne.5).or.(ll.eq.1))then
	! 		! call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
	! 		! stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL),IDEG)

	! 		ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))

	! 		do iq=1,imax
	! 			ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)&
	! 				*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
	! 		end do
	! 	else
	! 		! call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
	! 		! stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stenciltc(1:IDEG,1:IMAX,LL),IDEG)

	! 		ILOCAL_RECON3(i)%invmat_stenciltC(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))

	! 		do iq=1,imax
	! 			ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)&
	! 				*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
	! 		end do
	! 	end if

	! 	! svd_diagonal = ZERO
	! 	! matrix_copy_u = ZERO
	! 	! matrix_copy_u(1:ideg,1:imax) = ILOCAL_RECON3(I)%invmat_stencilt(1:ideg,1:imax,1)
	! 	! ! WRITE(*,*) N, "copied"
	! 	! CALL svdcmp(matrix_copy_u, ideg, imax, idegfree, 100, svd_diagonal, matrix_copy_v)
	! 	! ! WRITE(*,*) N, "SVD done"
	! 	! max_eigv = 0.0
	! 	! min_eigv = 1000000000000.0
	! 	! do iter = 1,imax
	! 	! 	if (abs(svd_diagonal(iter)) >= 0.000000001) then
	! 	! 		if (abs(svd_diagonal(iter)) > max_eigv) then
	! 	! 			max_eigv = abs(svd_diagonal(iter))
	! 	! 		end if
	! 	! 		if (abs(svd_diagonal(iter)) < min_eigv) then
	! 	! 			min_eigv = abs(svd_diagonal(iter))
	! 	! 		end if
	! 	! 	end if
	! 	! end do
	! 	! condition_number = max_eigv / min_eigv

	! 	! write(7000+N, '(G0, X, G0, X, G0, X, G0, X, G0, X, G0)') IELEM(N,i )%XXC, IELEM(N,i)%YYC, IELEM(N,i)%ZZC, condition_number, max_eigv, min_eigv

	! 	if (ielem(n,i)%ggs.ne.1)then	!ggs
	! 		IF (ITESTCASE.EQ.4)THEN	!TEST
	! 			IF (LL.EQ.1)THEN		!stencils
	! 				if (idum.eq.1)then		!for wall only	
	! 					DO IHGT=1,3; DO IHGJ=1,3
	! 						AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
	! 					END DO; END DO

	! 					idum2=0
		
	! 					BASEFACEVAL=zero
	! 					BASEFACGVAL=zero
	! 					PERMUTATION=zero
	! 					PERMUTATIONg=zero

	! 					DO j=1,IELEM(N,I)%IFCA		!for all faces
	! 						if (ielem(n,i)%ibounds(J).gt.0)then		!for bounded only
	! 							if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then		!!for bounded only 2                         
	! 								facex=J;iconsidered=i
	! 								CALL coordinates_face_inner(N,Iconsidered,facex,vext,NODES_LIST)

	! 								if (ielem(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
	! 									N_NODE=4
	! 								else
	! 									N_NODE=3
	! 								end if

	! 								CORDS(1:3)=zero
	! 								CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
								
	! 								AY=cords(2)
	! 								AX=cords(1)
	! 								AZ=cords(3)
				
	! 								VEXT(1,1)=AX;VEXT(1,2)=AY;VEXT(1,3)=AZ
	! 								VEXT(1,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(1,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
								
	! 								AX=VEXT(1,1);AY=VEXT(1,2);AZ=VEXT(1,3)
				
	! 								ANGLE1=IELEM(N,I)%FACEANGLEX(j)
	! 								ANGLE2=IELEM(N,I)%FACEANGLEY(j)
									
	! 								NX=(COS(ANGLE1)*SIN(ANGLE2))
	! 								NY=(SIN(ANGLE1)*SIN(ANGLE2))
	! 								NZ=(COS(ANGLE2))
	! 								NNX=(NX*AINVJT(1,1))+(NY*AINVJT(2,1))+(NZ*AINVJT(3,1))
	! 								NNY=(NX*AINVJT(1,2))+(NY*AINVJT(2,2))+(NZ*AINVJT(3,2))
	! 								NNZ=(NX*AINVJT(1,3))+(NY*AINVJT(2,3))+(NZ*AINVJT(3,3))
				
	! 								DO IQ=1, IDEG
	! 									IF (POLY.EQ.1) THEN
	! 										XDER(IQ)=DFX(AX,AY,AZ,IQ,I);  YDER(IQ)=DFY(AX,AY,AZ,IQ,I);  ZDER(IQ)=DFZ(AX,AY,AZ,IQ,I)
	! 									END IF
	! 									IF (POLY.EQ.2) THEN
	! 										XDER(IQ)=DLX(AX,AY,AZ,IQ,I);  YDER(IQ)=DLY(AX,AY,AZ,IQ,I);  ZDER(IQ)=DLZ(AX,AY,AZ,IQ,I)
	! 									END IF
	! 									IF (POLY.EQ.4) THEN
	! 										XDER(IQ)=TL3DX(AX,AY,AZ,IQ,I);  YDER(IQ)=TL3DY(AX,AY,AZ,IQ,I);  ZDER(IQ)=TL3DZ(AX,AY,AZ,IQ,I)
	! 									END IF
	! 								END DO

	! 								icompwrt=0
	! 								BASEFACEVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,ielem(n,i)%Iorder,I,ielem(n,i)%IDEGFREE,icompwrt)

	! 								if (thermal.eq.0) then
	! 									BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=((NNX*XDER(1:ielem(n,i)%IDEGFREE))+(NNY*YDER(1:ielem(n,i)%IDEGFREE))+(NNZ*ZDER(1:ielem(n,i)%IDEGFREE)))
	! 								else
	! 									icompwrt=0
	! 									BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,ielem(n,i)%Iorder,I,ielem(n,i)%IDEGFREE,icompwrt)
	! 								end if

	! 								DO IQ=1,IDEG
	! 									ILOCAL_RECON3(I)%WALLCOEFF(IQ)=BASEFACEVAL(IQ)
	! 									ILOCAL_RECON3(I)%WALLCOEFG(IQ)=BASEFACGVAL(IQ)
	! 									PERMUTATION(IQ)=IQ
	! 									PERMUTATIONG(IQ)=IQ
	! 								END DO
	! 							end if!for bounded only 2        
	! 						end if!for bounded only
	! 					end do!for all faces
		
	! 					GGGG=-TOLBIG
	! 					G0=0
	! 					DO IQ=1,IDEG
	! 						IF (ABS(BASEFACGVAL(IQ)) >GGGG)THEN
	! 							GGGG=ABS(BASEFACGVAL(IQ))
	! 							G0=IQ
	! 						END IF
	! 					END DO

	! 					SSSS=-TOLBIG
	! 					K0=0
	! 					DO IQ=1,IDEG
	! 						IF (ABS(BASEFACEVAL(IQ)) >SSSS)THEN
	! 							SSSS=ABS(BASEFACEVAL(IQ))
	! 							K0=IQ
	! 						END IF
	! 					END DO

	! 					! now swap basis functions and thus coefficients
	! 					PERMUTATION(1)=K0; PERMUTATION(K0)=1; SSSS=BASEFACEVAL(1)
	! 					BASEFACEVAL(1)=BASEFACEVAL(K0); BASEFACEVAL(K0)=SSSS
	! 					PERMUTATIONG(1)=G0; PERMUTATION(G0)=1; GGGG=BASEFACGVAL(1)
	! 					BASEFACGVAL(1)=BASEFACGVAL(G0); BASEFACGVAL(G0)=GGGG
	! 					ILOCAL_RECON3(I)%K0=K0
	! 					ILOCAL_RECON3(I)%G0=G0
		
	! 					LSQM = ZERO
	! 					DO LQ=1,IMAX
	! 						LCOU=0
	! 						DO IQ=1,IDEG
	! 							IF (IQ.EQ.G0) CYCLE
	! 							LCOU=LCOU+1
	! 							LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
	! 								-ILOCAL_RECON3(I)%STENCILS(LL,LQ,G0)*ILOCAL_RECON3(I)%WALLCOEFG(IQ)/ILOCAL_RECON3(I)%WALLCOEFG(G0)
	! 						END DO
	! 					END DO
	! 					ILOCAL_RECON3(I)%TEMPSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)

	! 					VELLSQMAT=ZERO
	! 					DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
	! 						!now store the least square matrix
	! 						VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
	! 					END DO;	END DO;	END DO

	! 					LSCQM=ZERO
	! 					LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)
	! 					Qff=ZERO; Rff=ZERO; QTff=ZERO; INVRff=ZERO
	! 					IVGT=IDEG
		
	! 					CALL QRDECOMPOSITION(LSCQM,Qff,Rff,IVGT-1)
	! 					CALL TRANSPOSEMATRIX(Qff,QTff,IVGT-1)
						
	! 					CALL INVERT(Rff,INVRff,IVGT)

	! 					! final inverted R^(-1)*Q^(-1)
	! 					ILOCAL_RECON3(I)%TEMPSQMAT(1:IDEG-1,1:IDEG-1) =&
	! 					MATMUL(INVRff(1:IDEG-1,1:IDEG-1),QTff(1:IDEG-1,1:IDEG-1))

	! 					! definy the matrix defining the least-square reconstruction in this case for velocity
	! 					LSQM = ZERO
	! 					DO LQ=1,IMAX
	! 						LCOU=0
	! 						DO IQ=1,IDEG
	! 							IF (IQ.EQ.K0) CYCLE
	! 							LCOU=LCOU+1
	! 							LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
	! 								-ILOCAL_RECON3(I)%STENCILS(LL,LQ,K0)*ILOCAL_RECON3(I)%WALLCOEFF(IQ)/ILOCAL_RECON3(I)%WALLCOEFF(K0)
	! 						END DO
	! 					END DO
	! 					ILOCAL_RECON3(I)%VELLSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)

	! 					VELLSQMAT=ZERO
	! 					DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
	! 						!now store the least square matrix
	! 						VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
	! 					END DO;	END DO;	END DO

	! 					LSCQM=ZERO
	! 					LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)

	! 					Qff=ZERO; Rff=ZERO; QTff=ZERO; INVRff=ZERO
	! 					IVGT=IDEG
	! 					CALL QRDECOMPOSITION(LSCQM,Qff,Rff,IVGT-1)
	! 					CALL TRANSPOSEMATRIX(Qff,QTff,IVGT-1)

	! 					CALL INVERT(Rff,INVRff,IVGT)
	! 					! final inverted R^(-1)*Q^(-1)
	! 					ILOCAL_RECON3(I)%VELINVLSQMAT(1:IDEG-1,1:IDEG-1) = &
	! 					MATMUL(INVRff(1:IDEG-1,1:IDEG-1),QTff(1:IDEG-1,1:IDEG-1))
		
	! 				end if!for wall only
	! 			end if!stencils
	! 		end if!for test
	! 	end if!ggs
	! end do !for all stencils

	! ! deallocate(matrix_copy_u,matrix_copy_v)
	! ! deallocate(svd_diagonal)
		
	! deallocate(LSQM,QFF,RFF,QTFF,INVRFF)	
	! deallocate(VELLSQMAT)
	! deallocate(LSCQM)	
	! deallocate(INTBS,BASEFACEVAL,BASEFACGVAL,PERMUTATION,PERMUTATIONG,xder,yder,zder)	
	! deallocate(WLSQR)	
	! deallocate(stencil)	
	! deallocate(invmat)	
	! deallocate(STENCILS)	

END SUBROUTINE PRESTORE_RECONSTRUCTION_MovingMesh_3D





SUBROUTINE PRESTORE_RECONSTRUCTION_MovingMesh_2D(N,iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z, node_position_index)
	!> @brief
	!> This subroutine prestores the pseudoinverse reconstruction least square matrices in 2d
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,iconsi
	INTEGER::I,J,K,llco,ll,ii,igf,IGF2,IFD2,idum,idum2,iq,jq,lq,IHGT,IHGJ,iqp,iqp2,NND,k0,g0,lcou,lcc,iqqq,ICOND1,ICOND2,N_NODE
	integer::ideg,ideg2,imax,imax2,ivgt,jxx,ixx,lxx1,kxx,icompwrt,number_of_dog,ELTYPE,inumo,inumo2,iconsidered,facex,inum,ai,aj
	real::ssss,gggg,UPTEMP,LOTEMP,X_STENCIL,Y_STENCIL,Z_STENCIL,DIST_STEN,DIST_STEN2
	real::ax,ay,az,ANGLE1,angle2,nx,ny,nz,nnx,nnz,nny,x1,y1,z1,maxai,minai
	REAL,DIMENSION(1:8,1:dimensiona)::NODES_LIST
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	real,DIMENSION(1:dims,1:dims)::AINVJT
	real,dimension(1:dimensiona)::cords
	real,allocatable,dimension(:)::INTBS,BASEFACEVAL,BASEFACGVAL,PERMUTATION,PERMUTATIONG,xder,yder,zder
	real,allocatable,dimension(:,:)::stencil,invmat,WLSQR,LSQM,LSCQM,QFF,RFF,QTFF,INVRFF,VELLSQMAT
	real,allocatable,dimension(:,:,:)::stencilS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN Y
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN Z
	INTEGER,INTENT(IN)::node_position_index

	allocate(LSCQM(1:idegfree,1:idegfree),QFF(1:idegfree,1:idegfree),RFF(1:idegfree,1:idegfree),QTFF(1:idegfree,1:idegfree),INVRFF(1:idegfree,1:idegfree))
	allocate(VELLSQMAT(1:idegfree-1,1:idegfree-1))
	allocate(LSQM(1:IMAXDEGFREE,1:IDEGFREE-1))
	allocate(INTBS(1:idegfree),BASEFACEVAL(1:idegfree),BASEFACGVAL(1:idegfree),PERMUTATION(1:idegfree),PERMUTATIONG(1:idegfree),xder(1:idegfree),yder(1:idegfree),zder(1:idegfree))
	allocate(invmat(1:idegfree,1:idegfree))
	allocate(WLSQR(1:20,1:numneighbours-1))
	allocate(stencil(1:numneighbours-1,1:idegfree))
	allocate(stencilS(1:20,1:numneighbours-1,1:idegfree))

	i=iconsi

	IDUM=0;
	if (ielem(n,i)%interior.eq.1)then
		DO j=1,IELEM(N,I)%IFCA
			if (ielem(n,i)%ibounds(J).gt.0)then
				if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
					IDUM=1
				end if
			END IF
		END DO
	end if
			
	INTBS=zero;JXX=1;IXX=i;LXX1=1;number_of_dog=ielem(n,i)%idegfree;kxx=ielem(n,i)%iorder;ELTYPE=ielem(n,i)%ishape
	Icompwrt=0

	INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	INTEG_BASIS(I)%VALUE(1:ielem(n,i)%IDEGFREE)=INTBS(1:ielem(n,i)%IDEGFREE)
									
	!the indicator matrix for the large polynomial
	IF (IWENO.EQ.1)THEN
		CALL INDICATORMATRIX(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	END IF

	if (ees.eq.5)then
		INTBS=zero;JXX=1;IXX=i;LXX1=1
		number_of_dog=idegfree2;kxx=IORDER2;ELTYPE=ielem(n,i)%ishape;

		icompwrt=0

		INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
		
		INTEG_BASIS(I)%VALUEc(1:number_of_dog)=INTBS(1:number_of_dog)

		!the indicator matrix for the smaller polynomial
		IF (IWENO.EQ.1)THEN
			CALL INDICATORMATRIX2(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
		END IF
	end if
			
	LLCO=IELEM(N,I)%ADMIS
					
	IMAX=IELEM(N,I)%inumneighbours-1;INUM=IELEM(N,I)%inumneighbours;IDEG=IELEM(N,I)%iDEGFREE;INUMO=ielem(n,i)%iorder
	imax2=numneighbours2-1;inum2=numneighbours2;ideg2=iDEGFREE2;inumo2=iorder2
	DIST_STEN=ZERO
	!now we start with looping all the admissible stencils
	DO LL=1,LLCO	!ADMIS

		if((ees.ne.5).OR.(ll.eq.1))then
			IMAX=IELEM(N,I)%inumneighbours-1;INUM=IELEM(N,I)%inumneighbours;IDEG=IELEM(N,I)%iDEGFREE;INUMO=ielem(n,i)%iorder
			icompwrt=0;number_of_dog=IDEG
		else
			imax=numneighbours2-1;inum=numneighbours2;ideg=iDEGFREE2;inumo=IORDER2;number_of_dog=IDEG
			icompwrt=1
		end if

		DO K=1,imax !for all neighbours
			ixx=i;kxx=INUMO

			IF (WEIGHT_LSQR.EQ.1)THEN
				WLSQR(LL,K)=1.0D0/((SQRT(((ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1))**2)+((ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1))**2))))
			ELSE
				WLSQR(LL,K)=1.0D0
			END IF
			X_STENCIL=(ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1))**2
			Y_STENCIL=(ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1))**2

			DIST_STEN2=SQRT(X_STENCIL+Y_STENCIL)

			DIST_STEN=MAX(DIST_STEN,DIST_STEN2)

			IF (WEIGHT_LSQR.EQ.1)THEN
				WLSQR(ll,K)=1.0D0/SQRT(X_STENCIL+Y_STENCIL)
			ELSE
				WLSQR(ll,K)=1.0D0
			END IF

			IELEM(N,I)%STENCIL_DIST=DIST_STEN/(ILOX_VOLUME(1,1)**(1/2))

			if (fastest.eq.1)then
				x1 = ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1)
				y1 = ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1)

				if((ees.ne.5).OR.(ll.eq.1))then
					icompwrt=0
					ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*basis_rec2d(N,x1,y1,ielem(n,i)%iorder,IXX,ielem(n,i)%idegfree,icompwrt)
					ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
				else
					icompwrt=1
					ILOCAL_RECON3(I)%STENCILSc(LL,K,1:ideg)=WLSQR(ll,K)*basis_rec2d(N,x1,y1,inumo,IXX,ideg,icompwrt)
					ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
				end if
			else
				IXX=i;jxx=k+1;lxx1=ll
				ELTYPE=ILOX_ishape(ll,k+1)

				IF (GREENGO.EQ.0)THEN

					if (idum.eq.1)then
						if((ees.ne.5).OR.(ll.eq.1))then
							icompwrt=0
							ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
							ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
						else
							icompwrt=1
							ILOCAL_RECON3(I)%STENCILSc(LL,K,1:ideg)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
							ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
							icompwrt=0
						end if
					else
						if((ees.ne.5).OR.(ll.eq.1))then
							icompwrt=0
						else
							icompwrt=1
						end if
						STENCILS(LL,K,1:ideg)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
						icompwrt=0
					end if

				else
					if((ees.ne.5).OR.(ll.eq.1))then
						icompwrt=0
					else
						icompwrt=1
					end if
					STENCILS(LL,K,1:ideg)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
					icompwrt=0
				end if
			end if		!fastest
		END DO	!IMAXEDEGFREE

		IF ((GREENGO.EQ.0))THEN
			if (idum.eq.1)then
				invmat=zero
				LSCQM(:,:)=zero
				DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
					if((ees.ne.5).OR.(ll.eq.1))then
						LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
							+((ILOCAL_RECON3(I)%STENCILS(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)))
					else
						LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
							+((ILOCAL_RECON3(I)%STENCILSc(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILSc(LL,LQ,IQ)))
					end if
				END DO;END DO;END DO  
			else
				invmat=zero
				LSCQM(:,:)=zero
				DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
					LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
						+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
				END DO;END DO;END DO 
			end if
		ELSE
			invmat=zero
			LSCQM(:,:)=zero
			DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
				LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
					+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
			END DO;END DO;END DO 
		END IF 

		QFF(:,:)=zero; RFF(:,:)=zero; QTFF(:,:)=zero; RFF(:,:)=zero;  INVRFF(:,:)=zero
		CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IDEG)
		CALL TRANSPOSEMATRIX(QFF,QTFF,IDEG)
		IVGT=IDEG+1
		CALL INVERT(RFF,INVRFF,IVGT)
		invmat(1:IDEg,1:ideg)=MATMUL(INVRFF(1:ideg,1:IDEG),QTFF(1:IDEG,1:IDEG))

		IF (GREENGO.EQ.0)THEN
			if (idum.eq.1)then
				if((ees.ne.5).OR.(ll.eq.1))then    
					stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:ideg)
				else
					stencil(1:imax,1:ideg)=ILOCAL_RECON3(I)%STENCILSc(LL,1:imax,1:ideg)
				end if
			else
				stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)
			end if
		ELSE
			stencil(1:imax,1:ideg)=STENCILS(LL,1:imax,1:ideg)
		END IF

		if((ees.ne.5).OR.(ll.eq.1))then
			! call gemm(                                               &
			!           invmat,                                        &
			!           stencil,                                       &
			!           ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),      &
			!           'N',                                           & ! transposition flag for invmat
			!           'T'                                            & ! transposition flag for stencil
			!           )
							
			! call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
			! stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL),IDEG)
		
			ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))

			do iq=1,imax
				ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)&
					*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do
		else
			! call gemm(                                                      &
			!           invmat(1:IDEG,1:IDEG),                                &
			!           stencil(1:imax,1:ideg),                               &
			!           ILOCAL_RECON3(I)%invmat_stenciltc(1:ideg,1:imax,LL),  &
			!           'N',                                                  & ! transposition flag for invmat
			!           'T'                                                   & ! transposition flag for stencil
			!           )
						
			! call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
			! stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stenciltc(1:IDEG,1:IMAX,LL),IDEG)

			ILOCAL_RECON3(i)%invmat_stenciltC(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))

			do iq=1,imax
				ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)&
					*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do
		end if

		if (initcond.eq.0)then
			maxai=zero
			minai=tolbig
				
			do ai=1,imax
				do aj=1,ideg
					if (abs(ILOCAL_RECON3(I)%invmat_stencilt(aj,ai,LL)).ne.zero)then
						maxai=max(maxai,abs(ILOCAL_RECON3(I)%invmat_stencilt(aj,ai,LL)))
						minai=min(minai,abs(ILOCAL_RECON3(I)%invmat_stencilt(aj,ai,LL)))
					end if
				end do
			end do
			ilocal_recon3(i)%cond(ll)=maxai/minai        
		end if

		if (ielem(n,i)%ggs.ne.1)then	!ggs
			IF (ITESTCASE.EQ.4)THEN	!TEST
				IF (LL.EQ.1)THEN		!stencils
					if (idum.eq.1)then		!for wall only

						DO IHGT=1,2; DO IHGJ=1,2
							AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
						END DO; END DO

						idum2=0
						BASEFACEVAL=zero
						BASEFACGVAL=zero
						PERMUTATION=zero
						PERMUTATIONg=zero
			
						DO j=1,IELEM(N,I)%IFCA		!for all faces
							if (ielem(n,i)%ibounds(J).gt.0)then		!for bounded only
								if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then		!!for bounded only 2                         
									facex=J;iconsidered=i
									CALL coordinates_face_inner_MovingMesh_2D(N,Iconsidered,facex,vext,NODES_LIST, node_position_index)
									N_NODE=2
									CORDS(1:2)=zero
									CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
									
									AY=cords(2)
									AX=cords(1)

									VEXT(1,1)=AX
									VEXT(1,2)=AY
									VEXT(1,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2),VEXT(1,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
								
									AX=VEXT(1,1);AY=VEXT(1,2)
					
									ANGLE1=IELEM(N,I)%FACEANGLEX(j)
									ANGLE2=IELEM(N,I)%FACEANGLEY(j)
									
									NX=ANGLE1
									NY=ANGLE2
									
									NNX=(NX*AINVJT(1,1))+(NY*AINVJT(2,1))
									NNY=(NX*AINVJT(1,2))+(NY*AINVJT(2,2))
					
									IF (POLY.EQ.4)THEN
										DO IQ=1, IDEG
											XDER(IQ)=TL2dX(AX,AY,IQ,I);  YDER(IQ)=TL2dY(AX,AY,IQ,I);
										END DO
									ELSE
										DO IQ=1, IDEG
											XDER(IQ)=DF2dX(AX,AY,IQ,I);  YDER(IQ)=DF2dY(AX,AY,IQ,I);
										END DO
									END IF

									ICONSIDERED=I
									Icompwrt=0
									BASEFACEVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC2D(N,AX,AY,ielem(n,i)%Iorder,Iconsidered,ielem(n,i)%IDEGFREE,ICOMPWRT)

									if (thermal.eq.0)then
										Icompwrt=0
										BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=((NNX*XDER(1:ielem(n,i)%IDEGFREE))+(NNY*YDER(1:ielem(n,i)%IDEGFREE)))
									ELSE
										Icompwrt=0
										BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC2D(N,AX,AY,ielem(n,i)%Iorder,Iconsidered,ielem(n,i)%IDEGFREE,ICOMPWRT)
									end if

									DO IQ=1,IDEG
										ILOCAL_RECON3(I)%WALLCOEFF(IQ)=BASEFACEVAL(IQ)
										ILOCAL_RECON3(I)%WALLCOEFG(IQ)=BASEFACGVAL(IQ)
										PERMUTATION(IQ)=IQ
										PERMUTATIONG(IQ)=IQ
									END DO
								end if!for bounded only 2        
							end if!for bounded only
						end do!for all faces
			
						GGGG=-TOLBIG
						G0=0
						DO IQ=1,IDEG
							IF (ABS(BASEFACGVAL(IQ)) >GGGG)THEN
								GGGG=ABS(BASEFACGVAL(IQ))
								G0=IQ
							END IF
						END DO

						SSSS=-TOLBIG
						K0=0
						DO IQ=1,IDEG
							IF (ABS(BASEFACEVAL(IQ)) >SSSS)THEN
								SSSS=ABS(BASEFACEVAL(IQ))
								K0=IQ
							END IF
						END DO

						! now swap basis functions and thus coefficients
						PERMUTATION(1)=K0; PERMUTATION(K0)=1; SSSS=BASEFACEVAL(1)
						BASEFACEVAL(1)=BASEFACEVAL(K0); BASEFACEVAL(K0)=SSSS

						PERMUTATIONG(1)=G0; PERMUTATION(G0)=1; GGGG=BASEFACGVAL(1)
						BASEFACGVAL(1)=BASEFACGVAL(G0); BASEFACGVAL(G0)=GGGG
						! ILOCAL_RECON3(II)%K0=K0
						! ILOCAL_RECON3(II)%G0=G0
						ILOCAL_RECON3(I)%K0=K0
						ILOCAL_RECON3(I)%G0=G0

						LSQM = ZERO
						DO LQ=1,IMAX
							LCOU=0
							DO IQ=1,IDEG
							IF (IQ.EQ.G0) CYCLE
							LCOU=LCOU+1
							LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
								-ILOCAL_RECON3(I)%STENCILS(LL,LQ,G0)*ILOCAL_RECON3(I)%WALLCOEFG(IQ)/ILOCAL_RECON3(I)%WALLCOEFG(G0)
							END DO
						END DO
						ILOCAL_RECON3(I)%TEMPSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)

						VELLSQMAT=ZERO
						DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
							!now store the least square matrix
							VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
						END DO;	END DO;	END DO
						LSCQM=ZERO
						LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)

						QFF=ZERO; RFF=ZERO; QTFF=ZERO; INVRFF=ZERO
						IVGT=IDEG
						CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IVGT-1)
						CALL TRANSPOSEMATRIX(QFF,QTFF,IVGT-1)
						CALL INVERT(RFF,INVRFF,IVGT)
						! final inverted R^(-1)*Q^(-1)
						ILOCAL_RECON3(I)%TEMPSQMAT(1:IDEG-1,1:IDEG-1) =&
								MATMUL(INVRFF(1:IDEG-1,1:IDEG-1),QTFF(1:IDEG-1,1:IDEG-1))

						! definy the matrix defining the least-square reconstruction in this case for velocity
						LSQM = ZERO
						DO LQ=1,IMAX
							LCOU=0
							DO IQ=1,IDEG
								IF (IQ.EQ.K0) CYCLE
								LCOU=LCOU+1
								LSQM(LQ,LCOU)=ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)&
										-ILOCAL_RECON3(I)%STENCILS(LL,LQ,K0)*ILOCAL_RECON3(I)%WALLCOEFF(IQ)/ILOCAL_RECON3(I)%WALLCOEFF(K0)
							END DO
						END DO
						ILOCAL_RECON3(I)%VELLSQ(1:IMAX,1:IDEG-1)=LSQM(1:IMAX,1:IDEG-1)
						VELLSQMAT=ZERO

						DO IQ=1,IDEG-1; DO JQ=1,IDEG-1;	DO LCC=1,IMAX
							!now store the least square matrix
							VELLSQMAT(JQ,IQ)= VELLSQMAT(JQ,IQ)+(LSQM(LCC,JQ)*LSQM(LCC,IQ))
						END DO;	END DO;	END DO

						LSCQM=ZERO
						LSCQM(1:IDEG-1,1:IDEG-1)=VELLSQMAT(1:IDEG-1,1:IDEG-1)

						QFF=ZERO; RFF=ZERO; QTFF=ZERO; INVRFF=ZERO
						IVGT=IDEG
						CALL QRDECOMPOSITION(LSCQM,QFF,RFF,IVGT-1)
						CALL TRANSPOSEMATRIX(QFF,QTFF,IVGT-1)
						CALL INVERT(RFF,INVRFF,IVGT)

						! final inverted R^(-1)*Q^(-1)
						ILOCAL_RECON3(I)%VELINVLSQMAT(1:IDEG-1,1:IDEG-1) = &
								MATMUL(INVRFF(1:IDEG-1,1:IDEG-1),QTFF(1:IDEG-1,1:IDEG-1))

					end if!for wall only
				end if!stencils
			end if!for test
		end if!ggs
	END DO
	
	deallocate(LSQM,QFF,RFF,QTFF,INVRFF)
	deallocate(VELLSQMAT)
	deallocate(LSCQM)
	deallocate(INTBS,BASEFACEVAL,BASEFACGVAL,PERMUTATION,PERMUTATIONG,xder,yder,zder)
	deallocate(WLSQR)
	deallocate(stencil)
	deallocate(invmat)
	deallocate(STENCILS)

END SUBROUTINE PRESTORE_RECONSTRUCTION_MovingMesh_2D


 



SUBROUTINE LOCALISE_STENCIL_MovingMesh_2D(N,Iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z,node_position_index)
	!> @brief
	!> This subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::n,iconsi,node_position_index
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	INTEGER::I,J,l,IXFF,IXSST,ikg,ismp,inv,ineedt,ikg2,IN_STEN,K,itarget,NOJCOUNT
	REAL::RIN_STEN

	INEEDT=IRECEXR(1)%TOT
	i=iconsi
	IKG2=0

	DO ISMP=1,TYPESTEN
		IKG=0
		if ((ees.ne.5).or.(ismp.eq.1))then
			itarget=ielem(n,i)%iNUMNEIGHBOURS                    
		else
			itarget=NUMNEIGHBOURS2
		end if 
		DO L=1,itarget
			IF (ILOCALSTENCIL(N,I,ISMP,L).GT.0)THEN
				IKG=IKG+1
			END IF
		END DO
		IF (IKG.EQ.itarget)THEN
			if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
				IKG2=IKG2+1
				DO L=1,itarget
					ILOX_IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
					j=(XMPIL(ILOX_IHEXG(IKG2,L)))
					ILOX_IHEXL(IKG2,L)=J
					ILOX_ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
					CALL COMPUTE_CENTRE_MovingMesh_2d(j,cords,node_position_index)
					ILOX_XXC(IKG2,L)=CORDS(1)
					ILOX_YYC(IKG2,L)=CORDS(2)
					
					ILON_NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)
					DO K=1,ielem(n,j)%nonodes    
						! ILON_x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1)
						! ILON_y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2)
						ILON_x(IKG2,L,K)=local_nodes(ielem(n,j)%nodes(K))%positions(node_position_index,1)
						ILON_y(IKG2,L,K)=local_nodes(ielem(n,j)%nodes(K))%positions(node_position_index,2)
					END DO
				END DO
			else
				IKG2=IKG2+1
				DO L=1,itarget
					ILOX_IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
					ILOX_IHEXB(IKG2,L)=XMPIE(ILOCALSTENCIL(N,I,ISMP,L))
					IF (ILOX_IHEXB(IKG2,L).EQ.N)THEN
						j=(XMPIL(ILOX_IHEXG(IKG2,L)))
						IF (ILOX_IHEXG(IKG2,L).EQ.IELEM(N,J)%IHEXGL)THEN
							ILOX_IHEXL(IKG2,L)=J
							ILOX_ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
						END IF
						CALL COMPUTE_CENTRE_MovingMesh_2d(j,cords,node_position_index)
						ILOX_XXC(IKG2,L)=CORDS(1)
						ILOX_YYC(IKG2,L)=CORDS(2)

						ILON_NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)
						DO K=1,ielem(n,j)%nonodes    
							! ILON_x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1)
							! ILON_y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2)
							ILON_x(IKG2,L,K)=local_nodes(ielem(n,j)%nodes(K))%positions(node_position_index,1)
							ILON_y(IKG2,L,K)=local_nodes(ielem(n,j)%nodes(K))%positions(node_position_index,2)
						END DO
					ELSE
						DO IXFF=1,INEEDT
							IF (IEXCORDR(IXFF)%PROCID.EQ.ILOX_IHEXB(IKG2,L))THEN
								DO IXSST=1,IRECEXR(IXFF)%MUCHINEED(1)
									IF ((ILOX_IHEXG(IKG2,L)).EQ.IRECEXR1(IXFF)%WHATINEED(IXSST))THEN
										ILOX_IHEXL(IKG2,L)=IXSST
										ILOX_IHEXN(IKG2,L)=IXFF
										ILOX_ISHAPE(IKG2,L)=IRECEXR1(IXFF)%ISHAPE(IXSST)
										EXIT
									END IF
								END DO
							END IF
						END DO

						SELECT CASE(ILOX_ISHAPE(IKG2,L))
						CASE(5)
							IN_STEN=4
							RIN_STEN=4.D0

						CASE(6)
							IN_STEN=3
							RIN_STEN=3.D0
						END SELECT
						ILON_X(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,1)
						ILON_Y(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,2)
						
						ILOX_XXC(IKG2,L)=sum(ILON_X(IKG2,L,1:IN_STEN))/RIN_STEN
						ILOX_YYC(IKG2,L)=sum(ILON_Y(IKG2,L,1:IN_STEN))/RIN_STEN
						
					END IF
				END DO
			END IF
		END IF
	END DO	

END SUBROUTINE LOCALISE_STENCIL_MovingMesh_2D






SUBROUTINE  LOCALISE_STENCIL_STEP2_MovingMesh_2D(N, I, ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z, node_position_index)
	!> @brief
	!> This subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
	IMPLICIT NONE
	INTEGER::J,K,L,KK,PRK,JJ,kmaxe,ineedt,jx2,jx,in1,facexx,ixxfff,IHGT,IHGJ,ITARGET,IDUM,IN_STEN,NJ,ELEM_DEC,ELtype
	integer,intent(in)::n, i, node_position_index
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
	INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X
	real,dimension(2)::tempcentres,rel2
	real::dumv1,dumv2,detjc,dist1,DISTFD,X1X,Y1Y,X2X,Y2Y
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	INTEGER,DIMENSION(4)::NOJCOUNT
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
	REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
	REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
	REAL,DIMENSION(1:dimensiona,1:dimensiona)::VVA1
	REAL,DIMENSION(1)::DETA

	KMAXE=XMPIELRANK(N)
	INEEDT=IRECEXR(1)%TOT

	IF (IPERIODICITY.EQ.1)THEN
	  	DO JJ=1,IELEM(N,I)%ADMIS
	  		if ((EES.ne.5).or.(jj.eq.1))then
                itarget=ielem(n,i)%iNUMNEIGHBOURS            
            else
                itarget=NUMNEIGHBOURS2
            end if
	    	DO J=2,Itarget
	    		IF(ABS(ILOX_XXC(JJ,J)-ILOX_XXC(JJ,1)).GT.XPER*oo2)THEN
		    		ILOX_XXC(JJ,J)=ILOX_XXC(JJ,J)+(XPER*SIGN(1.0,ILOX_XXC(JJ,1)-XPER*oo2))
		    		DO KK=1,4
		    			ILON_X(JJ,J,KK)=ILON_X(JJ,J,KK)+(XPER*SIGN(1.0,ILOX_XXC(JJ,1)-XPER*oo2))
		    		END DO
	    		END IF
	    		IF(ABS(ILOX_YYC(JJ,J)-ILOX_YYC(JJ,1)).GT.YPER*oo2)THEN
		    		ILOX_YYC(JJ,J)=ILOX_YYC(JJ,J)+(YPER*SIGN(1.0,ILOX_YYC(JJ,1)-YPER*oo2))
		    		DO KK=1,4
		    			ILON_Y(JJ,J,KK)=ILON_Y(JJ,J,KK)+(YPER*SIGN(1.0,ILOX_YYC(JJ,1)-YPER*oo2))
		    		END DO
	    		END IF
	    	END DO
		END DO
	END IF
	
	!$OMP MASTER

		IF (CODE_PROFILE.EQ.30)THEN
			DO NJ=1,IELEM(N,I)%nonodes
				X1X=ILON_x(1,1,NJ)
				Y1Y=ILON_Y(1,1,NJ)
				NOJCOUNT(NJ)=0

				IF (ILOCAL_RECON3(I)%LOCAL.eq.1)then
					DO L=2,itarget
						j=(XMPIL(ILOX_IHEXG(1,L)))
						DO K=1,ielem(n,j)%nonodes    
							X2X=ILON_x(1,L,K)
							Y2Y=ILON_y(1,L,K)

							DISTFD=SQRT(((X1X-X2X)**2)+((Y1Y-Y2Y)**2))

							IF (DISTFD.LT.TOLSMALL)THEN
								NOJCOUNT(NJ)=NOJCOUNT(NJ)+1
								IELEM(N,I)%NODES_NEIGHBOURS(NJ,NOJCOUNT(NJ))=L
							END IF
						END DO
					END DO
				END IF
				!---------------------- MIXED----------!
				IF (ILOCAL_RECON3(I)%LOCAL.NE.1)then
					DO L=2,itarget
						IF (ILOX_IHEXB(1,L).EQ.N)THEN
							j=(XMPIL(ILOX_IHEXG(1,L)))
							DO K=1,ielem(n,j)%nonodes    
								X2X=ILON_x(1,L,K)
								Y2Y=ILON_y(1,L,K)

								DISTFD=SQRT(((X1X-X2X)**2)+((Y1Y-Y2Y)**2))

								IF (DISTFD.LT.TOLSMALL)THEN
									NOJCOUNT(NJ)=NOJCOUNT(NJ)+1
									IELEM(N,I)%NODES_NEIGHBOURS(NJ,NOJCOUNT(NJ))=L
								END IF
							END DO
						ELSE
							SELECT CASE(ILOX_ISHAPE(1,L))
							CASE(5)
							IN_STEN=4
							CASE(6)
							IN_STEN=3
							END SELECT
							DO K=1,IN_STEN
								X2X=ILON_X(1,L,K)
								Y2Y=ILON_Y(1,L,K)

								DISTFD=SQRT(((X1X-X2X)**2)+((Y1Y-Y2Y)**2))

								IF (DISTFD.LT.TOLSMALL)THEN
									NOJCOUNT(NJ)=NOJCOUNT(NJ)+1
									IELEM(N,I)%NODES_NEIGHBOURS(NJ,NOJCOUNT(NJ))=L
								END IF

							END DO
						END IF
					END DO
				END IF

			END DO
			! WRITE(630+N,*)"ELEMENT NUMBER GLOBAL",IELEM(N,I)%IHEXGL
			! ALLOCATE(IELEM(N,I)%NOJECOUNT(IELEM(N,I)%nonodes))
			
			! DO NJ=1,IELEM(N,I)%nonodes
			! 	WRITE(630+N,*)"NODE NUMBER",NJ
			! 	IELEM(N,I)%NOJECOUNT(NJ)=NOJCOUNT(NJ)

			! 	DO J=1,NOJCOUNT(NJ)
			! 	!IF (IELEM(N,I)%NODES_NEIGHBOURS(NJ,J).GT.0)THEN
			! 		WRITE(630+N,*)J,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)
			! 	!END IF
			! 	END DO
			! END DO

		END IF
	!$OMP END MASTER
  
    VEXT=0.0d0
    NODES_LIST=0.0d0
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=0.0d0
       
    jx=IELEM(N,I)%NONODES
	do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILON_x(1,1,K)
	    NODES_LIST(k,2)=ILON_y(1,1,K)
	   
	    VEXT(K,:)=NODES_LIST(k,:)
	END DO
	CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
    
    SELECT CASE(ielem(n,i)%ishape)

      CASE(5)
      	if (IELEM(N,I)%MODE.eq.0)then
			CALL QUADRATUREQUAD(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
			DUMV1=QUADVOLUME(N,VEXT,QPOINTS,WEQUA3D)
			CALL COMPUTE_CENTRE_MovingMesh_2D(i,cords,node_position_index)
			vext(1,1:dims)=cords(1:dims)
      	else
			VEXT(1:3,1:2)=ELEM_LISTD(1,1:3,1:2)
			! DUMV1=TRIANGLEVOLUME(N)
	  		call COMPUTEJACOBIANS2(N,VEXT,VVA1,DETA)
		end if
      
      CASE(6)
      	VEXT(1:3,1:2)=ELEM_LISTD(1,1:3,1:2)
		! DUMV1=TRIANGLEVOLUME(N)
		call COMPUTEJACOBIANS2(N,VEXT,VVA1,DETA)

    END SELECT

	ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2)=VVa1(1:2,1:2)
		  	    
 	IDUM=0
	if (ielem(n,i)%interior.eq.1)then
		DO j=1,IELEM(N,I)%IFCA
			if (ielem(n,i)%ibounds(J).gt.0)then
				if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
					IDUM=1
				end if
			END IF
		END DO
	end if
		   
	detjc=deta(1)
	ILOCAL_RECON3(I)%VEXT_REF(1:2)=VEXT(1,1:2)
		    
	if ((poly.eq.4).OR.(DG.EQ.1))then
		ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2)=0.0d0
		ILOCAL_RECON3(I)%INVCCJAC(1,1)=1.0d0
		ILOCAL_RECON3(I)%INVCCJAC(2,2)=1.0d0
		detjc=1.0
		VEXT(1,1)=ielem(n,i)%xxc
		VEXT(1,2)=ielem(n,i)%yyc
		ILOCAL_RECON3(I)%VEXT_REF(1:2)=VEXT(1,1:2)
	end if
		    	    
	DO JJ=1,IELEM(N,I)%ADMIS
        if ((EES.ne.5).or.(jj.eq.1))then
			itarget=ielem(n,i)%iNUMNEIGHBOURS
		else
			itarget=NUMNEIGHBOURS2
		end if
		  
		DO L=1,Itarget
			! TEMPCENTRES=ZERO
			! TEMPCENTRES(1)=ILOX_XXC(JJ,L)
			! TEMPCENTRES(2)=ILOX_YYC(JJ,L)
			! TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
			
			! ILOX_XXC(JJ,L)=TEMPCENTRES(1)
			! ILOX_YYC(JJ,L)=TEMPCENTRES(2)
				
			DO KK=1,4
				TEMPCENTRES=ZERO
				TEMPCENTRES(1)=ILON_X(JJ,L,KK)
				TEMPCENTRES(2)=ILON_Y(JJ,L,KK)
				
				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
				
				ILON_X(JJ,L,KK)=TEMPCENTRES(1)
				ILON_Y(JJ,L,KK)=TEMPCENTRES(2)
			END DO
				
		END DO	
	END DO			

	DO JJ=1,IELEM(N,I)%ADMIS
		if ((EES.ne.5).or.(jj.eq.1))then
			itarget=ielem(n,i)%iNUMNEIGHBOURS	
		else
			itarget=NUMNEIGHBOURS2	
		end if
		DO L=1,Itarget
			if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
				ILOX_VOLUME(JJ,L)=(IELEM(N,ILOX_IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				if ((EES.ne.5).or.(jj.eq.1))then
				 	if (idum.eq.1)then
						ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
					else
						ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
					end if
				
					ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
					ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				else
					ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
					ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
					ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				end if
			Else
				IF (ILOX_IHEXB(JJ,L).eq.N)THEN
					ILOX_VOLUME(JJ,L)=(IELEM(N,ILOX_IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				 	if ((EES.ne.5).or.(jj.eq.1))then
						if (idum.eq.1)then
							ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
						else
							ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
						end if
						ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
						ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
						ILOCAL_RECON3(I)%IHEXB(JJ,L)=ILOX_IHEXB(JJ,L)
					else
						ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
						ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
						ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
						ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ILOX_IHEXB(JJ,L)
					end if
				else 
				 	if ((EES.ne.5).or.(jj.eq.1))then
						ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
						ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
						ILOCAL_RECON3(I)%IHEXB(JJ,L)=ILOX_IHEXB(JJ,L)
						ILOCAL_RECON3(I)%IHEXN(JJ,L)=ILOX_IHEXN(JJ,L)
					else
						ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
						ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
						ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ILOX_IHEXB(JJ,L)
						ILOCAL_RECON3(I)%IHEXNc(JJ,L)=ILOX_IHEXN(JJ,L)
					end if
				
				    ELTYPE=ILOX_ISHAPE(jj,L)
				    ELEM_LISTD=0.0d0; VEXT=0.0d0; NODES_LIST=0.0d0
				    select case(ELTYPE)
				      case(5)
				      	ELEM_DEC=2; jx=4
				      case(6)
				    	ELEM_DEC=1; jx=3
				    end select

					do K=1,jx
						NODES_LIST(k,1)=ILON_x(jj,l,K)
						NODES_LIST(k,2)=ILON_y(jj,l,K)
						
						VEXT(K,:)=NODES_LIST(k,:)
					END DO
					CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
					dumv2=0.0d0
					do k=1,ELEM_DEC
					    VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
					    dumv2=dumv2+TRIANGLEVOLUME(N,vext)
					end do
					ILOX_VOLUME(JJ,L)=dumv2
					if ((EES.ne.5).or.(jj.eq.1))then
						if (idum.eq.1)then
							ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
						else
							ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
						end if
                    else
                        ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)                    
                    end if	

				end if
			END IF
		END DO	
	END DO

	if (ielem(n,i)%interior.eq.0)then
	  	CALL COMPUTE_CENTRE_MovingMesh_2D(i, cords, node_position_index)
	  	vext(1,1:dims)=cords(1:dims)
	  	do k=1,ielem(n,i)%ifca
		  	j=IELEM(N,i)%INEIGH(K)
		  	CALL COMPUTE_CENTRE_MovingMesh_2D(j,cords, node_position_index)
		    vext(2,1:dims)=cords(1:dims)
		    dist1=distance2(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    	IELEM(N,i)%DIH(K)=dist1
		    end if
	  	end do
    else
		CALL COMPUTE_CENTRE_MovingMesh_2d(i,cords, node_position_index)
		vext(1,1:dims)=cords(1:dims)
		do k=1,ielem(n,i)%ifca
			if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
				facexx=k
		  
				IXXFFf=2
		  
				call COMPUTE_CENTRE_MovingMesh_2dF(N,i,facexx,IXXFFf,CORDS,node_position_index)
				VEXT(2,1:dims)=cords(1:dims)
				
				dist1=distance2(n,vext)
				IF (RUNGEKUTTA.ge.2)THEN
					IELEM(N,i)%DIH(K)=dist1*2.0d0
				end if
			end if
			if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
				if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
					j=IELEM(N,i)%INEIGH(K)
					CALL COMPUTE_CENTRE_MovingMesh_2D(j,cords, node_position_index)
					vext(2,1:dims)=cords(1:dims)
					dist1=distance2(n,vext)
					IF (RUNGEKUTTA.ge.2)THEN
						IELEM(N,i)%DIH(K)=dist1
					end if
				else						!from another cpu 
		    		DO In1=1,IELEM(N,I)%iNUMNEIGHBOURS
			  			IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  			IELEM(N,i)%INDEXI(K)=In1
				    		IF (RUNGEKUTTA.ge.2)THEN
		    					vext(2,1)=ILOX_XXC(1,In1);vext(2,2)=ILOX_yyC(1,In1)
		     					dist1=distance2(n,vext)
		    					IELEM(N,i)%DIH(K)=dist1
				      		end if
			  			end if
		    		end do
				end if
			end if
			if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
				if (ielem(n,i)%ineighb(k).eq.n)then	
					j=IELEM(N,i)%INEIGH(K)
					CALL COMPUTE_CENTRE_MovingMesh_2D(j,cords, node_position_index)
					vext(2,1:dims)=cords(1:dims)  
					IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
						vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
					end if
					IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
						vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
					end if
		    
		    		dist1=distance2(n,vext)
		    		IF (RUNGEKUTTA.ge.2)THEN
		    			IELEM(N,i)%DIH(K)=dist1
		    		end if
				else	!periodic boundaries from another cpu
					DO In1=1,IELEM(N,I)%iNUMNEIGHBOURS
						IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
							IELEM(N,i)%INDEXI(K)=In1
							IF (RUNGEKUTTA.ge.2)THEN
								vext(2,1)=ILOX_XXC(1,In1);vext(2,2)=ILOX_yyC(1,In1);
								
								IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
									vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
								end if
								IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    						vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    					end if
		    
		    					dist1=distance2(n,vext)
		    					IELEM(N,i)%DIH(K)=dist1
							end if
			  			end if
		    		end do
				end if
			end if
	  	end do
	end if

    DO JJ=1,IELEM(N,I)%ADMIS
		if ((EES.ne.5).or.(jj.eq.1))then
            itarget=ielem(n,i)%iNUMNEIGHBOURS
        else
            itarget=NUMNEIGHBOURS2
        end if
		DO L=1,Itarget
			TEMPCENTRES=ZERO
			TEMPCENTRES(1)=ILOX_XXC(JJ,L)
			TEMPCENTRES(2)=ILOX_YYC(JJ,L)
			TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
			
			ILOX_XXC(JJ,L)=TEMPCENTRES(1)
			ILOX_YYC(JJ,L)=TEMPCENTRES(2)	
		END DO	
	END DO

END SUBROUTINE LOCALISE_STENCIL_STEP2_MovingMesh_2D









END MODULE PRESTORE_MOVINGMESH
