MODULE PRESTORE
USE DECLARATION
USE DERIVATIVES
USE LIBRARY
USE BASIS
USE TRANSFORM
USE LOCAL
USE LAPCK
USE DG_FUNCTIONS
IMPLICIT NONE


 CONTAINS
 
 
 SUBROUTINE PRESTORE_1(N)
!> @brief
!> This subroutine calls other subroutines to prestore the pseudoinverse reconstruction least square matrices
IMPLICIT NONE
INTEGER,INTENT(IN)::N
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
!$OMP DO
DO I=1,KMAXE

CALL LOCALISE_STENCIL(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)


CALL LOCALISE_STEN2(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
call CHECKGRADS(N,I)
CALL FIND_ROT_ANGLES(N,I)

CALL PRESTORE_RECONSTRUCTION3(N,i,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)

 if ((dg.eq.1).OR.(ADDA_type.EQ.2))then
 
 CALL PRESTORE_DG1(i)
    end if
 

END DO
!$OMP END DO




else





!$OMP DO
DO I=1,KMAXE
CALL LOCALISE_STENCIL2d(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
CALL LOCALISE_STEN2d(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
 call CHECKGRADS2d(N,I)
CALL FIND_ROT_ANGLES2D(N,I)
 CALL PRESTORE_RECONSTRUCTION2(N,i,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
 if ((dg.eq.1).OR.(ADDA_type.EQ.2))then
 
 CALL PRESTORE_DG1(i)
 
 
    end if
END DO
!$OMP END DO






end if


deallocate(ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y)

	IF (DIMENSIONA.EQ.3)THEN
	deallocate(ILOX_ZZC,ILON_Z)
	end if




END SUBROUTINE PRESTORE_1






SUBROUTINE PRESTORE_RECONSTRUCTION3(N,iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
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
                

iconsidered=i


		
        !firstly compute the integrals of the basis functions
		
		INTBS=zero;JXX=1;IXX=i;LXX1=1;number_of_dog=ielem(n,i)%idegfree;kxx=ielem(n,i)%iorder;ELTYPE=ielem(n,i)%ishape
		icompwrt=0

		INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)

		INTEG_BASIS(I)%VALUE(1:ielem(n,i)%IDEGFREE)=INTBS(1:ielem(n,i)%IDEGFREE)
                
        !the indicator matrix for the smaller polynomial
		IF (IWENO.EQ.1)THEN
		CALL INDICATORMATRIX(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
		END IF
		
					!the indicator matrix for the lower order polynomial
					if (ees.eq.5)then
					INTBS=zero;JXX=1;IXX=i;LXX1=1;number_of_dog=idegfree2;kxx=IORDER2;ELTYPE=ielem(n,i)%ishape
					icompwrt=1
					INTBS=CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
					INTEG_BASIS(I)%VALUEc(1:number_of_dog)=INTBS(1:number_of_dog)
						IF (IWENO.EQ.1)THEN
						CALL INDICATORMATRIX2(N,I,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
						icompwrt=0
						END IF
					end if

		
		LLCO=IELEM(N,I)%ADMIS;IMAX=IELEM(N,I)%inumneighbours-1;INUM=IELEM(N,I)%inumneighbours;IDEG=IELEM(N,I)%iDEGFREE
	   INUMO=ielem(n,i)%iorder;imax2=numneighbours2-1;inum2=numneighbours2;ideg2=iDEGFREE2;inumo2=iorder2


	   !this section computes the ratio of volumes inside each stencil
                DIST_STEN=ZERO
				
				DIST_STEN=-tolbig; DIST_STEN2=tolbig
	       	DO LL=1,ielem(n,i)%iNUMNEIGHBOURS
				if (ILOX_VOLUME(1,LL).lt.DIST_STEN2)then
				DIST_STEN2=ILOX_VOLUME(1,LL)
				end if
				if (ILOX_VOLUME(1,LL).gt.DIST_STEN)then
				DIST_STEN=ILOX_VOLUME(1,LL)
				end if
			end do
			
			IELEM(N,I)%STENCIL_DIST=MAX((DIST_STEN/DIST_STEN2),(DIST_STEN2/DIST_STEN))
				
		


		!now we start with looping all the admissible stencils
		DO LL=1,LLCO	!for all stencils
									if((ees.ne.5).or.(ll.eq.1))then
										IMAX=IELEM(N,I)%inumneighbours-1;INUM=IELEM(N,I)%inumneighbours;IDEG=IELEM(N,I)%iDEGFREE;INUMO=ielem(n,i)%iorder
										icompwrt=0;number_of_dog=IDEG
									else
										imax=numneighbours2-1;inum=numneighbours2;ideg=iDEGFREE2;inumo=IORDER2;number_of_dog=IDEG
										icompwrt=1
									end if

				DO K=1,imax	!for all neighbours
								ixx=i; kxx=INUMO


								X_STENCIL=(ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1))**2
								Y_STENCIL=(ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1))**2

								Z_STENCIL=(ILOX_ZZC(ll,k+1)-ILOX_ZZC(ll,1))**2

								DIST_STEN2=SQRT(X_STENCIL+Y_STENCIL+Z_STENCIL)

									IF (WEIGHT_LSQR.EQ.1)THEN
									WLSQR(ll,K)=1.0D0/((DIST_STEN2))
									ELSE
									WLSQR(ll,K)=1.0D0
									END IF






					if (fastest.eq.1)then	!this is when transformation is not active (Rarely used)
					x1 = ILOX_XXC(ll,k+1)-ILOX_XXC(ll,1)
					y1 = ILOX_YYC(ll,k+1)-ILOX_YYC(ll,1)
					z1 = ILOX_ZZC(ll,k+1)-ILOX_ZZC(ll,1)


					IF (GREENGO.EQ.0)THEN	!for the least squares green gauss gradients
							if (idum.eq.1)then
								if((ees.ne.5).or.(ll.eq.1))then
								icompwrt=0
								ILOCAL_RECON3(I)%STENCILS(LL,K,1:ielem(n,i)%idegfree)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,ielem(n,i)%iorder,I,ielem(n,i)%idegfree,icompwrt)
								ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
								else
								icompwrt=1
								ILOCAL_RECON3(I)%STENCILSc(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,I,IDEG,icompwrt)
								ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
								end if
							else
								if((ees.ne.5).or.(ll.eq.1))then

								icompwrt=0
								STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)
								else
								icompwrt=1
								STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)

								end if
							end if

					ELSE
							if((ees.ne.5).or.(ll.eq.1))then
							icompwrt=0
							STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)
							else
							icompwrt=1
							STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*basis_rec(N,x1,y1,z1,INUMO,IXX,IDEG,icompwrt)


							end if
					END IF



					else		!with coordinate transformation


					IXX=i;jxx=k+1;lxx1=ll
					ELTYPE=ILOX_ishape(ll,k+1)
						IF (GREENGO.EQ.0)THEN
								if (idum.eq.1)then	!smaller memory footprint for non boundary elements
										if((ees.ne.5).or.(ll.eq.1))then
										icompwrt=0




										ILOCAL_RECON3(I)%STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
										ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)
										else
										icompwrt=1
										ILOCAL_RECON3(I)%STENCILSc(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
										ilocal_recon3(i)%WEIGHTL(ll,k)=WLSQR(ll,K)

										end if
								else
										if((ees.ne.5).or.(ll.eq.1))then
										icompwrt=0
										STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
										else
										icompwrt=1
										STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)

										end if
								end if
						ELSE
								if((ees.ne.5).or.(ll.eq.1))then
							icompwrt=0
								STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
								else
								icompwrt=1
								STENCILS(LL,K,1:IDEG)=WLSQR(ll,K)*COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)

								end if

						END IF


					end if




		END DO	!close loop for all neighbours








		IF ((GREENGO.EQ.0))THEN
				if (idum.eq.1)then
					invmat=zero;LSCQM=zero
					DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
					if((ees.ne.5).or.(ll.eq.1))then
					LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
					+((ILOCAL_RECON3(I)%STENCILS(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILS(LL,LQ,IQ)))
					else
					LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
					+((ILOCAL_RECON3(I)%STENCILSc(LL,LQ,JQ)*ILOCAL_RECON3(I)%STENCILSc(LL,LQ,IQ)))
					end if
					END DO;END DO;END DO
				else
					invmat=zero;LSCQM=zero
					DO IQ=1,IDEG;DO JQ=1,IDEG;DO LQ=1,IMAX
					LSCQM(JQ,IQ)=LSCQM(JQ,IQ)&
					+((STENCILS(LL,LQ,JQ)*STENCILS(LL,LQ,IQ)))
					END DO;END DO;END DO
				end if
		ELSE

				invmat=zero;LSCQM=zero
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
				if((ees.ne.5).or.(ll.eq.1))then
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

		if((ees.ne.5).or.(ll.eq.1))then

! 		call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
! 		stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL),IDEG)


		ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))



					do iq=1,imax
					ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)&
					*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
					end do
		else
! 					call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
! 				stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stenciltc(1:IDEG,1:IMAX,LL),IDEG)


				ILOCAL_RECON3(i)%invmat_stenciltC(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))



					do iq=1,imax
					ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stenciltc(:,iq,LL)&
					*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
					end do
		end if


	
		if (ielem(n,i)%ggs.ne.1)then	!ggs

		IF (ITESTCASE.EQ.4)THEN	!TEST

	 

	 

IF (LL.EQ.1)THEN		!stencils
! 

  
  if (idum.eq.1)then		!for wall only
    
	 
	
!     	
    DO IHGT=1,3; DO IHGJ=1,3
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
								  CALL coordinates_face_inner(N,Iconsidered,facex,vext,NODES_LIST)


								   if (ielem(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
                                            N_NODE=4
                                    else
                                            N_NODE=3
                                    end if



								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    AY=cords(2)
								    AX=cords(1)
								    AZ=cords(3)
				
					
					    VEXT(1,1)=AX;VEXT(1,2)=AY;VEXT(1,3)=AZ
					    VEXT(1,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(1,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  
					    AX=VEXT(1,1);AY=VEXT(1,2);AZ=VEXT(1,3)
				
				
				
				
				
				
				ANGLE1=IELEM(N,I)%FACEANGLEX(j)
				ANGLE2=IELEM(N,I)%FACEANGLEY(j)
				
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
				NNX=(NX*AINVJT(1,1))+(NY*AINVJT(2,1))+(NZ*AINVJT(3,1))
				NNY=(NX*AINVJT(1,2))+(NY*AINVJT(2,2))+(NZ*AINVJT(3,2))
				NNZ=(NX*AINVJT(1,3))+(NY*AINVJT(2,3))+(NZ*AINVJT(3,3))
				
				DO IQ=1, IDEG
				IF (POLY.EQ.1) THEN
				XDER(IQ)=DFX(AX,AY,AZ,IQ,I);  YDER(IQ)=DFY(AX,AY,AZ,IQ,I);  ZDER(IQ)=DFZ(AX,AY,AZ,IQ,I)
				END IF
				IF (POLY.EQ.2) THEN
				XDER(IQ)=DLX(AX,AY,AZ,IQ,I);  YDER(IQ)=DLY(AX,AY,AZ,IQ,I);  ZDER(IQ)=DLZ(AX,AY,AZ,IQ,I)
				END IF
				IF (POLY.EQ.4) THEN
				XDER(IQ)=TL3DX(AX,AY,AZ,IQ,I);  YDER(IQ)=TL3DY(AX,AY,AZ,IQ,I);  ZDER(IQ)=TL3DZ(AX,AY,AZ,IQ,I)
				END IF
				
				END DO

				icompwrt=0

				BASEFACEVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,ielem(n,i)%Iorder,I,ielem(n,i)%IDEGFREE,icompwrt)

				!if (thermal.eq.0)then
				BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=((NNX*XDER(1:ielem(n,i)%IDEGFREE))+(NNY*YDER(1:ielem(n,i)%IDEGFREE))+(NNZ*ZDER(1:ielem(n,i)%IDEGFREE)))
				!ELSE
				!icompwrt=0
				!BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,ielem(n,i)%Iorder,I,ielem(n,i)%IDEGFREE,icompwrt)

! 				end if


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

!  

! now swap basis functions and thus coefficients
PERMUTATION(1)=K0; PERMUTATION(K0)=1; SSSS=BASEFACEVAL(1)
BASEFACEVAL(1)=BASEFACEVAL(K0); BASEFACEVAL(K0)=SSSS
PERMUTATIONG(1)=G0; PERMUTATION(G0)=1; GGGG=BASEFACGVAL(1)
BASEFACGVAL(1)=BASEFACGVAL(G0); BASEFACGVAL(G0)=GGGG
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
Qff=ZERO; Rff=ZERO; QTff=ZERO; INVRff=ZERO
	  IVGT=IDEG
	  
    CALL QRDECOMPOSITION(LSCQM,Qff,Rff,IVGT-1)
    CALL TRANSPOSEMATRIX(Qff,QTff,IVGT-1)
    
    CALL INVERT(Rff,INVRff,IVGT)
!   final inverted R^(-1)*Q^(-1)
ILOCAL_RECON3(I)%TEMPSQMAT(1:IDEG-1,1:IDEG-1) =&
 MATMUL(INVRff(1:IDEG-1,1:IDEG-1),QTff(1:IDEG-1,1:IDEG-1))

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

		Qff=ZERO; Rff=ZERO; QTff=ZERO; INVRff=ZERO
			IVGT=IDEG
			CALL QRDECOMPOSITION(LSCQM,Qff,Rff,IVGT-1)
			CALL TRANSPOSEMATRIX(Qff,QTff,IVGT-1)

			CALL INVERT(Rff,INVRff,IVGT)
		!   final inverted R^(-1)*Q^(-1)
		ILOCAL_RECON3(I)%VELINVLSQMAT(1:IDEG-1,1:IDEG-1) = &
		MATMUL(INVRff(1:IDEG-1,1:IDEG-1),QTff(1:IDEG-1,1:IDEG-1))
			
		
		
		
		
		
		
		end if!for wall only
		end if!stencils
		end if!for test
		end if!ggs
		end do !for all stencils
		

		DEALLOCATE(LSQM,QFF,RFF,QTFF,INVRFF)	
		DEALLOCATE(VELLSQMAT)
		DEALLOCATE(LSCQM)	
		DEALLOCATE(INTBS,BASEFACEVAL,BASEFACGVAL,PERMUTATION,PERMUTATIONG,xder,yder,zder)	
		DEALLOCATE(WLSQR)	
		DEALLOCATE(stencil)	
		DEALLOCATE(invmat)	
		DEALLOCATE(STENCILS)	




END SUBROUTINE PRESTORE_RECONSTRUCTION3


subroutine walls_higher(n)
!> @brief
!> This subroutine allocates the memory for wall bounded cells for their constrained least squares reconstruction
implicit none
integer,intent(in)::n
integer::i,j,k,imax,ideg,inumo,inum,kmaxe,idum,idum2

KMAXE=XMPIELRANK(N)

DO I=1,KMAXE

  ielem(n,i)%walls=0
 IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder

if (ielem(n,i)%ggs.ne.1)then

 IDUM=0;idum2=0
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
	        IDUM=1
		if (dimensiona.eq.3)then
		if (ielem(n,i)%types_faces(j).eq.5)then
		  idum2=idum2+qp_quad_n
		else
		  idum2=idum2+QP_TRIANGLE_n
		end if
		else
		  idum2=idum2+qp_line_n
		end if
	      end if
	  END IF
	END DO
  end if
  if (idum.eq.1)then
    allocate(ielem(n,i)%num_of_wall_gqp(1))
    ielem(n,i)%num_of_wall_gqp(1)=idum2
    ielem(n,i)%walls=1
	    
	    IF (FASTEST.NE.1)THEN
	    ALLOCATE (ILOCAL_RECON3(I)%VELINVLSQMAT(IDEG-1,IDEG-1))
	    ALLOCATE (ILOCAL_RECON3(I)%WALLCOEFF(ideg))
	    ALLOCATE (ILOCAL_RECON3(I)%VELLSQ(IMAX,IDEG-1))
	    ALLOCATE (ILOCAL_RECON3(I)%TEMPSQMAT(IDEG-1,IDEG-1))
	    ALLOCATE (ILOCAL_RECON3(I)%WALLCOEFG(ideg))
	    ALLOCATE (ILOCAL_RECON3(I)%TEMPSQ(IMAX,IDEg-1))
	  
	    ILOCAL_RECON3(I)%VELINVLSQMAT=zero
	    ILOCAL_RECON3(I)%WALLCOEFF=zero
	    ILOCAL_RECON3(I)%VELLSQ=zero
	    ILOCAL_RECON3(I)%WALLCOEFG=zero
	    ILOCAL_RECON3(I)%TEMPSQMAT=zero
	    ILOCAL_RECON3(I)%TEMPSQ=zero
	   END IF
	    
  end if
  
  end if

END DO



end subroutine walls_higher


SUBROUTINE PRESTORE_RECONSTRUCTION2(N,iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
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
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X



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




! 
! 
! 

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
!                invmat,                                               &
!                stencil,                                              &
!                ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),             &
!                'N',                                                  & ! transposition flag for invmat
!                'T'                                                   & ! transposition flag for stencil
!             )
            
            
! call DGEMM ('N','T',IDEG,IMAX,IDEG,ALPHA,invmat(1:ideg,1:ideg),IDEG,&
! stencil(1:imax,1:ideg),IMAX,BETA,ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL),IDEG)

            
            ILOCAL_RECON3(i)%invmat_stencilt(1:IDEG,1:IMAX,LL)=MATMUL(invmat(1:ideg,1:ideg),TRANSPOSE(stencil(1:imax,1:ideg)))





           do iq=1,imax
			ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)=ILOCAL_RECON3(I)%invmat_stencilt(:,iq,LL)&
			*ILOX_VOLUME(ll,iq+1)*WLSQR(ll,iq)
			end do




else
! call gemm(                                               &
!                invmat(1:IDEG,1:IDEG),                                               &
!                stencil(1:imax,1:ideg),                                              &
!                ILOCAL_RECON3(I)%invmat_stenciltc(1:ideg,1:imax,LL),             &
!                'N',                                                  & ! transposition flag for invmat
!                'T'                                                   & ! transposition flag for stencil
!             )
            
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
! 

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
								  CALL coordinates_face_inner2D(N,Iconsidered,facex,vext,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    AY=cords(2)
								    AX=cords(1)

					    VEXT(1,1)=AX;VEXT(1,2)=AY;
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

				!if (thermal.eq.0)then
				Icompwrt=0

				BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=((NNX*XDER(1:ielem(n,i)%IDEGFREE))+(NNY*YDER(1:ielem(n,i)%IDEGFREE)))
				!ELSE
				!Icompwrt=0
				!BASEFACGVAL(1:ielem(n,i)%IDEGFREE)=BASIS_REC2D(N,AX,AY,ielem(n,i)%Iorder,Iconsidered,ielem(n,i)%IDEGFREE,ICOMPWRT)

				!end if


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
!   final inverted R^(-1)*Q^(-1)
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
!   final inverted R^(-1)*Q^(-1)
ILOCAL_RECON3(I)%VELINVLSQMAT(1:IDEG-1,1:IDEG-1) = &
MATMUL(INVRFF(1:IDEG-1,1:IDEG-1),QTFF(1:IDEG-1,1:IDEG-1))
			
		
		
		
		
		
		
		end if!for wall only
		end if!stencils
		end if!for test
		end if!ggs
		end do
		
		
DEALLOCATE(LSQM,QFF,RFF,QTFF,INVRFF)
		DEALLOCATE(VELLSQMAT)
		DEALLOCATE(LSCQM)
		DEALLOCATE(INTBS,BASEFACEVAL,BASEFACGVAL,PERMUTATION,PERMUTATIONG,xder,yder,zder)
		DEALLOCATE(WLSQR)
		DEALLOCATE(stencil)
		DEALLOCATE(invmat)
		DEALLOCATE(STENCILS)


END SUBROUTINE PRESTORE_RECONSTRUCTION2


SUBROUTINE INDICATORMATRIX(N,iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine computes the indicator matrices for weno reconstructions
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iconsi
INTEGER::I,J,K,L,M,jx,jx2,IMAX,INUM,IDEG,INUMO,ELTYPE,ELEM_DEC,inump,ICONSIDERED
REAL::VOLTEMP
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
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
REAL,allocatable,DIMENSION(:,:)::WEFF
allocate(weff(1:IDEGFREE,1:IDEGFREE))
i=iconsi




	IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	 INUMO=ielem(n,i)%iorder
	 
	iconsidered=i
	 

	VEXT=ZERO
    NODES_LIST=ZERO
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=ZERO
          ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ZERO
      jx=IELEM(N,I)%NONODES
     
	  if (dimensiona.eq.3)then
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILON_X(1,1,k)
	    NODES_LIST(k,2)=ILON_y(1,1,k)
	    NODES_LIST(k,3)=ILON_z(1,1,k)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
	  
	  else

	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILON_X(1,1,k)
	    NODES_LIST(k,2)=ILON_y(1,1,k)
	    VEXT(K,1:2)=NODES_LIST(k,1:2)
	  END DO
	  CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
	  	  end if

	 SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
      CALL QUADRATUREHEXA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=HEXAVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_hexa
      CASE(2)
      CALL QUADRATURETETRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=tetraVOLUME(N,VEXT)
      inump=qp_tetra
      CASE(3)
      CALL QUADRATUREPYRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=pyraVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_pyra
      CASE(4)
      CALL QUADRATUREPRISM(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=prismVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_prism

      CASE(5)
      CALL QUADRATUREquad(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=quadVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_quad


      CASE(6)
      CALL QUADRATUREtriangle(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=TRIANGLEVOLUME(N,VEXT)
      inump=qp_triangle

      END SELECT      
      
	if (dimensiona.eq.3)then
	IF (ielem(n,i)%mode.EQ.1)THEN
	  do K=1,ELEM_DEC
	      VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	      CALL QUADRATURETETRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
	      voltemp=tetraVOLUME(N,VEXT)
	      inump=qp_tetra
	      CALL WENOTET(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
	      ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
	
	  END DO
	ELSE
	  

	    CALL WENOTET(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
	      ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)






	END IF
	else

	    do K=1,ELEM_DEC
            VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            CALL QUADRATUREtriangle(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
            voltemp=TRIANGLEVOLUME(N,VEXT)
            inump=qp_triangle
           
	    CALL WENOTET2d(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,ideg,iconsidered)
             
	      ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATOR(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
            END DO


	end if

deallocate(weff)


END SUBROUTINE INDICATORMATRIX


SUBROUTINE INDICATORMATRIX2(N,iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine computes the indicator matrices for cweno reconstructions of the lower-order polynomials
IMPLICIT NONE
INTEGER,INTENT(IN)::N,iconsi
INTEGER::I,J,K,L,M,jx,jx2,IMAX,INUM,IDEG,INUMO,ELTYPE,ELEM_DEC,inump,ICONSIDERED
REAL::VOLTEMP
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
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
REAL,allocatable,DIMENSION(:,:)::WEFF
allocate(weff(1:IDEGFREE,1:IDEGFREE))
i=iconsi

	IMAX=numneighbours2-1
	INUM=numneighbours2
	IDEG=iDEGFREE2
	 INUMO=IORDER2
	 
	
	 iconsidered=i

	VEXT=ZERO
    NODES_LIST=ZERO
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=ZERO
          ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ZERO
      jx=IELEM(N,I)%NONODES
     
	  if (dimensiona.eq.3)then
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILON_X(1,1,k)
	    NODES_LIST(k,2)=ILON_y(1,1,k)
	    NODES_LIST(k,3)=ILON_z(1,1,k)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
	  
	  else

	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILON_X(1,1,k)
	    NODES_LIST(k,2)=ILON_y(1,1,k)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
	  	  end if

	 SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
      CALL QUADRATUREHEXA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=HEXAVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_hexa
      CASE(2)
      CALL QUADRATURETETRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=tetraVOLUME(N,VEXT)
      inump=qp_tetra
      CASE(3)
      CALL QUADRATUREPYRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=pyraVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_pyra
      CASE(4)
      CALL QUADRATUREPRISM(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=prismVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_prism

      CASE(5)
      CALL QUADRATUREquad(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=quadVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      inump=qp_quad


      CASE(6)
      CALL QUADRATUREtriangle(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      voltemp=TRIANGLEVOLUME(N,VEXT)
      inump=qp_triangle

      END SELECT      
      
	if (dimensiona.eq.3)then
	IF (ielem(n,i)%mode.EQ.1)THEN
	  do K=1,ELEM_DEC
	      VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	      CALL QUADRATURETETRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
	      voltemp=tetraVOLUME(N,VEXT)
	      inump=qp_tetra
	      CALL WENOTET(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
	      ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
	
	  END DO
	ELSE
	  

	    CALL WENOTET(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
	      ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)






	END IF
	else
            
            do K=1,ELEM_DEC
            VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            CALL QUADRATUREtriangle(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
            voltemp=TRIANGLEVOLUME(N,VEXT)
            inump=qp_triangle

	    CALL WENOTET2d(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,ideg,iconsidered)

	      ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)=ILOCAL_RECON3(I)%INDICATORc(1:IDEG,1:IDEG)+WEFF(1:IDEG,1:IDEG)
            END DO
            
	end if

	
	
deallocate(weff)


END SUBROUTINE INDICATORMATRIX2

	



SUBROUTINE WENOTET(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N,INUMO,INUMP,IDEG,iconsidered
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
REAL,allocatable,DIMENSION(:,:),INTENT(INOUT)::WEFF
WEFF=zero


		IF (POLY.EQ.1)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL PH1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(2)
		CALL PH2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(3)
		CALL PH3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(4)
		CALL PH4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(5)
		CALL PH5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(6)
		CALL PH6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
		
	END SELECT
	END IF
	IF (POLY.EQ.2)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL PL1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(2)
		CALL PL2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(3)
		CALL PL3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(4)
		CALL PL4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(5)
		CALL PL5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)

		CASE(6)
		CALL PL6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,iconsidered)
		
	END SELECT
	END IF
	
	IF (POLY.EQ.4)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL TL3D1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(2)
		CALL TL3D2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(3)
		CALL TL3D3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(4)
		CALL TL3D4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(5)
		CALL TL3D5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(6)
		CALL TL3D6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
		
	END SELECT
	END IF
	
	
	

END SUBROUTINE WENOTET

SUBROUTINE PH1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=ZERO
		
 	DO I=1,IDEG
        	DO J=1,IDEG
			INTEG =ZERO
                scalerx=1.0d0
                
                
       			 DO K=1,inump
	
	    PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

	    INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH1
SUBROUTINE PH2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=ZERO
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =ZERO
  DO K=1,iNUMp
	PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)
	
	INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
  ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH2
SUBROUTINE PH3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
  DO K=1,inump
	  PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
	DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)
	

	INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
  ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH3
SUBROUTINE PH4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
DO K=1,inump
    PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
    DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

    INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH4
SUBROUTINE PH5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
DO I=1,IDEG
DO J=1,IDEG
scalerx=1.0d0
                
	INTEG =zero
	  DO K=1,inump
PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
  DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

		INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
	  ENDDO
	WEFF(I,J)=WEFF(I,J)+INTEG
ENDDO
ENDDO

END SUBROUTINE PH5
SUBROUTINE PH6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
REAL::PH,INTEG,scalerx
integer::i,j,k
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			PH=DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		    DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		    DFX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DFZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DFZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PH6

SUBROUTINE PL1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				   DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				   DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL1
SUBROUTINE PL2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL2
SUBROUTINE PL3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				 PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL3
SUBROUTINE PL4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL4
SUBROUTINE PL5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
		PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		  DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL5
SUBROUTINE PL6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			     PH=DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		    DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		    DLX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  DLZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*DLZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE PL6


SUBROUTINE TL3D1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				   TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				   TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D1
SUBROUTINE TL3D2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D2
SUBROUTINE TL3D3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				 PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D3
SUBROUTINE TL3D4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
				TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D4
SUBROUTINE TL3D5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
		PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		  TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D5
SUBROUTINE TL3D6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			     PH=TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		    TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
		    TL3DX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX5Y(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4Y2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4YZ(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Y2Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y3Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY4Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY5Z(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3YZ2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Y2Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY3Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY4Z2(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2YZ3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXY2Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY3Z3(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DX2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXYZ4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DY2Z4(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DXZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DYZ5(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED) + &
                  TL3DZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),I,ICONSIDERED)*TL3DZ6(QPOINTS(1,K),QPOINTS(2,K),QPOINTS(3,K),J,ICONSIDERED)

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL3D6





SUBROUTINE WENOTET2D(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N,INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
REAL,allocatable,DIMENSION(:,:),INTENT(INOUT)::WEFF
WEFF=zero

 
		IF (POLY.EQ.1)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL P2DH1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(2)
		CALL P2DH2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(3)
		CALL P2DH3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(4)
		CALL P2DH4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(5)
		CALL P2DH5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(6)
		CALL P2DH6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
		
	END SELECT
	
	END IF
	
	
	
	IF (POLY.EQ.4)THEN
           SELECT CASE(Inumo)
		 CASE(1)
		CALL TL2DH1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(2)
		CALL TL2DH2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(3)
		CALL TL2DH3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(4)
		CALL TL2DH4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(5)
		CALL TL2DH5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)

		CASE(6)
		CALL TL2DH6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
		
	END SELECT
	
	END IF
	
	
	

END SUBROUTINE WENOTET2D
! 
SUBROUTINE P2DH1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				   DF2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J,iCONSIDERED)
				   

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
!           		
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH1
SUBROUTINE P2DH2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH2
SUBROUTINE P2DH3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
        	
        	
        	
        	
			INTEG =zero
       			 DO K=1,inump
				 PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH3
SUBROUTINE P2DH4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				DF2DY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH4
SUBROUTINE P2DH5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
		PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
		  DF2DX5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
                 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH5
SUBROUTINE P2DH6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			PH=DF2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
		    DF2DX5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
		    DF2DX6(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX6(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX5Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX5Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DXY5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DXY5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  DF2DY6(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*DF2DY6(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE P2DH6

SUBROUTINE TL2DH1(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
	
				PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				   TL2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
				   

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
!           		
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH1
SUBROUTINE TL2DH2(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
				
				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH2
SUBROUTINE TL2DH3(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
        	
        	
        	
        	
			INTEG =zero
       			 DO K=1,inump
				 PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH3
SUBROUTINE TL2DH4(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
				PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
				TL2DY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH4
SUBROUTINE TL2DH5(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
		PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
		  TL2DX5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)
                 

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH5
SUBROUTINE TL2DH6(N,WEFF,INUMO,INUMP,VOLTEMP,QPOINTS,WEQUA3D,IDEG,ICONSIDERED)
!> @brief
!> This subroutine computes derivative operator for the smoothness indicator
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::INUMO,INUMP,IDEG,ICONSIDERED
REAL,INTENT(IN)::VOLTEMP
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS),INTENT(IN)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::I,J,K
REAL::PH,INTEG,scalerx
WEFF=zero
		
 	DO I=1,IDEG
        	DO J=1,IDEG
        	scalerx=1.0d0
                
			INTEG =zero
       			 DO K=1,inump
			PH=TL2DX(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
		    TL2DX5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX4Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
		    TL2DX6(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX6(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX5Y(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX5Y(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX4Y2(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX3Y3(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DX2Y4(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DXY5(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DXY5(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED) + &
                  TL2DY6(QPOINTS(1,K),QPOINTS(2,K),I,ICONSIDERED)*TL2DY6(QPOINTS(1,K),QPOINTS(2,K),J,ICONSIDERED)

				

				INTEG=INTEG+(PH*WEQUA3D(K)*voltemp)*(scalerx**(i+j-1))
          		 ENDDO
         		WEFF(I,J)=WEFF(I,J)+INTEG
        	ENDDO
        ENDDO

END SUBROUTINE TL2DH6



FUNCTION CALINTBASIS(N,IXX,JXX,KXX,LXX1,number_of_dog,ICOMPWRT,ELTYPE,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine computes basis functions for each element
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog,ICOMPWRT,ELTYPE
INTEGER,INTENT(IN):: IXX,JXX,KXX,LXX1
INTEGER::K
REAL,DIMENSION(1:number_of_dog)::CALINTBASIS
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


! 	kxx=ielem(n,ixx)%iorder
	CALINTBASIS = ZERO
     	
	CALINTBASIS(1:number_of_dog)=COMPBASEL(N,eltype,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
	
   END FUNCTION


FUNCTION COMPBASEL(N,ELTYPE,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine computes basis functions for each element
IMPLICIT NONE
INTEGER,INTENT(IN)::N,lxx1,jxx,ixx,icompwrt,KXX
INTEGER,INTENT(IN)::ELTYPE,number_of_dog
INTEGER::JX,K,ELEM_DEC
real,dimension(1:number_of_dog)::compbasel,s1
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
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
S1=ZERO

 SELECT CASE(ELTYPE)
      CASE(1)
      JX=8;   ELEM_DEC=6

		 if ((JXX.eq.1))then
		  
		      do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILON_Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		  if (ielem(n,ixx)%mode.eq.0)then
		  
		  COMPBASEL=COMPBASHEX(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)

		    else
		  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		  end do
		    COMPBASEL=S1
		    end if
		  
		  ELSE

		   do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILON_Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
			
		      END DO
		      
		      
		  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		    
		  end do
		    COMPBASEL=S1
		 END IF

      CASE(2)	

		JX=4;   ELEM_DEC=1
		do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILON_Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		  end do
		    COMPBASEL=S1
		
      CASE(3)

		JX=5;   ELEM_DEC=2
		do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILON_Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		  end do
		    COMPBASEL=S1
	CASE(4)
		JX=6;   ELEM_DEC=3

	      if ((JXX.eq.1))then
		  
		      do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILON_Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
		      END DO
		      CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
		  if (ielem(n,ixx)%mode.eq.0)then
		  
		  COMPBASEL=COMPBASPR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)

		    else
		  
		  do K=1,ELEM_DEC
		   VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
! 		   
		  end do
		    COMPBASEL=S1
		    end if
		  
		  ELSE
			    
		   do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			  NODES_LIST(k,3)=ILON_Z(LXX1,JXX,K)
			VEXT(K,:)=NODES_LIST(k,:)
			
		      END DO
		  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
		  
		  do K=1,ELEM_DEC
		    VEXT(1:4,1)=ELEM_LISTD(k,1:4,1)
		    VEXT(1:4,2)=ELEM_LISTD(k,1:4,2)
		    VEXT(1:4,3)=ELEM_LISTD(k,1:4,3)
		   
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		    
		  end do
		    COMPBASEL=S1
		 END IF

	      CASE(5)
		JX=4;   ELEM_DEC=2
		     if ((JXX.eq.1))then
		     
		     
		   do K=1,JX
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
		    vext(k,1:2)=NODES_LIST(k,1:2)
		    
		    END DO
		    CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
		    if (ielem(n,ixx)%mode.eq.0)then
		  COMPBASEL=COMPBASquad(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		    else
		    
		  do K=1,ELEM_DEC
		    VEXT(1:3,1)=ELEM_LISTD(k,1:3,1)
		    VEXT(1:3,2)=ELEM_LISTD(k,1:3,2)
		    
		    
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTRi(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		    
		  end do
		    COMPBASEL=S1
		    
		    
		    end if
		    
		    
		    
		else
		  do K=1,jx
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
			 
			VEXT(K,:)=NODES_LIST(k,:)
			
		      END DO
		      
		      
		  CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
		  do K=1,ELEM_DEC
		    VEXT(1:3,1)=ELEM_LISTD(k,1:3,1)
		    VEXT(1:3,2)=ELEM_LISTD(k,1:3,2)
		    
		    
		    S1(1:number_of_dog)=S1(1:number_of_dog)+COMPBASTRi(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
		    
		  end do
		    COMPBASEL=S1	
		
		end if

	      CASE(6)

		JX=3;
		
		
		   do K=1,JX
			  NODES_LIST(k,1)=ILON_X(LXX1,JXX,K)
			  NODES_LIST(k,2)=ILON_Y(LXX1,JXX,K)
		    vext(k,1:2)=NODES_LIST(k,1:2)
		    
		    END DO
		    
		  COMPBASEL=COMPBAStri(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)


     END SELECT

END FUNCTION

FUNCTION COMPBASTR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
!> @brief
!> This subroutine computes basis functions for each triangle
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog,IXX,JXX,KXX,LXX1,icompwrt
REAL::VOL
REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(IN)::VEXT
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL::X1,Y1,Z1
integer::lc
REAL,dimension(1:number_of_dog)::INTEG,COMPBASTR


QPOINTS(1:3,1:QP_TETRA)=ZERO
WEQUA3D(1:QP_TETRA)=ZERO

VOL=TETRAVOLUME(N,VEXT)

CALL QUADRATURETETRA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)

INTEG=ZERO
DO Lc=1,QP_TETRA
	  x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc);z1=QPOINTS(3,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC(N,x1,y1,z1,kxx,IXX,number_of_dog,ICOMPWRT)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASTR = INTEG

END FUNCTION


FUNCTION COMPBASHEX(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
!> @brief
!> This subroutine computes basis functions for hexahedral
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog,IXX,JXX,KXX,LXX1,icompwrt
REAL::VOL
REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(IN)::VEXT
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL::X1,Y1,Z1
integer::lc
REAL,dimension(number_of_dog)::INTEG,COMPBASHEX

QPOINTS(1:3,1:QP_HEXA)=zero
WEQUA3D(1:QP_HEXA)=zero

CALL QUADRATUREHEXA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)


 VOL=hexaVOLUME(N,VEXT,QPOINTS,WEQUA3D)

INTEG=zero
DO Lc=1,QP_HEXA
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc);z1=QPOINTS(3,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC(N,x1,y1,z1,kxx,IXX,number_of_dog,ICOMPWRT)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASHEX = INTEG

END FUNCTION

FUNCTION COMPBASPR(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
!> @brief
!> This subroutine computes basis functions for prisms
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog,IXX,JXX,KXX,LXX1,icompwrt
REAL::VOL
REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(IN)::VEXT
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL::X1,Y1,Z1
integer::lc
REAL,dimension(number_of_dog)::INTEG,COMPBASPR
QPOINTS(1:3,1:QP_PRISM)=zero
WEQUA3D(1:QP_PRISM)=zero

CALL QUADRATUREPRISM(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)


 VOL=PRISMVOLUME(N,VEXT,QPOINTS,WEQUA3D)

INTEG=zero
DO Lc=1,QP_PRISM
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc);z1=QPOINTS(3,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC(N,x1,y1,z1,kxx,IXX,number_of_dog,ICOMPWRT)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASPR= INTEG

END FUNCTION

FUNCTION COMPBASQUAD(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
!> @brief
!> This subroutine computes basis functions for quadrilateral
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog,IXX,JXX,KXX,LXX1,icompwrt
REAL::VOL
REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(IN)::VEXT
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL::X1,Y1,Z1
integer::lc
REAL,dimension(number_of_dog)::INTEG,COMPBASQUAD

QPOINTS(1:2,1:QP_QUAD)=zero
WEQUA3D(1:QP_QUAD)=zero

CALL QUADRATUREQUAD(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)


 VOL=QUADVOLUME(N,VEXT,QPOINTS,WEQUA3D)

INTEG=zero
DO Lc=1,qp_quad
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC2D(N,x1,y1,kxx,IXX,number_of_dog,ICOMPWRT)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASQUAD= INTEG

END FUNCTION

FUNCTION COMPBASTRI(N,IXX,JXX,KXX,LXX1,number_of_dog,icompwrt,VEXT)
!> @brief
!> This subroutine computes basis functions for each triangle
IMPLICIT NONE
INTEGER,INTENT(IN)::N,number_of_dog,IXX,JXX,KXX,LXX1,icompwrt
REAL::VOL
REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(IN)::VEXT
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL::X1,Y1,Z1
integer::lc
REAL,dimension(number_of_dog)::INTEG,COMPBAStri

QPOINTS(1:2,1:QP_QUAD)=zero
WEQUA3D(1:QP_QUAD)=zero

CALL QUADRATURETRIANGLE(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)


 VOL=TRIANGLEVOLUME(N,VEXT)

INTEG=zero
DO Lc=1,qp_triangle
	x1=QPOINTS(1,Lc);y1=QPOINTS(2,Lc)
	INTEG(1:number_of_dog)=INTEG(1:number_of_dog)+(BASIS_REC2D(N,x1,y1,kxx,IXX,number_of_dog,ICOMPWRT)*&
WEQUA3D(Lc)*VOL)
END DO
	  COMPBASTRI= INTEG

END FUNCTION


SUBROUTINE INVERT(RFF,INVRFF,IVGT)
!> @brief
!> This subroutine inverts a special matrix
  IMPLICIT NONE
  INTEGER,INTENT(IN)::IVGT
  REAL,allocatable,DIMENSION(:,:),INTENT(IN) ::RFF
  REAL,allocatable,DIMENSION(:,:),INTENT(inOUT)::INVRFF
  INTEGER I,J,K,GT
 
  INVRff(1:IVGT-1,1:IVGT-1) = zero
 GT=IVGT-1
  DO I=GT,1,-1
    INVRFF(i,i) = 1./RFF(i,i)
	do j=i+1,GT
	  INVRFF(i,j) = zero
	  do k= 1,j-1
	   INVRFF(i,j) = INVRFF(i,j) - RFF(k,j)*INVRFF(i,k)
	  enddo
 	  INVRFF(i,j) =INVRFF(i,j) /RFF(j,j)
	enddo
  enddo

 End subroutine INVERT
 


 
 
END MODULE PRESTORE
