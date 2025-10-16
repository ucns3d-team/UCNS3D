MODULE MESHMOVEMENT_module
    USE MPIINFO
    USE DECLARATION
    USE FLUXES
    USE TRANSFORM
    USE TRANSFORM_MovingMesh
    IMPLICIT NONE

CONTAINS

subroutine NodeSwapXY(N)
    implicit NONE
    integer,intent(in)::N
    real::x,y
    integer::node_index

    !$omp do
    do node_index=1,kmaxn 
        ! write(*,*) "n = ", n, " node index = ", node_index
        ! if (inoder4(node_index)%itor.gt.0) then
            write(*,*) "n = ", n, " node index = ", node_index
            x = inoder4(node_index)%cord(1)
            y = inoder4(node_index)%cord(2)

            if ((x.eq.0.0).and.(y.eq.0.0)) then
                write(*,*) "swapping node (0.0, 0.0)\n"
            end if

            inoder4(node_index)%cord(1) = y
            inoder4(node_index)%cord(2) = X
        ! endif
    end do
    !$omp end do    

end subroutine





subroutine MOVE_NODES(time_step, index_from, index_to)
    implicit NONE
    integer,intent(in)::index_from,index_to
    real,intent(in)::time_step
    integer::node_index

    !$omp do
    do node_index=1,kmaxn 
        local_nodes(node_index)%positions(index_to,1:dimensiona) = local_nodes(node_index)%positions(index_from,1:dimensiona) &
                                                                 + time_step * local_nodes(node_index)%velocity(1:dimensiona)
    end do
    !$omp end do

end subroutine





subroutine COPY_BACK_LOCAL_NODES(index_from, index_to)
    implicit NONE
    integer,intent(in)::index_from,index_to
    integer::node_index

    !$omp do
    do node_index=1,kmaxn 
        local_nodes(node_index)%positions(index_to,1:dimensiona) = local_nodes(node_index)%positions(index_from,1:dimensiona)
    end do
    !$omp end do

end subroutine





subroutine find_node_velocities(position_index)
    implicit none
    integer,intent(in)::position_index
    real::x,y
    integer::node_index

    if ((governingequations.eq.3).and.(INITCOND.eq.3)) then
        !$omp do
        do node_index = 1, kmaxn
            x = local_nodes(node_index)%positions(position_index,1)
            y = local_nodes(node_index)%positions(position_index,2)
            local_nodes(node_index)%velocity(1) = 0.5 - Y
            local_nodes(node_index)%velocity(2) = x - 0.5
        end do
        !$omp end do
    else
        print*,"Not supported test case!"
        call abort
    endif

end subroutine





SUBROUTINE QuadraturePoint_velocity(N,IGQRULES,NODE_velocities,QPOINTS2D_velocity)
	!> @brief
	!> This subroutine computes the velocity quadrature points for a line and returns it in QPOINTS2D_velocity(DIM,QP)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,IGQRULES
	REAL,DIMENSION(1:2,1:DIMENSIONA),INTENT(IN)::NODE_velocities
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2),INTENT(INOUT)::QPOINTS2D_velocity
	real,dimension(1:2,1:2)::VVA,VVA1
	REAL,dimension(1)::DETA
	REAL,DIMENSION(1:4)::VVNXI
	real,dimension(1:ALLS)::VVwg
	real,dimension(1:ALLS)::VVR1
	real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
	REAL::R,S,TX,a,b,c,d,e,f,G,H,K
	INTEGER::Kk,J,ii,ij,ik,count1

	SELECT CASE(IGQRULES)
	
	  case(1)

		VVR1(1)=0.0d0

	  case(2)

		a=-0.5773502691896257
		b=0.5773502691896257
		vvnpox(1)=a	;  vvnpox(2)=b		
  
  		count1=0
 		do ii=1,2
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
  	
	  CASE(3)

		a=0.0d0
		b=-0.7745966692414834
		c=0.7745966692414834
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c
		
		count1=0
		do ii=1,3
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
 		
	  CASE(4)
		a=-0.3399810435848563
		b=0.3399810435848563
		c=-0.8611363115940526
		d=0.8611363115940526
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d
		
		count1=0
		do ii=1,4
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
			
	  CASE(5)

		a=0.0d0
		b=-0.5384693101056831
		c=0.5384693101056831
		d=-0.9061798459386640
		e=0.9061798459386640
		vvnpox(1)=a	;  vvnpox(2)=b ;   vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e
		
		count1=0
		do ii=1,5
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
		
	  CASE(6)

		a=0.6612093864662645
		b=-0.6612093864662645
		c=-0.2386191860831969
		d=0.2386191860831969
		e=-0.9324695142031521
		F=0.9324695142031521
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e ;  vvnpox(6)=f
		
		count1=0
		do ii=1,6
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do

	  CASE(7,8,9)

		a=0.0000000000000000
		b=-0.8360311073266358
		C=0.8360311073266358
		D=-0.9681602395076261
		E=0.9681602395076261
		F=-0.3242534234038089
		G=0.3242534234038089
		H=-0.6133714327005904
		K=0.6133714327005904
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e ;  vvnpox(6)=f ;  vvnpox(7)=G ;  vvnpox(8)=H ;  vvnpox(9)=K
	
		count1=0
 		do ii=1,9
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do	
			
	END SELECT

	QPOINTS2D_velocity(:,:)=0.0d0

	do kk=1,IGQRULES
		R=VVR1(kk)		
		QPOINTS2D_velocity(:,kk)=((NODE_velocities(1,1:2)+NODE_velocities(2,1:2))/2.0D0)+(R*((NODE_velocities(2,1:2)-NODE_velocities(1,1:2))/2.0D0))
	END DO

END SUBROUTINE QuadraturePoint_velocity





SUBROUTINE Find_QP_velocities(N)
    !> @brief
    !> Subroutine for computing and storing velocities of gaussian quadrature points at cell interfaces
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX
    INTEGER::ICONSIDERED,FACEX,POINTX
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D_velocity
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    ! REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
    REAL,DIMENSION(1:8,1:DIMENSIONA)::Node_velocities
    REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
    
    KMAXE=XMPIELRANK(N)
    
    if (dimensiona.eq.3)then
        ! DO I=1,KMAXE
        !     IF (IELEM(N,I)%ISHAPE.EQ.2)THEN
        !         ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !         IF (SRFG.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !         END IF
        !         IF (MRF.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_ORIGIN(1:3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_VELOCITY(1:3))
        !         END IF
        !     ELSE
        !         ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
        !         IF (SRFG.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_QUAD,3))
        !         END IF
        !         IF (MRF.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_QUAD,3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_ORIGIN(1:3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_VELOCITY(1:3))
        !         END IF
        !     END IF
        !     ICONSIDERED=I
        !     DO L=1,IELEM(N,I)%IFCA
        !         IDUMMY=0
        !         if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
        !             IF (IELEM(N,I)%IBOUNDS(l).GT.0) THEN	!CHECK FOR BOUNDARIES
        !                 if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50)) then	!PERIODIC IN OTHER CPU
        !                     IDUMMY=1
        !                 END IF
        !             END IF	
        !             if (ielem(n,i)%types_faces(L).eq.5) then
        !                 iqp=qp_quad
        !                 NND=4
        !                 IF (IDUMMY.EQ.0)THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 ELSE
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     ! do K=1,nnd
        !                     !     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     ! END DO
        !                 END IF
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3
        !                 IF (IDUMMY.EQ.0)THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 ELSE
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     ! do K=1,nnd
        !                     !     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     ! END DO
        !                 END IF
                            
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         else
        !             if (ielem(n,i)%types_faces(L).eq.5)then
        !                 iqp=qp_quad
        !                 NND=4
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO 
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3 
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO  
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         end if
                
        !         do NGP=1,iqp			!for gqp
        !             ! ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !             IF (SRFG.EQ.1) THEN
        !                 ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !                 POX(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)-SRF_ORIGIN(1:3)
        !                 POY(1:3)=SRF_VELOCITY(1:3)
        !                 ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)=VECT_FUNCTION(POX,POY)
        !             END IF
                        
        !             IF (MRF.EQ.1) THEN
        !                 ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !                 POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
        !                 POY(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)
        !                 ICONSIDERED=I
        !                 FACEX=L
        !                 POINTX=NGP
        !                 CALL MRFSWITCH(N,ICONSIDERED,FACEX,POINTX,pox,poy)
                        
        !             END IF 
        !         END DO	!NGP
        !     END DO
    
        !     DO L=1,IELEM(N,I)%IFCA
        !         IDUMMY=0
        !         if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
        !             IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
        !                 if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
        !                     IDUMMY=1
        !                 END IF
        !             END IF	
        !             if (ielem(n,i)%types_faces(L).eq.5) then
        !                 iqp=qp_quad
        !                 NND=4
        !                 IF (IDUMMY.EQ.0) THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO	!NGP
        !                 ELSE
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 END IF
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3
        !                 IF (IDUMMY.EQ.0)THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        
        !                 ELSE ! 2 dimensions (Michael's edit: ?????)
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 END IF
                                                
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         else
        !             if (ielem(n,i)%types_faces(L).eq.5)then
        !                 iqp=qp_quad
        !                 NND=4
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO 
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3 
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO  
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         end if
            
        !         do NGP=1,iqp			!for gqp
        !             ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !         END DO	!NGP
        !     END DO
        ! END DO
          
    else ! 2 dimensions
        
        DO I=1,KMAXE
            
            ICONSIDERED=I
        
            DO L=1,IELEM(N,I)%IFCA
                IDUMMY=0
            
                ! if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                !     IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
                !         if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
                !             IDUMMY=1
                !         END IF
                !     END IF
            
                !     IQP=QP_LINE
                !     NND=2
                !     IF (IDUMMY.EQ.0)THEN
                !         DO K=1,NND
                !             ! VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                !             Node_velocities(k,1:2) = local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(1:2)
                !             !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                !                 ! VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                !             !END IF
                !         END DO
                !     ELSE
                !         ! facex=l;
                !         ! CALL coordinates_face_PERIOD2D1(n,iconsidered,facex,VEXT,NODES_LIST)
                !         ! DO K=1,NND
                !         !     !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                !         !         VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                !         !     !END IF
                !         ! END DO
                !     END IF
                !     CALL QUADRATURELINE(N,IGQRULES,Node_velocities,QPOINTS2D_velocity,WEQUA2D)
                ! ELSE
                    IQP=QP_LINE
                    NND=2
                    DO K=1,NND
                        Node_velocities(k,1:2) = local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(1:2)
                        ! print*, "velocity", local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(1), local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(2), "expected", (0.5 - local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(1,2)), (local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(1,1) - 0.5)
                        ! !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                        !     VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                        ! !END IF
                    END DO
                    CALL QUADRATURELINE(N,IGQRULES,Node_velocities,QPOINTS2D_velocity,WEQUA2D)
                ! END IF
            
                DO NGP=1,iqp !for gqp
                    ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2)=QPOINTS2D_velocity(1:2,NGP) ! Storing surface quadrature points
                END DO !NGP
            END DO
        END DO
    END IF
END SUBROUTINE Find_QP_velocities
  




SUBROUTINE CALCULATE_FLUXESHI_MovingMesh_2D(N)
    !> @brief
    !> This subroutine computes the fluxes for linear-advection equation in 2D
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    REAL::GODFLUX2,sum_detect,lamxl,lamyl
    INTEGER::I,L,K,NGP,KMAXE,IQP, NEIGHBOR_INDEX, NEIGHBOR_FACE_INDEX
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_TEMP,WEIGHTS_DG !Quadrature weights for interfaces
    REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    REAL::ANGLE1,ANGLE2,NX,NY,NZ,NORMALVECT, normal_qp_velocity, net_normal_velocity
    INTEGER::facex,POINTX,ICONSIDERED
    REAL,DIMENSION(1:NOF_VARIABLES)::CLEFT,CRIGHT,HLLCFLUX,RHLLCFLUX
    REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG
    real,dimension(1:dimensiona)::vleft,vright, qp_velocity
  
    IF (DG.EQ.1)THEN
        print *, "Discontinous Galerkin not supported in moving mesh mode"
        call abort
        allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
    END IF
    
    KMAXE = XMPIELRANK(N)
    
    CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
    WEIGHTS_TEMP(1:QP_LINE_N) = WEQUA2D(1:QP_LINE_N)
      
    !$OMP DO
    DO I=1,KMAXE
    
        if (initcond.eq.3)then
            lamxl =-ielem(n,i)%yyc + 0.5d0
            lamyl = ielem(n,i)%xxc - 0.5d0
        else
            lamxl=lamx
            lamyl=lamy
        end if
        
        if (dg.eq.1)then
            RHS(I)%VALDG = ZERO
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
        else
            RHS(I)%VAL=ZERO
        end if
        
        ICONSIDERED=I
        
        IF (DG.EQ.1) THEN
            DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
        END IF
          
        IF (IELEM(N,I)%INTERIOR.EQ.0) THEN ! Element is interior
  
            DO L=1,IELEM(N,I)%IFCA
                GODFLUX2=ZERO
                NX=IELEM(N,I)%FACEANGLEX(L)
                NY=IELEM(N,I)%FACEANGLEY(L)
                facex=l
                
                NORMALVECT=(NX*LAMXl)+(NY*LAMYl)

                IQP=QP_LINE_N
                
                NEIGHBOR_INDEX = IELEM(N,I)%INEIGH(L)
                NEIGHBOR_FACE_INDEX = IELEM(N,I)%INEIGHN(L)
                  
                DO NGP=1,IQP
                    POINTX=NGP
  
                    IF (DG.EQ.1) THEN
                        CLEFT  = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                    ELSE !FV
                        ! CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
                        ! CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                        CLEFT(1)=U_C(I)%VAL(1,1)
                        CRIGHT(1)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1)
                    END IF
                    if (CLEFT(1).ne.CLEFT(1)) print*, "CLEFT(1) NaN"
                    if (CRIGHT(1).ne.CRIGHT(1)) print*, "CRIGHT(1) NaN"
                    vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2)
                    vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(IELEM(N,I)%INEIGHN(L),NGP,1:2)
                    if ((vleft(1) .ne.vleft(1) ).or.(vleft(2) .ne.vleft(2) )) print*, "vleft NaN"
                    if ((vright(1).ne.vright(1)).or.(vright(2).ne.vright(2))) print*, "vright NaN"
                    qp_velocity = 0.5*(vleft + vright)
                    normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                    net_normal_velocity = NORMALVECT - normal_qp_velocity
                    ! print*,"net normal velocity", net_normal_velocity
                    CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,net_normal_velocity,HLLCFLUX)
                    if (HLLCFLUX(1).ne.HLLCFLUX(1)) print*, "NaN HLLC flux"
                    IF (DG.EQ.1) THEN
                        ! Riemann flux at interface quadrature points times basis
                        RHLLCFLUX(1)=HLLCFLUX(1)
                        DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
                    ELSE !FV
                        GODFLUX2 = GODFLUX2 + (HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                    END IF
                END DO
                  
                ! if (RHS(I)%VAL(1).ne.RHS(I)%VAL(1)) print*, "NaN flux"
                IF (DG /= 1) RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2
  
            END DO
  
        ELSE IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
  
            DO L=1,IELEM(N,I)%IFCA
                FACEX = L
                NX=IELEM(N,I)%FACEANGLEX(L)
                NY=IELEM(N,I)%FACEANGLEY(L)
                NORMALVECT=(NX*LAMXl)+(NY*LAMYl)
                IQP=QP_LINE_N
                
                GODFLUX2=ZERO
                DO NGP=1,IQP
                    POINTX = NGP

                    IF (DG == 1) THEN
                        CLEFT = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                    ELSE
                        ! CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
                        CLEFT(1)=U_C(I)%VAL(1,1)
                    END IF
                    
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                IF (DG == 1) THEN
                                    CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE !FV
                                    CRIGHT(1) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                                    ! CRIGHT(1)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1)
                                END IF
                            ELSE !NOT PERIODIC ONES IN MY CPU
                                CRIGHT(1:nof_variables)=CLEFT(1:nof_variables)
                            END IF

                            ! if (CLEFT(1).ne.CLEFT(1)) print*, "CLEFT(1) NaN at an interface"
                            ! if (CRIGHT(1).ne.CRIGHT(1)) print*, "CRIGHT(1) NaN at an interface"
                            vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2)
                            vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(IELEM(N,I)%INEIGHN(L),NGP,1:2)
                            ! if ((vleft(1) .ne.vleft(1) ).or.(vleft(2) .ne.vleft(2) )) print*, "vleft NaN at an interface"
                            ! if ((vright(1).ne.vright(1)).or.(vright(2).ne.vright(2))) print*, "vright NaN at an interface"
                            qp_velocity = 0.5*(vleft+vright)
                            normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                            net_normal_velocity = NORMALVECT - normal_qp_velocity
                            ! print*,"net normal velocity", net_normal_velocity
                            CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,net_normal_velocity,HLLCFLUX)
                            ! if (HLLCFLUX(1).ne.HLLCFLUX(1)) print*, "NaN HLLC flux at an interface"
                            IF (DG.EQ.1) THEN
                                !Riemann flux at interface quadrature points times basis
                                RHLLCFLUX(1)=HLLCFLUX(1)
                                DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX) 
                            ELSE !FV
                                ! GODFLUX2 = GODFLUX2 + (HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                            END IF
                        ELSE
                            IF (DG == 1) THEN
                                CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                            ELSE !FV
                                ! CRIGHT(1) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                                CRIGHT(1)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1)
                            END IF

                            ! if (CLEFT(1).ne.CLEFT(1)) print*, "CLEFT(1) NaN at an interface"
                            ! if (CRIGHT(1).ne.CRIGHT(1)) print*, "CRIGHT(1) NaN at an interface"
                            vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2)
                            vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(IELEM(N,I)%INEIGHN(L),NGP,1:2)
                            ! if ((vleft(1) .ne.vleft(1) ).or.(vleft(2) .ne.vleft(2) )) print*, "vleft NaN at an interface"
                            ! if ((vright(1).ne.vright(1)).or.(vright(2).ne.vright(2))) print*, "vright NaN at an interface"
                            qp_velocity = 0.5*(vleft+vright)
                            normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                            net_normal_velocity = NORMALVECT - normal_qp_velocity
                            ! print*,"net normal velocity", net_normal_velocity
                            CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,net_normal_velocity,HLLCFLUX)
                            ! if (HLLCFLUX(1).ne.HLLCFLUX(1)) print*, "NaN HLLC flux at an interface"
                            IF (DG.EQ.1) THEN
                                !Riemann flux at interface quadrature points times basis
                                RHLLCFLUX(1)=HLLCFLUX(1)
                                DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX) 
                            ELSE !FV
                                GODFLUX2 = GODFLUX2 + (HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                            END IF
                        END IF
                        
                    ELSE !IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                            
                                IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                else
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                end if
                            END IF
                        ELSE
                            IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                            ELSE
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                            end if
                        END IF

                        ! vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:dimensiona)
                        ! ! vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(L,NGP,1:dimensiona)
                        ! ! qp_velocity = 0.5*(vleft+vright)
                        ! qp_velocity = vleft
                        ! normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                        ! CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT - normal_qp_velocity,HLLCFLUX)
    
                        ! IF (DG.EQ.1) THEN
                        !     !Riemann flux at interface quadrature points times basis
                        !     RHLLCFLUX(1)=HLLCFLUX(1)
                        !     DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX) 
                        ! ELSE !FV
                        !     ! GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                        ! END IF
                    END IF
                      
                    ! vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:dimensiona)
                    ! ! vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(L,NGP,1:dimensiona)
                    ! ! qp_velocity = 0.5*(vleft+vright)
                    ! qp_velocity = vleft
                    ! normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                    ! CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT - normal_qp_velocity,HLLCFLUX)

                    ! IF (DG.EQ.1) THEN
                    !     !Riemann flux at interface quadrature points times basis
                    !     RHLLCFLUX(1)=HLLCFLUX(1)
                    !     DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX) 
                    ! ELSE !FV
                    !     ! GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                    ! END IF
    
                END DO

                ! if (RHS(I)%VAL(1).ne.RHS(I)%VAL(1)) print*, "NaN flux"
                IF (DG /= 1) RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2

            END DO
  
        END IF
          
        IF (DG == 1) DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
        IF (DG == 1) RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS
          
    END DO
    !$OMP END DO 
  
    IF (DG.EQ.1)THEN
        deallocate(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
    END IF
  
END SUBROUTINE CALCULATE_FLUXESHI_MovingMesh_2D





  





end module