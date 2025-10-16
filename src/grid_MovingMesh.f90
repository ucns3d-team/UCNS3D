MODULE TRANSFORM_MovingMesh
	USE MPIINFO
	USE DECLARATION
	USE TRANSFORM

	IMPLICIT NONE

CONTAINS	


SUBROUTINE VOLUME_CALCULATOR_MovingMesh_2D(i, node_position_index)
	!> @brief
	!> This subroutine computes the volume of elements in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::i, node_position_index
	!$ integer::OMP_IN_PARALLEL,OMP_GET_THREAD_NUM
	INTEGER::K,KMAXE,jx,JX2,ELTYPE,ELEM_DEC
	real::DUMV1,DUMV2,dumv3,DUMV5
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:8,1:DIMENSIONA)::NODES_LIST
	REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
	REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D

    ! print *, "I am inside VOLUME_CALCULATOR_MovingMesh_2D for", i, "current volume", IELEM(N,I)%moving_VOLUME(node_position_index)

	ELTYPE=IELEM(N,I)%ISHAPE
	ELEM_DEC=IELEM(N,I)%VDEC
	! IELEM(N,I)%TOTVOLUME=0.0d0

	do K=1,IELEM(N,I)%NONODES
        !print *, "trying to access local_node", IELEM(N,I)%NODES_local(K)
		! NODES_LIST(k,1:2)=LOCAL_NODES(IELEM(N,I)%NODES_local(K))%positions(node_position_index,1:2)
        ! print *, "trying to access local_node", IELEM(N,I)%NODES(K)
		NODES_LIST(k,1:2)=LOCAL_NODES(IELEM(N,I)%NODES(K))%positions(node_position_index,1:2)
        ! print *, k, "node position copied"
		vext(k,1:2)=NODES_LIST(k,1:2)
	END DO

    ! print *, "node positions copied"

	call DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)

	SELECT CASE(ielem(n,i)%ishape)

	  CASE(5)

		CALL QUADRATUREQUAD(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
		DUMV1=QUADVOLUME(N,VEXT,QPOINTS,WEQUA3D)
		DUMV2=0.0d0
		
		do K=1,ELEM_DEC
			VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
			DUMV2=DUMV2+TRIANGLEVOLUME(N,VEXT)
		END DO

		IELEM(N,I)%Moving_VOLUME(node_position_index)=DUMV2
     
	  CASE(6)

		DUMV1=TRIANGLEVOLUME(N,VEXT)
		DUMV2=0.0d0
		do K=1,ELEM_DEC
			VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
			DUMV2=DUMV2+TRIANGLEVOLUME(N,VEXT)
		END DO

		IELEM(N,I)%moving_VOLUME(node_position_index)=DUMV2
     
	END SELECT

    ! print *, "volume = ", IELEM(N,I)%moving_VOLUME(node_position_index)

END SUBROUTINE VOLUME_CALCULATOR_MovingMesh_2D





SUBROUTINE SURFACE_CALCULATOR_MovingMesh_2D(iconsidered, node_position_index)
	!> @brief
	!> This subroutine computes the length of edges of elements in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::iconsidered, node_position_index
	INTEGER::I,K,jx,JX2,ELTYPE,ELEM_DEC,nnd,j
	real::DUMV1, DUMV2, dumv3, DUMV5, DUMR
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT

	i=iconsidered
    DO J=1,IELEM(N,I)%IFCA
        NND=2
        
        DO K=1,NND
            ! VEXT(K,1:DIMS) = LOCAL_NODES(IELEM(N,I)%NODES_local_FACES(J,K))%positions(node_position_index, 1:DIMS)
            VEXT(K,1:DIMS) = LOCAL_NODES(IELEM(N,I)%NODES_FACES(J,K))%positions(node_position_index, 1:DIMS)
        END DO
        
        IELEM(N,I)%SURF(J)=LINEAREA(N,vext)
    END DO

END SUBROUTINE SURFACE_CALCULATOR_MovingMesh_2D





subroutine coordinates_face_inner_MovingMesh_2D(n, i, facex, VEXT, NODES_LIST, node_position_index)
	!> @brief
	!> This subroutine retrieves the nodes of edges of elements in 2D
	IMPLICIT NONE
	integer,intent(in)::n, i, facex, node_position_index
	REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(INOUT)::VEXT
	REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(INOUT)::NODES_LIST
	integer::nnd
	integer::k

	nnd=2
	       
	do K=1,nnd
		NODES_LIST(k,1:2)=local_nodes(IELEM(N,I)%NODES_FACES(facex,K))%positions(node_position_index, 1:2)
		VEXT(K,1:2)=NODES_LIST(k,1:2)
	END DO
	        
end subroutine coordinates_face_inner_MovingMesh_2D





SUBROUTINE EDGE_CALCULATOR_MovingMesh_2D(I, node_position_index)
	!> @brief
	!> This subroutine computes the radius of inscribed sphere or circle
	IMPLICIT NONE
	INTEGER,INTENT(IN)::i, node_position_index
	INTEGER::L,FACEX,N_NODE
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:8,1:DIMENSIONA)::NODES_LIST
	REAL::EDGEL,DIST

	IELEM(N,I)%MINEDGE=(2.0D0*IELEM(N,I)%moving_volume(node_position_index))/(SUM(IELEM(N,I)%SURF(1:IELEM(N,I)%IFCA)))
	
	DO L=1,IELEM(N,I)%IFCA
		FACEX=L
		N_NODE=2
		CALL coordinates_face_inner_MovingMesh_2D(N,I,facex,vext,NODES_LIST, node_position_index)
		
		VEXT(2,1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
		VEXT(1,1)=IELEM(N,I)%XXC;VEXT(1,2)=IELEM(N,I)%YYC; 
		DIST=DISTANCE2(N,VEXT)
		
		IELEM(N,I)%MINEDGE=MIN(DIST,IELEM(N,I)%MINEDGE)
	END DO
	
END SUBROUTINE EDGE_CALCULATOR_MovingMesh_2D





SUBROUTINE FIND_ROT_ANGLES_MovingMesh_2d(N, I, node_position_index)
	!> @brief
	!> This subroutine determines the normal vectors for each edge
	IMPLICIT NONE
	real::Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,DELXYA,DELyzA,DELzxA,DELXYb,DELyzb,DELzxb,DELXYc,DELyzc,DELzxc,nx,ny,nz,ANGLEFACEX,ANGLEFACEY
	REAL::X5,X6,X7,X8,Y5,Y6,Y7,Y8,Z5,Z6,Z7,Z8,XX,YY,ZZ
	REAL::DELXA,DELXB,DELYA,DELYB,DELZA,DELZB
	INTEGER::K,KMAXE,J,kk,kk2,ixf4,IXFV,n_node
	INTEGER,INTENT(IN)::N, I, node_position_index
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT

	DO K=1,IELEM(N,I)%IFCA
		! if (ielem(n,i)%types_faces(k).eq.5)then
			kk2=2;n_node=kk2
		! else
		! 	  kk2=3;n_node=kk2
		! end if
		IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
			IF ((IELEM(N,I)%INEIGHG(K).GT.0).AND.(IELEM(N,I)%IBOUNDS(K).GT.0))THEN 	!PERIODIC NEIGHBOUR

 			    ! XX=IELEM(N,I)%XXC  ;YY=IELEM(N,I)%YYC; !ZZ=IELEM(N,I)%ZZC

				! DO Kk=1,n_node
				!     IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN
				!        vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:2)
				!     ELSE
				! 		vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:2)
				!     END IF

				!     IF(ABS(vext(kk,1)-xx).GT.XPER*oo2)THEN
				!       	vext(kk,1)=vext(kk,1)+(XPER*SIGN(1.0d0,xx-XPER/2.0D0))
				!     end if
				!     IF(ABS(vext(kk,2)-yy).GT.yPER*oo2)THEN
				!       	vext(kk,2)=vext(kk,2)+(yPER*SIGN(1.0d0,yy-yPER/2.0D0))
				!     end if
			    ! end do
			Else
				DO Kk=1,n_node
					IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN
				        ! vext(kk,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,kk))%positions(node_position_index,1:2
                        vext(kk,1:2)=local_nodes(ielem(n,i)%NODES_FACES(k,kk))%positions(node_position_index,1:2)
				    ELSE
						! vext(KK,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
                        vext(KK,1:2)=local_nodes(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
				    END IF
				END DO
			END IF
		ELSE
			DO Kk=1,n_node
				IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN
				    ! vext(kk,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,kk))%positions(node_position_index,1:2)
                    vext(kk,1:2)=local_nodes(ielem(n,i)%NODES_FACES(k,kk))%positions(node_position_index,1:2)
				ELSE
					! vext(KK,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
                    vext(KK,1:2)=local_nodes(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
				END IF
			END DO
		END IF

		CALL ANGLE2D(VEXT,ANGLEFACEX,ANGLEFACEY)

		IELEM(N,I)%FACEANGLEX(K)=anglefacex
		IELEM(N,I)%FACEANGLEY(K)=anglefacey
	END DO

END SUBROUTINE FIND_ROT_ANGLES_MovingMesh_2d





SUBROUTINE FIND_ANGLES_MovingMesh(N, node_position_index)
	INTEGER,INTENT(IN)::N, node_position_index
	INTEGER::I,KMAXE

	KMAXE=XMPIELRANK(N)

	IF (DIMENSIONA.EQ.3)THEN
		!$OMP DO
		do i=1,kmaxe
			! CALL FIND_ROT_ANGLES(N,I)
		end do
		!$OMP END DO
	Else
		!$OMP DO
		do i=1,kmaxe
			CALL FIND_ROT_ANGLES_MovingMesh_2D(N, I, node_position_index)
		end do
		!$OMP END DO
	END IF

END SUBROUTINE FIND_ANGLES_MovingMesh





SUBROUTINE CENTRE_MovingMesh_2D(iconsidered, node_position_index)
	!> @brief
	!> This subroutine computes the cell centres
	IMPLICIT NONE
	INTEGER,INTENT(IN)::iconsidered, node_position_index
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	INTEGER::I
	i=iconsidered

    CALL COMPUTE_CENTRE_MovingMesh_2D(i,CORDS, node_position_index)
    IELEM(N,I)%XXC=CORDS(1)
    IELEM(N,I)%YYC=CORDS(2)

END SUBROUTINE CENTRE_MovingMesh_2D





SUBROUTINE COMPUTE_CENTRE_MovingMesh_2dF(N, I, facex,N_NODE,cords, node_position_index)
	!> @brief
	!> This subroutine retrieves the nodes of the vertices of edges of 2D elements
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N, I, facex, node_position_index
	INTEGER,INTENT(INOUT)::N_NODE
	REAL,dimension(1:dimensiona),INTENT(INOUT)::CORDS
	REAL,DIMENSION(1:8,1:DIMENSIONA)::NODES_LIST
	integer::k

    do K=1,N_NODE
      	NODES_LIST(k,1:2)=local_nodes(IELEM(N,I)%NODES_FACES(facex,K))%positions(node_position_index,1:2)
    END DO
    CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
   
END SUBROUTINE





SUBROUTINE COMPUTE_CENTRE_MovingMesh_2d(Iconsidered, CORDS, node_position_index)
	!> @brief
	!> This subroutine retrieves the nodes of the vertices of 2D elements
	IMPLICIT NONE
	INTEGER,INTENT(IN)::Iconsidered,node_position_index
	integer::k,i,j,n_node
	real,dimension(1:dimensiona),INTENT(INOUT)::cords
	REAL,DIMENSION(1:8,1:dimensiona)::NODES_LIST

	i=iconsidered

    N_NODE=IELEM(N,I)%NONODES
    do K=1,IELEM(N,I)%NONODES
      	! NODES_LIST(k,1:2)=local_nodes(IELEM(N,I)%NODES_local(K))%positions(node_position_index,1:2)
        NODES_LIST(k,1:2)=local_nodes(IELEM(N,I)%NODES(K))%positions(node_position_index,1:2)
    END DO
    CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
   
END SUBROUTINE





SUBROUTINE GEOMETRY_CALC_MovingMesh(n, node_position_index)
	!> @brief
	!> This subroutine computes the volume, surface, centre and min edge for each element
	IMPLICIT NONE
	INTEGER,INTENT(IN)::n, node_position_index
	INTEGER::KMAXE,i
	real::DUMV5

	KMAXE=XMPIELRANK(N)

	if (DIMENSIONA.EQ.3)THEN
		!$OMP DO 
		DO I=1,KMAXE
			! CALL VOLUME_CALCULATOR3(I)
			! call SURFACE_CALCULATOR3(i)
			! CALL CENTRE3D(i)
			! call EDGE_CALCULATOR3d(i)
		END DO
		!$OMP END DO
	ELSE
		!$OMP DO
		DO I=1,KMAXE
			CALL VOLUME_CALCULATOR_MovingMesh_2D(I, node_position_index)
			call SURFACE_CALCULATOR_MovingMesh_2D(I, node_position_index)
			call CENTRE_MovingMesh_2D(I, node_position_index)
			call EDGE_CALCULATOR_MovingMesh_2D(I, node_position_index)
		END DO
		!$OMP END DO 
	END IF

	!$OMP BARRIER 
	!$OMP MASTER
		DUMV5=ZERO
		DO I=1,KMAXE
			DUMV5=DUMV5+IELEM(N,I)%moving_volume(node_position_index)
		END DO
		CALL MPI_ALLREDUCE(DUMV5,Moving_TOTALVOLUME(node_position_index),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
	!$OMP END MASTER
	!$OMP BARRIER 

END SUBROUTINE GEOMETRY_CALC_MovingMesh


END MODULE