MODULE TRANSFORM_MovingMesh
	USE MPIINFO
	USE OMP_LIB
	USE DECLARATION
	USE TRANSFORM

	IMPLICIT NONE

CONTAINS	

function cross_product_2D(v1, v2)
	implicit none
	real,dimension(1:2)::v1, v2
	real::cross_product_2D

	cross_product_2D = (v1(1)*v2(2)) - (v1(2)*v2(1))

end function





! SUBROUTINE VOLUME_CALCULATOR_MovingMesh_2D(i, node_position_index)
! 	!> @brief
! 	!> This subroutine computes the volume of elements in 2D
! 	IMPLICIT NONE
! 	INTEGER,INTENT(IN)::i, node_position_index
! 	!$ integer::OMP_IN_PARALLEL,OMP_GET_THREAD_NUM
! 	INTEGER::K,KMAXE,jx,JX2,ELTYPE,ELEM_DEC
! 	real::DUMV1,DUMV2,dumv3,DUMV5
! 	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
! 	REAL,DIMENSION(1:8,1:DIMENSIONA)::NODES_LIST
! 	REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
! 	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
! 	REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D

!     ! print *, "I am inside VOLUME_CALCULATOR_MovingMesh_2D for", i, "current volume", IELEM(N,I)%moving_VOLUME(node_position_index)

! 	ELTYPE=IELEM(N,I)%ISHAPE
! 	ELEM_DEC=IELEM(N,I)%VDEC
! 	! IELEM(N,I)%TOTVOLUME=0.0d0

! 	do K=1,IELEM(N,I)%NONODES
!         !print *, "trying to access local_node", IELEM(N,I)%NODES_local(K)
! 		! NODES_LIST(k,1:2)=LOCAL_NODES(IELEM(N,I)%NODES_local(K))%positions(node_position_index,1:2)
!         ! print *, "trying to access local_node", IELEM(N,I)%NODES(K)
! 		NODES_LIST(k,1:2)=LOCAL_NODES(IELEM(N,I)%NODES(K))%positions(node_position_index,1:2)
!         ! print *, k, "node position copied"
! 		vext(k,1:2)=NODES_LIST(k,1:2)
! 	END DO

!     ! print *, "node positions copied"

! 	call DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)

! 	SELECT CASE(ielem(n,i)%ishape)

! 	  CASE(5)

! 		CALL QUADRATUREQUAD(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
! 		DUMV1=QUADVOLUME(N,VEXT,QPOINTS,WEQUA3D)
! 		DUMV2=0.0d0
		
! 		do K=1,ELEM_DEC
! 			VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
! 			DUMV2=DUMV2+TRIANGLEVOLUME(N,VEXT)
! 		END DO

! 		IELEM(N,I)%Moving_VOLUME(node_position_index)=DUMV2
     
! 	  CASE(6)

! 		DUMV1=TRIANGLEVOLUME(N,VEXT)
! 		DUMV2=0.0d0
! 		do K=1,ELEM_DEC
! 			VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
! 			DUMV2=DUMV2+TRIANGLEVOLUME(N,VEXT)
! 		END DO

! 		IELEM(N,I)%moving_VOLUME(node_position_index)=DUMV2
     
! 	END SELECT

!     ! print *, "volume = ", IELEM(N,I)%moving_VOLUME(node_position_index)

! END SUBROUTINE VOLUME_CALCULATOR_MovingMesh_2D





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
        DO K=1, NND
            ! VEXT(K,1:DIMS) = LOCAL_NODES(IELEM(N,I)%NODES_local_FACES(J,K))%positions(node_position_index, 1:DIMS)
            VEXT(K,1:DIMS) = LOCAL_NODES(IELEM(N,I)%NODES_FACES(J,K))%positions(node_position_index, 1:DIMS)
        END DO
        
        IELEM(N,I)%SURF(J) = LINEAREA(N,vext)
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
	do K=1, nnd
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

	IELEM(N,I)%MINEDGE = (2.0D0*IELEM(N,I)%moving_volume(node_position_index))/(SUM(IELEM(N,I)%SURF(1:IELEM(N,I)%IFCA)))
	
	DO L=1, IELEM(N,I)%IFCA
		FACEX=L
		N_NODE=2
		CALL coordinates_face_inner_MovingMesh_2D(N,I,facex,vext,NODES_LIST, node_position_index)
		
		VEXT(2,1:2)= CORDINATES2(N,NODES_LIST,N_NODE)
		VEXT(1,1)  = IELEM(N,I)%XXC
		VEXT(1,2)  = IELEM(N,I)%YYC; 
		DIST = DISTANCE2(N,VEXT)
		
		IELEM(N,I)%MINEDGE = MIN(DIST, IELEM(N,I)%MINEDGE)
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

 			    XX = IELEM(N,I)%XXC
				YY = IELEM(N,I)%YYC
				! ZZ = IELEM(N,I)%ZZC

				DO Kk=1,n_node
				    IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN
				       !vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:2)
						vext(kk,1:2) = local_nodes(ielem(n,i)%NODES_FACES(k,kk))%positions(node_position_index, 1:2)
				    ELSE
						!vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:2)
						vext(kk,1:2) = local_nodes(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%positions(node_position_index, 1:2)
				    END IF

				    IF(ABS(vext(kk,1)-xx).GT.XPER*oo2)THEN
				      	vext(kk,1) = vext(kk,1)+(XPER*SIGN(1.0d0,xx-XPER/2.0D0))
				    end if
				    IF(ABS(vext(kk,2)-yy).GT.yPER*oo2)THEN
				      	vext(kk,2) = vext(kk,2)+(yPER*SIGN(1.0d0,yy-yPER/2.0D0))
				    end if
			    end do
			Else
				DO Kk=1, n_node
					IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN
				        ! vext(kk,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,kk))%positions(node_position_index,1:2
                        vext(kk,1:2) = local_nodes(ielem(n,i)%NODES_FACES(k,kk))%positions(node_position_index,1:2)
				    ELSE
						! vext(KK,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
                        vext(KK,1:2) = local_nodes(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
				    END IF
				END DO
			END IF
		ELSE
			DO Kk=1, n_node
				IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN
				    ! vext(kk,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,kk))%positions(node_position_index,1:2)
                    vext(kk,1:2) = local_nodes(ielem(n,i)%NODES_FACES(k,kk))%positions(node_position_index,1:2)
				ELSE
					! vext(KK,1:2)=local_nodes(ielem(n,i)%NODES_local_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
                    vext(KK,1:2) = local_nodes(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%positions(node_position_index,1:2)
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
		do i=1, kmaxe
			print*,"IND_ROT_ANGLES_MovingMesh_3D not implemented yet!"
			call abort()
			! CALL FIND_ROT_ANGLES_MovingMesh_3D(N, I, node_position_index)
		end do
		!$OMP END DO
	Else
		!$OMP DO
		do i=1, kmaxe
			CALL FIND_ROT_ANGLES_MovingMesh_2D(N, I, node_position_index)
		end do
		!$OMP END DO
	END IF

END SUBROUTINE FIND_ANGLES_MovingMesh





! SUBROUTINE CENTRE_MovingMesh_2D(iconsidered, node_position_index)
! 	!> @brief
! 	!> This subroutine computes the cell centres
! 	IMPLICIT NONE
! 	INTEGER,INTENT(IN)::iconsidered, node_position_index
! 	REAL,DIMENSION(1:DIMENSIONA)::CORDS
! 	INTEGER::I
! 	i=iconsidered

!     CALL COMPUTE_CENTRE_MovingMesh_2D(i,CORDS, node_position_index)
!     IELEM(N,I)%XXC=CORDS(1)
!     IELEM(N,I)%YYC=CORDS(2)

! END SUBROUTINE CENTRE_MovingMesh_2D





SUBROUTINE CENTRE_MovingMesh_2D(cell_index, node_position_index)
	!> @brief
	!> This subroutine computes the cell centres
	IMPLICIT NONE
	INTEGER,INTENT(IN)::cell_index, node_position_index
	REAL,DIMENSION(1:DIMENSIONA)::CORDS, temp_cords, trinagle_centre
	real,dimension(1:dimensiona)::v1, v2
	INTEGER::I, j
	integer::num_nodes, node_index, node_index_1, node_index_2, node_index_3
	real::area, area_sum
	
	cords(:) = zero

	num_nodes = IELEM(N,cell_index)%nonodes
	do i = 1, num_nodes
		node_index = ielem(N, cell_index)%nodes_counterclockwise(i)
		cords(:) = cords(:) + (local_nodes(node_index)%positions(node_position_index, 1:dimensiona)/real(num_nodes))
	end do

	if (num_nodes.gt.3) then
		temp_cords(:) = cords(:)
		cords(:) = zero
		area_sum = zero
		do i = 1, num_nodes
			j = i+1
			if (j.gt.num_nodes) then
				j = j - num_nodes
			end if

			node_index_1 = ielem(N, cell_index)%nodes_counterclockwise(i)
			node_index_2 = ielem(N, cell_index)%nodes_counterclockwise(j)

			v1(:) = local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona)
			v2(:) = temp_cords(1:dimensiona) - local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona)
			
			area = cross_product_2D(v1, v2)*0.5
			area_sum = area_sum + area

			trinagle_centre(:) = (local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona) &
							   + local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) &
							   + temp_cords(1:dimensiona)) &
							   / 3.0
			
			cords(:) = cords(:) + (trinagle_centre(:)*area)
		end do

		cords(:) = cords(:)/area_sum

	end if

    IELEM(N, cell_index)%XXC=CORDS(1)
    IELEM(N, cell_index)%YYC=CORDS(2)

END SUBROUTINE CENTRE_MovingMesh_2D





SUBROUTINE VOLUME_CALCULATOR_MovingMesh_2D(cell_index, node_position_index)
  !> @brief
  !> This subroutine computes the cell centres
	IMPLICIT NONE
	INTEGER,INTENT(IN)::cell_index, node_position_index
	REAL,DIMENSION(1:DIMENSIONA)::temp_cords
	real,dimension(1:dimensiona)::v1, v2
	INTEGER::I, j
	integer::num_nodes, node_index, node_index_1, node_index_2, node_index_3
	real::area, area_sum
	
	num_nodes = IELEM(N,cell_index)%nonodes

	if (num_nodes.gt.3) then
		temp_cords(:) = zero
		do i = 1, num_nodes
			node_index = ielem(N, cell_index)%nodes_counterclockwise(i)
			temp_cords(:) = temp_cords(:) + (local_nodes(node_index)%positions(node_position_index, 1:dimensiona)/real(num_nodes))
		end do

		area_sum = zero
		do i = 1, num_nodes
			j = i+1
			if (j.gt.num_nodes) then
				j = j - num_nodes
			end if

			node_index_1 = ielem(N, cell_index)%nodes_counterclockwise(i)
			node_index_2 = ielem(N, cell_index)%nodes_counterclockwise(j)

			v1(:) = local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona)
			v2(:) = temp_cords(1:dimensiona) - local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona)
			
			area = cross_product_2D(v1, v2)*0.5
			area_sum = area_sum + area
		end do

	else
		node_index_1 = ielem(N, cell_index)%nodes_counterclockwise(1)
		node_index_2 = ielem(N, cell_index)%nodes_counterclockwise(2)
		node_index_3 = ielem(N, cell_index)%nodes_counterclockwise(3)
			
		v1(:) = local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona)
		v2(:) = local_nodes(node_index_3)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona)

		area_sum = cross_product_2D(v1, v2)*0.5

	end if

	IELEM(N, cell_index)%moving_VOLUME(node_position_index) = area_sum

END SUBROUTINE VOLUME_CALCULATOR_MovingMesh_2D





SUBROUTINE CENTREandVOLUME_MovingMesh_2D(cell_index, node_position_index)
	!> @brief
	!> This subroutine computes the cell centres
	IMPLICIT NONE
	INTEGER,INTENT(IN)::cell_index, node_position_index
	REAL,DIMENSION(1:DIMENSIONA)::CORDS, temp_cords, trinagle_centre
	real,dimension(1:dimensiona)::v1, v2
	INTEGER::I, j
	integer::num_nodes, node_index, node_index_1, node_index_2, node_index_3
	real::area, area_sum
	
	cords(:) = zero

	num_nodes = IELEM(N,cell_index)%nonodes
	do i = 1, num_nodes
		node_index = ielem(N, cell_index)%nodes_counterclockwise(i)
		cords(:) = cords(:) + (local_nodes(node_index)%positions(node_position_index, 1:dimensiona)/real(num_nodes))
	end do

	if (num_nodes.gt.3) then
		temp_cords(:) = cords(:)
		cords(:) = zero
		area_sum = zero
		do i = 1, num_nodes
			j = i+1
			if (j.gt.num_nodes) then
				j = j - num_nodes
			end if

			node_index_1 = ielem(N, cell_index)%nodes_counterclockwise(i)
			node_index_2 = ielem(N, cell_index)%nodes_counterclockwise(j)

			v1(:) = local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona)
			v2(:) = temp_cords(1:dimensiona) - local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona)
			
			area = cross_product_2D(v1, v2)*0.5
			area_sum = area_sum + area

			trinagle_centre(:) = (local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona) &
							   + local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) &
							   + temp_cords(1:dimensiona)) &
							   / 3.0
			
			cords(:) = cords(:) + (trinagle_centre(:)*area)
		end do

		cords(:) = cords(:)/area_sum
	else
		node_index_1 = ielem(N, cell_index)%nodes_counterclockwise(1)
		node_index_2 = ielem(N, cell_index)%nodes_counterclockwise(2)
		node_index_3 = ielem(N, cell_index)%nodes_counterclockwise(3)
			
		v1(:) = local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_1)%positions(node_position_index, 1:dimensiona)
		v2(:) = local_nodes(node_index_3)%positions(node_position_index, 1:dimensiona) - local_nodes(node_index_2)%positions(node_position_index, 1:dimensiona)

		area_sum = cross_product_2D(v1, v2)*0.5

	end if

    IELEM(N, cell_index)%XXC=CORDS(1)
    IELEM(N, cell_index)%YYC=CORDS(2)

	IELEM(N, cell_index)%moving_VOLUME(node_position_index) = area_sum

END SUBROUTINE CENTREandVOLUME_MovingMesh_2D



SUBROUTINE QUADRATURELINE_WeightsOnly(N,IGQRULES,WEQUA2D)
  !> @brief
  !> This subroutine computes the quadrature points for a line and returns it in QPOINTS2D(DIM,QP)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,IGQRULES
	REAL,DIMENSION(1:NUMBEROFPOINTS2),INTENT(INOUT)::WEQUA2D
	! real,dimension(1:2,1:2)::VVA,VVA1
	! REAL,dimension(1)::DETA
	! REAL,DIMENSION(1:4)::VVNXI
	real,dimension(1:ALLS)::VVwg
	! real,dimension(1:ALLS)::VVR1,VVR2,VVR3
	real,dimension(1:igqrules)::vvwpox ! ,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
	REAL::R,S,TX,a,b,c,d,e,f,G,H,K
	REAL::a1,b1,c1,d1,e1,f1,G1,H1,K1
	INTEGER::Kk,J,ii,ij,ik,count1

	WEQUA2D=0.0d0
	! QPOINTS2D=0.0d0

	SELECT CASE(IGQRULES)
	  
	  case(1)
		vvwg(1) = 2.0d0
		! VVR1(1)=0.0d0
		! VVR2(1)=0.0d0	
  
	  case(2)
		a=-0.5773502691896257
		b=0.5773502691896257
		a1=1.0d0
		b1=1.0d0
		
		! vvnpox(1)=a	;  vvnpox(2)=b	
		vvwpox(1)=a1 ; vvwpox(2)=b1	
	
		count1=0
		do ii=1,2
			count1=count1+1
			! VVR1(count1)=vvnpox(ii)
			vvwg(count1)=vvwpox(ii)
		end do
		
	  CASE(3)
		a=0.0d0
		b=-0.7745966692414834
		c=0.7745966692414834
		a1=0.8888888888888888
		b1=0.5555555555555556
		c1=0.5555555555555556

		! vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c
		vvwpox(1)=a1 ; vvwpox(2)=b1	; vvwpox(3)=c1
		
		count1=0
		do ii=1,3
			count1=count1+1
			! VVR1(count1)=vvnpox(ii)
			vvwg(count1)=vvwpox(ii)
		end do
		   
	  CASE(4)
		a=-0.3399810435848563
		b=0.3399810435848563
		c=-0.8611363115940526
		d=0.8611363115940526
		a1=0.6521451548625461
		b1=0.6521451548625461
		c1=0.3478548451374538
		d1=0.3478548451374538

		! vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d
		vvwpox(1)=a1 ; vvwpox(2)=b1	; vvwpox(3)=c1 ; vvwpox(4)=d1
		
		count1=0
		do ii=1,4
			count1=count1+1
			! VVR1(count1)=vvnpox(ii)
			vvwg(count1)=vvwpox(ii)
		end do
			  
	  CASE(5)
		a=0.0d0
		b=-0.5384693101056831
		c=0.5384693101056831
		d=-0.9061798459386640
		e=0.9061798459386640
		a1=0.5688888888888889
		b1=0.4786286704993665
		c1=0.4786286704993665
		d1=0.2369268850561891
		E1=0.2369268850561891

		! vvnpox(1)=a	;  vvnpox(2)=b ;   vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e
		vvwpox(1)=a1 ; vvwpox(2)=b1	;  vvwpox(3)=c1 ; vvwpox(4)=d1 ; vvwpox(5)=e1
		
		count1=0
		do ii=1,5
			count1=count1+1
			! VVR1(count1)=vvnpox(ii)
			vvwg(count1)=(vvwpox(ii))
		end do
	  
	  CASE(6)
		a=0.6612093864662645
		b=-0.6612093864662645
		c=-0.2386191860831969
		d=0.2386191860831969
		e=-0.9324695142031521
		F=0.9324695142031521
		a1=0.3607615730481386
		b1=0.3607615730481386
		c1=0.4679139345726910
		d1=0.4679139345726910
		e1=0.1713244923791704
		F1=0.1713244923791704

		! vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e ;  vvnpox(6)=f
		vvwpox(1)=a1 ; vvwpox(2)=b1	; vvwpox(3)=c1 ; vvwpox(4)=d1 ; vvwpox(5)=e1 ; vvwpox(6)=f1
		
		count1=0
		do ii=1,6
			count1=count1+1
			! VVR1(count1)=vvnpox(ii)
			vvwg(count1)=vvwpox(ii)
		end do
  
	  CASE(7,8,9)
		A1=0.3302393550012598	;a=0.0000000000000000
		B1=0.1806481606948574	;b=-0.8360311073266358
		C1=0.1806481606948574	;C=0.8360311073266358
		D1=0.0812743883615744	;D=-0.9681602395076261
		E1=0.0812743883615744	;E=0.9681602395076261
		F1=0.3123470770400029	;F=-0.3242534234038089
		G1=0.3123470770400029	;G=0.3242534234038089
		H1=0.2606106964029354	;H=-0.6133714327005904
		K1=0.2606106964029354	;K=0.6133714327005904
	
		! vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e ;  vvnpox(6)=f ;  vvnpox(7)=G ;  vvnpox(8)=H ;  vvnpox(9)=K
		vvwpox(1)=a1 ; vvwpox(2)=b1	; vvwpox(3)=c1 ; vvwpox(4)=d1 ; vvwpox(5)=e1 ; vvwpox(6)=f1 ; vvwpox(7)=G1 ; vvwpox(8)=H1 ; vvwpox(9)=K1
	
		count1=0
		do ii=1,9
			count1=count1+1
			! VVR1(count1)=vvnpox(ii)
			vvwg(count1)=vvwpox(ii)
		end do	
			  
	END SELECT
  
	vvwg(:)=vvwg(:)*0.5d0
  
	do kk=1,qp_LINE
		WEQUA2D(kk)=vvwg(kk)	
	END DO
  
END SUBROUTINE QUADRATURELINE_WeightsOnly





function VolumeChangeFromVelocity_2D(N, cell_index, d_t, qp_weights)
	implicit none
	integer,intent(in)::N, cell_index
	real,intent(in)::d_t
	real,dimension(1:qp_line_n),intent(in)::qp_weights
	real::VolumeChangeFromVelocity_2D

	integer::face_index, qp_index, num_qp
	real,dimension(1:dimensiona)::qp_velocity
	real::nx, ny, qp_normal_velocity

	num_qp = qp_line_n
	VolumeChangeFromVelocity_2D = zero

	do face_index = 1, IELEM(N,cell_index)%IFCA
		nx = IELEM(N,cell_index)%FACEANGLEX(face_index)
        NY = IELEM(N,cell_index)%FACEANGLEY(face_index)
		do qp_index = 1, num_qp
			qp_velocity(:) = ILOCAL_RECON3(cell_index)%QPOINTS_velocity(face_index, qp_index, 1:dimensiona)
			qp_normal_velocity = (nx*qp_velocity(1)) + (ny*qp_velocity(2))

			VolumeChangeFromVelocity_2D = VolumeChangeFromVelocity_2D + ((qp_normal_velocity*d_t)*qp_weights(qp_index)*IELEM(N,cell_index)%SURF(face_index))
		end do
	end do	

end function





SUBROUTINE COMPUTE_CENTRE_MovingMesh_2dF(N, I, facex, N_NODE, cords, node_position_index)
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
		! NODES_LIST(k,1:2)=local_nodes(IELEM(N,I)%nodes_counterclockwise(K))%positions(node_position_index,1:2)
    END DO
    CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
	! cords = compute_cell_centre(N_NODE, NODES_LIST)
   
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
		DO I=1, KMAXE
			print*,"3D part of GEOMETRY_CALC_MovingMesh not implemented yet"
			call abort()
			! CALL VOLUME_CALCULATOR3(I)
			! call SURFACE_CALCULATOR3(i)
			! CALL CENTRE3D(i)
			! call EDGE_CALCULATOR3d(i)
		END DO
		!$OMP END DO
	ELSE
		!$OMP DO
		DO I=1, KMAXE
			! CALL VOLUME_CALCULATOR_MovingMesh_2D(I, node_position_index)
			call SURFACE_CALCULATOR_MovingMesh_2D(I, node_position_index)
			! call CENTRE_MovingMesh_2D(I, node_position_index)
			call CENTREandVOLUME_MovingMesh_2D(I, node_position_index)
			call EDGE_CALCULATOR_MovingMesh_2D(I, node_position_index)
		END DO
		!$OMP END DO 
	END IF

	!$OMP BARRIER 
	!$OMP MASTER
		DUMV5 = ZERO
		DO I=1, KMAXE
			DUMV5 = DUMV5+IELEM(N,I)%moving_volume(node_position_index)
		END DO
		CALL MPI_ALLREDUCE(DUMV5,Moving_TOTALVOLUME(node_position_index),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
	!$OMP END MASTER
	!$OMP BARRIER 

END SUBROUTINE GEOMETRY_CALC_MovingMesh





SUBROUTINE GEOMETRY_CALC_MovingMesh_v2(n, node_position_index)
  !> @brief
  !> This subroutine computes the volume, surface, centre and min edge for each element
	IMPLICIT NONE
	INTEGER,INTENT(IN)::n, node_position_index
	INTEGER::KMAXE,i
	real::DUMV5

	KMAXE=XMPIELRANK(N)

	!$OMP BARRIER 

	if (DIMENSIONA.EQ.3)THEN
		!$OMP DO 
		DO I=1,KMAXE
			print*,"3D part of GEOMETRY_CALC_MovingMesh_2 not implemented yet"
			call abort()
			! CALL VOLUME_CALCULATOR3(I)
			! call SURFACE_CALCULATOR3(i)
			! CALL CENTRE3D(i)
			! call EDGE_CALCULATOR3d(i)
		END DO
		!$OMP END DO
	ELSE
		!$OMP DO
		DO I=1,KMAXE
			! CALL VOLUME_CALCULATOR_MovingMesh_2D(I, node_position_index)
			call SURFACE_CALCULATOR_MovingMesh_2D(I, node_position_index)
			call CENTRE_MovingMesh_2D(I, node_position_index)
			! call CENTREandVOLUME_MovingMesh_2D(I, node_position_index)
			call EDGE_CALCULATOR_MovingMesh_2D(I, node_position_index)
		END DO
		!$OMP END DO 
	END IF

	!$OMP BARRIER 
	!$OMP MASTER
		DUMV5=ZERO
		DO I=1, KMAXE
			DUMV5 = DUMV5+IELEM(N,I)%moving_volume(node_position_index)
		END DO
		CALL MPI_ALLREDUCE(DUMV5,Moving_TOTALVOLUME(node_position_index),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
	!$OMP END MASTER
	!$OMP BARRIER 

END SUBROUTINE GEOMETRY_CALC_MovingMesh_v2





subroutine reorder_nodes(N)
	implicit none
	integer,intent(inout)::N
	integer:: M
	integer::cell_index, node_index, my_node_index, previous_node_index, next_node_index, my_edge_index, next_edge_index
	integer::num_cells_local, num_nodes
	integer::i, i_minus, i_plus, j, to_fill
	integer::swap_helper
	real,dimension(1:dimensiona)::v1, v2
	real::v1_cross_v2

	num_cells_local = XMPIELRANK(N)

	M = omp_get_thread_num()

	if (dimensiona.ne.2) then
		print*,"node reoredering supported only in 2D"
		call abort()
	end if

	!$omp do
	do cell_index = 1, num_cells_local

		if (ielem(N, cell_index)%ishape.eq.5) then
			num_nodes = 4
		else if (ielem(N, cell_index)%ishape.eq.6) then
			num_nodes = 3
		else
			print*,"invalid cell", cell_index, "on cpu", n, "thread", M
			call abort()
		end if

		if (num_nodes.ne.ielem(N, cell_index)%nonodes) then
			print*,"missmatch between cd in cell", cell_index, "on cpu", n, "thread", M
		end if

		allocate(ielem(N, cell_index)%nodes_counterclockwise(num_nodes))
		ielem(N, cell_index)%nodes_counterclockwise(:) = -1

		my_edge_index = 1
		ielem(N, cell_index)%nodes_counterclockwise(1) = ielem(N, cell_index)%nodes_faces(1,1)
		ielem(N, cell_index)%nodes_counterclockwise(2) = ielem(N, cell_index)%nodes_faces(1,2)
		do to_fill = 3, num_nodes
			do next_edge_index = 1, num_nodes
				if (next_edge_index.ne.my_edge_index) then
					if (ielem(N, cell_index)%nodes_faces(next_edge_index,1).eq.ielem(N, cell_index)%nodes_counterclockwise(to_fill-1)) then
						ielem(N, cell_index)%nodes_counterclockwise(to_fill) = ielem(N, cell_index)%nodes_faces(next_edge_index,2)
						my_edge_index = next_edge_index
						exit
					end if
					if (ielem(N, cell_index)%nodes_faces(next_edge_index,2).eq.ielem(N, cell_index)%nodes_counterclockwise(to_fill-1)) then
						ielem(N, cell_index)%nodes_counterclockwise(to_fill) = ielem(N, cell_index)%nodes_faces(next_edge_index,1)
						my_edge_index = next_edge_index
						exit
					end if
				end if
			end do
		end do

		do i = 1, num_nodes
			if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.-1) then
				print*, "assembling counterclockwise nodes in cell", cell_index, "on CPU", N, "thread", M, "has failed"
				call abort()
			end if
		end do

		! swap order if needed
		my_node_index = ielem(N, cell_index)%nodes_counterclockwise(2)
		previous_node_index = ielem(N, cell_index)%nodes_counterclockwise(1)
		next_node_index = ielem(N, cell_index)%nodes_counterclockwise(3)

		v1 = local_nodes(my_node_index)%positions(1, 1:dimensiona) - local_nodes(previous_node_index)%positions(1, 1:dimensiona)
		v2 = local_nodes(next_node_index)%positions(1, 1:dimensiona) - local_nodes(my_node_index)%positions(1, 1:dimensiona)

		v1_cross_v2 = cross_product_2D(v1,v2)

		if (v1_cross_v2.lt.zero) then
			do i = 1, num_nodes/2
				j = num_nodes + 1 - i
				swap_helper = ielem(N, cell_index)%nodes_counterclockwise(j)
				ielem(N, cell_index)%nodes_counterclockwise(j) = ielem(N, cell_index)%nodes_counterclockwise(i)
				ielem(N, cell_index)%nodes_counterclockwise(i) = swap_helper
			end do
		end if

		do i = 1, num_nodes
			i_minus = i-1
			if (i_minus.eq.0) then
				i_minus = num_nodes
			end if
			i_plus = i+1
			if (i_plus.gt.num_nodes) then
				i_plus = 1
			end if

			my_node_index = ielem(N, cell_index)%nodes_counterclockwise(i)
			previous_node_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
			next_node_index = ielem(N, cell_index)%nodes_counterclockwise(i_plus)

			v1 = local_nodes(my_node_index)%positions(1, 1:dimensiona) - local_nodes(previous_node_index)%positions(1, 1:dimensiona)
			v2 = local_nodes(next_node_index)%positions(1, 1:dimensiona) - local_nodes(my_node_index)%positions(1, 1:dimensiona)
	
			v1_cross_v2 = cross_product_2D(v1,v2)

			if (v1_cross_v2.lt.zero) then
				if (num_nodes.eq.3) then
					print*,"concave triangle cell???", cell_index, "cpu", N, "thread", M
				else if (num_nodes.eq.4) then
					print*,"concave quadrilateral cell???", cell_index, "cpu", N, "thread", M
				else
					print*,"invalid and concave cell??????"
				end if
			end if
		end do

	end do
	!$omp end do

end subroutine


END MODULE