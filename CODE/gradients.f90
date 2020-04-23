MODULE GRADIENTS
USE LIBRARY
USE FLOW_OPERATIONS
USE BLAS95
IMPLICIT NONE


 CONTAINS
 
 
SUBROUTINE ALLGRADS_INNER(N,ICONSIDERED)
!> @brief
!> This subroutine calls the gradient approximation subroutines for every interior cell
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I
I=iconsidered
NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
imax=NUMBER_OF_NEI-1


IF (FASTEST.ne.1)THEN !GREEN GAUSS EVERYTHING


SELECT CASE(IELEM(N,I)%GGS)

    CASE(0)	!LEAST SQUARES EVERYTHING
    SELECT case (ITESTCASE)
 
 
    CASE(1,2,3)
      
    CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
    CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                IF (TURBULENCE.EQ.1)THEN
                CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                CALL COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                CALL COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                
                End if
   
    CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
    
   
   END SELECT
   
   
 CASE(1)
      SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      END SELECT
      
 CASE(2)
 
 
 
 
 
     SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
       CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END SELECT
      
      
      
 END SELECT
 
 
else
 
 SELECT case (ITESTCASE)
 
  CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)

      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      END IF
 
 
 END SELECT
 
 
 END IF



END SUBROUTINE ALLGRADS_INNER

SUBROUTINE ALLGRADS_MIX(N,ICONSIDERED)
!> @brief
!> This subroutine calls the gradient approximation subroutines for every non-interior cell
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I
I=iconsidered
NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
imax=NUMBER_OF_NEI-1


IF (FASTEST.ne.1)THEN !GREEN GAUSS EVERYTHING
SELECT CASE(IELEM(N,I)%GGS)

    CASE(0)	!LEAST SQUARES EVERYTHING
    SELECT case (ITESTCASE)
 
 
    CASE(1,2,3)
      
    CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
    CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    IF (ielem(n,i)%walls.NE.1)THEN      
                    CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    ELSE
                    CALL COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    END IF
                    IF (TURBULENCE.EQ.1)THEN
                    CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    CALL COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    
                    
                    IF (ielem(n,i)%walls.NE.1)THEN
                    CALL COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    ELSE
                    CALL COMPUTE_GRADIENTS_wall_turb_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    
                    
                    END IF
                    END IF
    CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    
   END SELECT
   
   
 CASE(1)
      SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      END SELECT
      
 CASE(2)
 
 
 
 
 
     SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      IF (ielem(n,i)%walls.NE.1)THEN      
    CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    ELSE
    CALL COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    END IF
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      END SELECT
      
      
      
 END SELECT
 
 
else
 
 SELECT case (ITESTCASE)
 
  CASE(1,2,3)
      
!       
	CALL    COMPUTE_GRADIENTS_MIX_MEAN_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      
      CALL    COMPUTE_GRADIENTS_MIX_MEAN_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL    COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      END IF
 
 
 END SELECT
 
 
 END IF



END SUBROUTINE ALLGRADS_MIX



SUBROUTINE ALLGRADS_MIX_AV(N,ICONSIDERED)
!> @brief
!> This subroutine calls the average gradient approximation subroutines for every non interior cell 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I
I=iconsidered
NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
imax=NUMBER_OF_NEI-1

 call COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS_AV(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)


END SUBROUTINE ALLGRADS_MIX_AV


SUBROUTINE ALLGRADS_INNER_AV(N,ICONSIDERED)
!> @brief
!> This subroutine calls the average gradient approximation subroutines for every interior cell 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I
I=iconsidered
NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
imax=NUMBER_OF_NEI-1

 call COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS_AV(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)


END SUBROUTINE ALLGRADS_INNER_AV



SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all 
!> @brief
!> This subroutine computes the gradients of the conserved variables of each interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				
			SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			
			
			DO K=1,3
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,3
			ILOCAL_RECON5(1)%GRADIENTS(1,K,1:nof_variables)=sOLS_F(1:nof_variables,K)
			END DO
			
end subroutine COMPUTE_GRADIENTS_INNER_MEAN_GGS
 
SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the primitive variables of each interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME


	    leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	    call cons2prim(n)
	  SOLS1(1:nof_variables)=leftv(1:nof_variables)
	  sols1(5)=leftv(5)/leftv(1)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				
			leftv(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			call cons2prim(n)
			SOLS2(1:nof_variables)=leftv(1:nof_variables)
			sols2(5)=leftv(5)/leftv(1)
			DO K=1,3
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,3
			ILOCAL_RECON3(I)%GRADs(1:3,k)=sOLS_F(2:4,K)
			ILOCAL_RECON3(I)%GRADs(4,k)=sOLS_F(5,K)
			END DO
			
			
			
			
			
			
			
end subroutine COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS



SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS_AV(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the averaged primitive variables of each interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L,IND1

if (rungekutta.eq.4)then
ind1=7
else
ind1=5
end if


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME


	    leftv(1:nof_variables)=U_C(I)%VAL(IND1,1:nof_variables)
	    call cons2prim(n)
	  SOLS1(1:nof_variables)=leftv(1:nof_variables)
	  sols1(5)=leftv(5)/leftv(1)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				
			leftv(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(IND1,1:nof_variables)
			call cons2prim(n)
			SOLS2(1:nof_variables)=leftv(1:nof_variables)
			sols2(5)=leftv(5)/leftv(1)
			DO K=1,3
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,3
			ILOCAL_RECON3(I)%GRADsAV(1:3,k)=sOLS_F(2:4,K)
			ILOCAL_RECON3(I)%GRADsAV(4,k)=sOLS_F(5,K)
			END DO
			
			
			
			
			
			
			
end subroutine COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS_AV



 
SUBROUTINE COMPUTE_GRADIENTS_MIX_MEAN_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the conserved variables of each non-interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL,TEMP_VERT
REAL::OOV2
INTEGER::I,J,K,L


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	  leftv(1:nof_variables)=SOLS1(1:nof_variables)

DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2);nz=NORMAL_ALL(3)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner(N,ICONSIDERED,FACEX)
				  CORDS=CORDINATES3(N,NODES_LIST,N_NODE)
				
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2);poz(1)=CORDS(3)
				  
				   
				  
				  
				  
				  
				  LEFTV(1:nof_variables)=SOLS1(1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS(N,B_CODE,iconsidered)
				  
				  SOLS2(1:nof_variables)=RIGHTV(1:nof_variables)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				    
			     END IF
			END IF
			
			DO K=1,3
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO


			DO K=1,3
			ILOCAL_RECON5(1)%GRADIENTS(1,K,1:nof_variables)=sOLS_F(1:nof_variables,K)
			END DO
			
			
end subroutine COMPUTE_GRADIENTS_MIX_MEAN_GGS






SUBROUTINE COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS_AV(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the averaged primitive variables of each non-interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L,IND1

if (rungekutta.eq.4)then
ind1=7
else
ind1=5
end if

I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME



	  
	  leftv(1:nof_variables)=U_C(I)%VAL(IND1,1:nof_variables)
	    call cons2prim(n)
	  SOLS1(1:nof_variables)=leftv(1:nof_variables)
	  sols1(5)=leftv(5)/leftv(1)
	  
	  
	  leftv(1:nof_variables)=U_C(I)%VAL(IND1,1:nof_variables)
	  
	  
	  
	  
DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 b_code=0

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2);nz=NORMAL_ALL(3)
			
			    
			    
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(IND1,1:nof_variables)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner(N,ICONSIDERED,FACEX)
				  CORDS(1:3)=zero
 				  CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
			  
				  Poy(1)=cords(2)
				  Pox(1)=cords(1)
				  poz(1)=cords(3)
				  
 				  leftv(1:nof_variables)=U_C(I)%VAL(IND1,1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
 				  CALL BOUNDARYS(N,B_CODE,iconsidered)
				  
				  SOLS2(1:nof_variables)=RIGHTV(1:nof_variables)
				  
				  
				  				  
				  END IF
			    ELSE
				    SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(IND1,1:nof_variables)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				    
			     END IF
			END IF
			
			  leftv(1:nof_variables)=sols2(1:nof_variables)
			call cons2prim(n)
			SOLS2(1:nof_variables)=leftv(1:nof_variables)
			sols2(5)=leftv(5)/leftv(1)
			
			
			
			
			
			DO K=1,3
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO


			
			DO K=1,3
			ILOCAL_RECON3(I)%GRADsAV(1:3,k)=sOLS_F(2:4,K)
			ILOCAL_RECON3(I)%GRADsAV(4,k)=sOLS_F(5,K)
			END DO
			
			
end subroutine COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS_AV


SUBROUTINE COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the primitive variables of each non-interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L

I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME



	  
	  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	    call cons2prim(n)
	  SOLS1(1:nof_variables)=leftv(1:nof_variables)
	  sols1(5)=leftv(5)/leftv(1)
	  
	  
	  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	  
	  
	  
	  
DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 b_code=0

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2);nz=NORMAL_ALL(3)
			
			    
			    
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner(N,ICONSIDERED,FACEX)
				  CORDS(1:3)=zero
 				  CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
			  
				  Poy(1)=cords(2)
				  Pox(1)=cords(1)
				  poz(1)=cords(3)
				  
 				  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
 				  CALL BOUNDARYS(N,B_CODE,iconsidered)
				  
				  SOLS2(1:nof_variables)=RIGHTV(1:nof_variables)
				  
				  
				  				  
				  END IF
			    ELSE
				    SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				    
			     END IF
			END IF
			
			  leftv(1:nof_variables)=sols2(1:nof_variables)
			call cons2prim(n)
			SOLS2(1:nof_variables)=leftv(1:nof_variables)
			sols2(5)=leftv(5)/leftv(1)
			
			
			
			
			
			DO K=1,3
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO


			
			DO K=1,3
			ILOCAL_RECON3(I)%GRADs(1:3,k)=sOLS_F(2:4,K)
			ILOCAL_RECON3(I)%GRADs(4,k)=sOLS_F(5,K)
			END DO
			
			
end subroutine COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS



 
SUBROUTINE COMPUTE_GRADIENTS_INNER_turb_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all 
!> @brief
!> This subroutine computes the gradients of the turbulnce variables of each interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				
			SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)
			
			
			DO K=1,3
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,3
			ILOCAL_RECON5(1)%GRADIENTS2(1,K,1:turbulenceequations+passivescalar)=sOLS_F(1:turbulenceequations+passivescalar,K)
			END DO
			
			
			IF (ICOUPLETURB.EQ.0)THEN
			DO VAR2=1,turbulenceequations+passivescalar
			ILOCAL_RECON3(I)%ULEFTTURB(VAR2,:,:)=U_Ct(I)%VAL(1,VAR2)
			END DO
			END IF
			
end subroutine COMPUTE_GRADIENTS_INNER_turb_GGS
 
SUBROUTINE COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)/U_C(I)%VAL(1,1)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				
			SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
			
			
			DO K=1,3
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			  DO VAR2=1,turbulenceequations+passivescalar
			    ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:3,VAR2)=SOLs_f(var2,1:3)
			    ILOCAL_RECON3(I)%GRADs(4+var2,1:3)=SOLs_f(var2,1:3)
			 END DO

			

			
			
			
end subroutine COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS
 
SUBROUTINE COMPUTE_GRADIENTS_MIX_turb_GGS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each non-interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL,TEMP_VERT
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)/U_C(I)%VAL(1,1)
	  

DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2);nz=NORMAL_ALL(3)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner(N,ICONSIDERED,FACEX)
				  CORDS=CORDINATES3(N,NODES_LIST,N_NODE)
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2);poz(1)=CORDS(3)
				  
				   
				  
				  
				  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  
				  cturbl(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS(N,B_CODE,iconsidered)
				  
				  SOLS2(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)/rightv(1)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)&
					/SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1)
				      END IF
				    
			     END IF
			END IF
			
			DO K=1,3
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			

			 DO K=1,3
			ILOCAL_RECON5(1)%GRADIENTS2(1,K,1:turbulenceequations+passivescalar)=sOLS_F(1:turbulenceequations+passivescalar,K)
			END DO
			
			IF (ICOUPLETURB.EQ.0)THEN
			DO VAR2=1,turbulenceequations+passivescalar
			ILOCAL_RECON3(I)%ULEFTTURB(VAR2,:,:)=U_Ct(I)%VAL(1,VAR2)
			END DO
			END IF
			
end subroutine COMPUTE_GRADIENTS_MIX_turb_GGS

SUBROUTINE COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all 
!> @brief
!> This subroutine computes the gradients of the turbulnece variables of each non-interior cell using the Green-Gauss algorithm 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,3)::SOLS_F
REAL,DIMENSION(3)::NORMAL_ALL,TEMP_VERT
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)/U_C(I)%VAL(1,1)
	  

DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=(COS(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(2)=(SIN(ANGLE1)*SIN(ANGLE2))
				NORMAL_ALL(3)=(COS(ANGLE2))
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2);nz=NORMAL_ALL(3)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/&
				  U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner(N,ICONSIDERED,FACEX)
				  CORDS=CORDINATES3(N,NODES_LIST,N_NODE)
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2);poz(1)=CORDS(3)
				  
				   
				  
				  
				  
				  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  cturbl(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS(N,B_CODE,iconsidered)
				  
				  SOLS2(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)/rightv(1)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/&
				    U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1)
				      END IF
				    
			     END IF
			END IF
			
			DO K=1,3
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

					 DO VAR2=1,turbulenceequations+passivescalar
			    ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:3,VAR2)=SOLs_f(var2,1:3)
			    ILOCAL_RECON3(I)%GRADs(4+var2,1:3)=SOLs_f(var2,1:3)
			 END DO
			 
			 
			 
			
end subroutine COMPUTE_GRADIENTS_MIX_turb_GGS_viscous

SUBROUTINE COMPUTE_GRADIENTS_MEAN_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check all
!> @brief
!> This subroutine computes the gradients of the conserved variables of each cell using the least-squares
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
   REAL,DIMENSION(nof_variables)::SOLS1
   REAL,DIMENSION(nof_variables,20)::SOLS2
   REAL,DIMENSION(imax,nof_variables,20)::MATRIX_1
   REAL,DIMENSION(NUMBER_OF_DOG,NOF_VARIABLES,20)::MATRIX_2
   REAL,DIMENSION(NUMBER_OF_DOG,NOF_VARIABLES,20)::SOL_M
   INTEGER::I,VAR2,iq,ll



   I=ICONSIDERED
   SOLS1=ZERO;
   SOLS2=ZERO
   matrix_2=zero
   sol_m=zero
 MATRIX_1=ZERO;MATRIX_2=ZERO
  
   SOLS1(1:nof_variables)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:nof_variables)
   if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
    
      DO LL=1,IELEM(N,I)%ADMIS;
      
      
        if ((ees.ne.5).or.(ll.eq.1))then
             DO IQ=1,imax
            SOLS2(1:nof_variables,ll)=U_C(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:nof_variables)
            MATRIX_1(IQ,1:nof_variables,ll)=(SOLS2(1:nof_variables,ll)-SOLS1(1:nof_variables))
           
            END DO
        
        eLSE
        
!         NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
! NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
! imax=NUMBER_OF_NEI-1
        
               DO IQ=1,numneighbours2-1
            SOLS2(1:nof_variables,ll)=U_C(ILOCAL_RECON3(I)%IHEXLC(LL,IQ+1))%VAL(1,1:nof_variables)
            MATRIX_1(IQ,1:nof_variables,ll)=(SOLS2(1:nof_variables,ll)-SOLS1(1:nof_variables))
            
            
            
            
            
         END DO
        
        END IF

        

!          DO VAR2=1,nof_variables
!             matrix_2(1:NUMBER_OF_DOG,var2)=matmul(MATRIX_1(1:imax,VAR2),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
!             SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(LL,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(1:NUMBER_OF_DOG,VAR2))
!          END DO
    end do
     DO LL=1,IELEM(N,I)%ADMIS;
        if ((ees.ne.5).or.(ll.eq.1))then
        call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1(:,:,ll),                                                &
            SOL_M(:,:,ll)                                                    &
         )
         ELSE
         call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(1:IDEGFREE2,1:numneighbours2-1,LL),                &
            MATRIX_1(1:numneighbours2-1,1:nof_variables,ll),                                                &
            SOL_M(1:IDEGFREE2,1:nof_variables,ll)                                                    &
         )
         
         END IF
      end do
      
      
       DO LL=1,IELEM(N,I)%ADMIS;
       if ((ees.ne.5).or.(ll.eq.1))then
      ILOCAL_RECON5(1)%GRADIENTS(LL,1:NUMBER_OF_DOG,1:nof_variables)=SOL_M(1:NUMBER_OF_DOG,1:nof_variables,ll)
       
      
      ELSE
      ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,1:nof_variables)=SOL_M(1:IDEGFREE2,1:nof_variables,ll)
       
      END IF
      end do

   else
     

      DO LL=1,IELEM(N,I)%ADMIS;
     
      
        if ((ees.ne.5).or.(ll.eq.1))then
   
                DO IQ=1,imax
                
                
                    IF (ILOCAL_RECON3(I)%IHEXB(LL,IQ+1).EQ.N)THEN
                        
                    SOLS2(1:nof_variables,ll)=U_C(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:nof_variables)
                    else
                    SOLS2(1:nof_variables,ll)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1),1:nof_variables)
                    end if
                    
                
                    
                    MATRIX_1(IQ,1:nof_variables,ll)=(SOLS2(1:nof_variables,ll)-SOLS1(1:nof_variables))
                    

                    
                    
                END DO
                
                
                
                
                
         ELSE
                    DO IQ=1,numneighbours2-1
                
                
                    IF (ILOCAL_RECON3(I)%IHEXBC(LL,IQ+1).EQ.N)THEN
                        
                    SOLS2(1:nof_variables,ll)=U_C(ILOCAL_RECON3(I)%IHEXLC(LL,IQ+1))%VAL(1,1:nof_variables)
                    else
                    SOLS2(1:nof_variables,ll)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXNC(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXLC(LL,IQ+1),1:nof_variables)
                    end if
                    
                    
                    MATRIX_1(IQ,1:nof_variables,ll)=(SOLS2(1:nof_variables,ll)-SOLS1(1:nof_variables))
                    

                    
                    
                    
                END DO
                
               
         
         
         END IF
         
         
        end do
        
         DO LL=1,IELEM(N,I)%ADMIS;
        if ((ees.ne.5).or.(ll.eq.1))then
        call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1(:,:,ll),                                                &
            SOL_M(:,:,ll)                                                    &
         )
         ELSE
         call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(1:IDEGFREE2,1:numneighbours2-1,LL),                &
            MATRIX_1(1:numneighbours2-1,1:nof_variables,ll),                                                &
            SOL_M(1:IDEGFREE2,1:nof_variables,ll)                                                    &
         )
         
         END IF
      end do
      
      
       DO LL=1,IELEM(N,I)%ADMIS;
       if ((ees.ne.5).or.(ll.eq.1))then
      ILOCAL_RECON5(1)%GRADIENTS(LL,1:NUMBER_OF_DOG,1:nof_variables)=SOL_M(1:NUMBER_OF_DOG,1:nof_variables,ll)
        
      ELSE
       
      ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,1:nof_variables)=SOL_M(1:IDEGFREE2,1:nof_variables,ll)
        
      END IF
      end do

   END IF

END SUBROUTINE COMPUTE_GRADIENTS_MEAN_LSQ

! 
SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the primitve variables of each interior cell using the least-squares
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(imax,nof_variables)::MATRIX_1
REAL,DIMENSION(NOF_VARIABLES,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(NUMBER_OF_DOG,NOF_VARIABLES)::SOL_M
INTEGER::I,VAR2,iq,lq,ll

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO


		ll=1    
	      
		
		if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
		MATRIX_1=ZERO;MATRIX_2=ZERO
		LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:5)
		CALL CONS2PRIM(N)
		
	       SOLS1(2:4)=LEFTV(2:4)
	       SOLS1(1)=LEFTV(5)/LEFTV(1)
	       
               DO IQ=1,imax
                LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:5)
		CALL CONS2PRIM(N)
	       SOLS2(2:4)=LEFTV(2:4)
	       SOLS2(1)=LEFTV(5)/LEFTV(1)
  	        MATRIX_1(iq,1:4)=((SOLS2(1:4)-SOLS1(1:4)))
  	        
		END DO
		
! 		DO VAR2=1,nof_variables-1
! 		
! 
! 		   matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		 
! 		  
! 		
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(1,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		
! 		END DO

                    call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )



		DO VAR2=2,4
		ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		ILOCAL_RECON5(1)%GRADIENTSTEMP(1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,1)
		
		ELSE
		
		MATRIX_1=ZERO;MATRIX_2=ZERO
		LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:5)
		CALL CONS2PRIM(N)
		
	       SOLS1(2:4)=LEFTV(2:4)
	       SOLS1(1)=LEFTV(5)/LEFTV(1)
	       
               DO IQ=1,imax
		  IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		  LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:5)
		  
		  else
		  LEFTV(1:5)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1:5)
		  end if
		  CALL CONS2PRIM(N)
	       SOLS2(2:4)=LEFTV(2:4)
	       SOLS2(1)=LEFTV(5)/LEFTV(1)
  	        MATRIX_1(iq,1:4)=((SOLS2(1:4)-SOLS1(1:4)))
		END DO
		
		
		
! 		
                   call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )


		DO VAR2=2,4
		ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		ILOCAL_RECON5(1)%GRADIENTSTEMP(1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,1)
		
		
		END IF
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS




SUBROUTINE COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check all
!> @brief
!> This subroutine computes the gradients of the primitive variables of each non-interior cell using the least-squares
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables-1)::SOLS1,SOLS2
REAL,DIMENSION(4,imax)::MATRIX_1
REAL,DIMENSION(4,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(4)::MATRIX_3
REAL,DIMENSION(NUMBER_OF_DOG,4)::SOL_M
INTEGER::I,VAR2,ii,k0,g0,ttk,ivvm,iq,lq
real::attt
integer::ll
ll=1
I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO

	    ll=1
	    K0=ILOCAL_RECON3(I)%K0
	    G0=ILOCAL_RECON3(I)%G0
	    
	     MATRIX_1=ZERO;MATRIX_2=ZERO;sol_m=zero;
		LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:5)
		CALL CONS2PRIM(N)
		
	       SOLS1(2:4)=LEFTV(2:4)
	       SOLS1(1)=LEFTV(5)/LEFTV(1)
	    
	   
	   
	   
	   
	      DO IQ=1,imax
	      if (ilocal_Recon3(i)%local.eq.1)then
	       LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:5)
	      else
		IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		LEFTV(1:5)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:5)
	    else
		LEFTV(1:5)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1:5)
	    END IF
	      end if
   
		CALL CONS2PRIM(N)
	       SOLS2(2:4)=LEFTV(2:4)
	       SOLS2(1)=LEFTV(5)/LEFTV(1)
  	        MATRIX_1(1:4,IQ)=(ILOCAL_RECON3(I)%VOLUME(1,IQ+1)*(SOLS2(1:4)-SOLS1(1:4)))
  	        
  	        MATRIX_1(2:4,IQ)=MATRIX_1(2:4,IQ)+((SOLS1(2:4)*ILOCAL_RECON3(I)%STENCILS(LL,IQ,K0))/ILOCAL_RECON3(I)%WALLCOEFF(K0))
  	        
  	        
		END DO
		matrix_3(1:4)=-sols1(1:4)
		matrix_3(1)=zero
		
		DO VAR2=1,nof_variables-1
		  MATRIX_2=ZERO
		  IF (VAR2.gt.1)THEN
		  DO IQ=1,imax
		     
		      do lq=1,NUMBER_OF_DOG-1
		      MATRIX_2(VAR2,lq)=matrix_2(var2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%VELLSQ(IQ,LQ)
		      end do
		      
		  END DO
		  ELSE
		  DO IQ=1,imax
		     do lq=1,NUMBER_OF_DOG-1
		      MATRIX_2(VAR2,lq)=matrix_2(var2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%TEMPSQ(IQ,LQ)
		      end do
		  END DO
		  end if	  
		if (var2.eq.1)then
		SOL_M(1:NUMBER_OF_DOG-1,VAR2)=MATMUL(ILOCAL_RECON3(I)%TEMPSQMAT(1:NUMBER_OF_DOG-1,1:NUMBER_OF_DOG-1),MATRIX_2(VAR2,1:NUMBER_OF_DOG-1))
		else
		SOL_M(1:NUMBER_OF_DOG-1,VAR2)=MATMUL(ILOCAL_RECON3(I)%VELINVLSQMAT(1:NUMBER_OF_DOG-1,1:NUMBER_OF_DOG-1),MATRIX_2(VAR2,1:NUMBER_OF_DOG-1))
		
		end if
		
	     END DO
		DO VAR2=2,4
		
		
		 ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,1:IDEGFREE)=-TOLBIG
		    IVVM=0
		    DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.EQ.K0) CYCLE
					  IVVM=IVVM+1
					    ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,TTK)=SOL_M(IVVM,VAR2)
		  END DO
		  ATTT=ZERO
		  ATTT=-SOLS1(VAR2)
			  DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.NE.K0) &
				  ATTT=ATTT-ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,TTK)*&
						    ILOCAL_RECON3(I)%WALLCOEFF(TTK)
			  END DO
			    ATTT=ATTT/ILOCAL_RECON3(I)%WALLCOEFF(K0)
			    ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,K0)=ATTT
		
		END DO
		
		
		ILOCAL_RECON5(1)%GRADIENTSTEMP(1:NUMBER_OF_DOG)=-TOLBIG
		    IVVM=0
		    DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.EQ.G0) CYCLE
					  IVVM=IVVM+1
					    ILOCAL_RECON5(1)%GRADIENTSTEMP(TTK)=SOL_M(IVVM,1)
		    END DO
		    ATTT=ZERO
			  DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.NE.G0) &
				  ATTT=ATTT-ILOCAL_RECON5(1)%GRADIENTSTEMP(TTK)*&
						    ILOCAL_RECON3(I)%WALLCOEFG(TTK)
			  END DO
			    ATTT=ATTT/ILOCAL_RECON3(I)%WALLCOEFG(G0)
			    ILOCAL_RECON5(1)%GRADIENTSTEMP(G0)=ATTT
			



END SUBROUTINE COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS

SUBROUTINE COMPUTE_GRADIENTS_wall_turb_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each non-interior cell using the least-squares
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(1:turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(1:turbulenceequations+passivescalar,imax)::MATRIX_1
REAL,DIMENSION(1:turbulenceequations+passivescalar,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(1:turbulenceequations+passivescalar)::MATRIX_3
REAL,DIMENSION(NUMBER_OF_DOG,1:turbulenceequations+passivescalar)::SOL_M
INTEGER::I,VAR2,ii,k0,g0,ttk,ivvm,iq,lq
real::attt
integer::ll
ll=1

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO

	    K0=ILOCAL_RECON3(I)%K0
	    
	    
	    
	     MATRIX_1=ZERO;MATRIX_2=ZERO;sol_m=zero;
		
		
		
		sols1(1:turbulenceequations+passivescalar)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:turbulenceequations+passivescalar)/&
		U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		
		
	       
	    
	   
	   
	   
	   
	      DO IQ=1,imax
	      if (ilocal_Recon3(i)%local.eq.1)then
	      
	      sols2(1:turbulenceequations+passivescalar)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:turbulenceequations+passivescalar)/&
	      U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1)
	      else
	      
		 IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		sols2(1:turbulenceequations+passivescalar)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:turbulenceequations+passivescalar)/&
		U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1)
	    else
		sols2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),6:5+turbulenceequations+passivescalar)/&
		IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1)
	    END IF
	      end if
   
		
  	        MATRIX_1(1:turbulenceequations+passivescalar,IQ)=(ILOCAL_RECON3(I)%VOLUME(1,IQ+1)*(SOLS2(1:turbulenceequations+passivescalar)-SOLS1(1:turbulenceequations+passivescalar)))
  	        MATRIX_1(1:turbulenceequations+passivescalar,IQ)=MATRIX_1(1:turbulenceequations+passivescalar,IQ)+((SOLS1(1:turbulenceequations+passivescalar)*ILOCAL_RECON3(I)%STENCILS(LL,IQ,K0))/ILOCAL_RECON3(I)%WALLCOEFF(K0))
		END DO
		matrix_3(1:turbulenceequations+passivescalar)=-sols1(1:turbulenceequations+passivescalar)
		
		if (turbulencemodel.eq.2)then
		matrix_3(2)=60.0D0*VISC/(BETA_I1*(IELEM(N,ICONSIDERED)%WallDist**2))
		end if
		
		DO VAR2=1,turbulenceequations+passivescalar
		  MATRIX_2=ZERO
		  
		  DO IQ=1,imax
		     
		      do lq=1,NUMBER_OF_DOG-1
		      MATRIX_2(VAR2,lq)=matrix_2(var2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%VELLSQ(IQ,LQ)
		      end do
		      
		  END DO
		  	  
		
		SOL_M(1:NUMBER_OF_DOG-1,VAR2)=MATMUL(ILOCAL_RECON3(I)%VELINVLSQMAT(1:NUMBER_OF_DOG-1,1:NUMBER_OF_DOG-1),MATRIX_2(VAR2,1:NUMBER_OF_DOG-1))
		
		
		
	     END DO
		
		DO VAR2=1,turbulenceequations+passivescalar
		
		
		 ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:NUMBER_OF_DOG,VAR2)=-TOLBIG
		    IVVM=0
		    DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.EQ.K0) CYCLE
					  IVVM=IVVM+1
					    ILOCAL_RECON5(1)%GRADIENTSTURB(1,TTK,VAR2)=SOL_M(IVVM,VAR2)
		  END DO
		  ATTT=ZERO
		  ATTT=-SOLS1(VAR2)
			  DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.NE.K0) &
				  ATTT=ATTT-ILOCAL_RECON5(1)%GRADIENTSTURB(1,TTK,VAR2)*&
						    ILOCAL_RECON3(I)%WALLCOEFF(TTK)
			  END DO
			    ATTT=ATTT/ILOCAL_RECON3(I)%WALLCOEFF(K0)
			    ILOCAL_RECON5(1)%GRADIENTSTURB(1,K0,VAR2)=ATTT
		
		END DO
	   
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_wall_turb_LSQ_VISCOUS





SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each cell using the least-squares
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR)::SOLS1,SOLS2
REAL,DIMENSION(imax,TURBULENCEEQUATIONS+PASSIVESCALAR)::MATRIX_1
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(NUMBER_OF_DOG,TURBULENCEEQUATIONS+PASSIVESCALAR)::SOL_M
INTEGER::I,VAR2,ll,iq

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO


	      SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
	      if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
	      
	      DO LL=1,IELEM(N,I)%ADMIS;
		MATRIX_1=ZERO;MATRIX_2=ZERO
		if ((ees.ne.5).or.(ll.eq.1))then
		DO IQ=1,imax
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
		else
                DO IQ=1,numneighbours2-1
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXLc(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
            
            
                end if
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! 
! 		
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(LL,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO
             if ((ees.ne.5).or.(ll.eq.1))then
            call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )
            else
            call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(:,:,LL),                &
            MATRIX_1(1:numneighbours2-1,1:TURBULENCEEQUATIONS+PASSIVESCALAR),                                                &
            SOL_M(1:IDEGFREE2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)                                                    &
         )
            
            
            end if


                 if ((ees.ne.5).or.(ll.eq.1))then
		ILOCAL_RECON5(1)%GRADIENTS2(LL,1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		else
		ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		
		end if
	        end do
		
	       else
	       
	       DO LL=1,IELEM(N,I)%ADMIS;
		MATRIX_1=ZERO;MATRIX_2=ZERO
		 if ((ees.ne.5).or.(ll.eq.1))then
		DO IQ=1,imax
		
		  IF (ILOCAL_RECON3(I)%IHEXB(LL,IQ+1).EQ.N)THEN
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  
		  else
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
		  end if
		   
		   MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
		else
		DO IQ=1,numneighbours2-1
		
		  IF (ILOCAL_RECON3(I)%IHEXBc(LL,IQ+1).EQ.N)THEN
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXLc(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  
		  else
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXNc(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXLc(LL,IQ+1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
		  end if
		   
		   MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
		
		
		
		end if
		
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! 		DO IQ=1,imax
! 		MATRIX_2(VAR2,1:NUMBER_OF_DOG)=MATRIX_2(VAR2,1:NUMBER_OF_DOG) + MATRIX_1(VAR2,IQ)*ILOCAL_RECON3(I)%STENCILS(LL,IQ,1:NUMBER_OF_DOG)
! 		END DO
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(LL,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO
                if ((ees.ne.5).or.(ll.eq.1))then
                call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )
                else
                call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(:,:,LL),                &
            MATRIX_1(1:numneighbours2-1,1:nof_variables),                                                &
            SOL_M(1:IDEGFREE2,1:nof_variables)                                                    &
         )
                
                
                end if

                

		  if ((ees.ne.5).or.(ll.eq.1))then
		ILOCAL_RECON5(1)%GRADIENTS2(LL,1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		else
		ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		
		end if
	            
	       end do
		
		
		END IF
		
		
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ

SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each interior cell using the least-squares
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR)::SOLS1,SOLS2
REAL,DIMENSION(imax,TURBULENCEEQUATIONS+PASSIVESCALAR)::MATRIX_1
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(NUMBER_OF_DOG,TURBULENCEEQUATIONS+PASSIVESCALAR)::SOL_M
INTEGER::I,VAR2,iq,lq,ll

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO
ideg=ielem(n,i)%idegfree
		 ll=1
		
		
		  
	      
		
		
		
		if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
		MATRIX_1=ZERO;MATRIX_2=ZERO
		SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)/&
		U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		
		
	       
	       
               DO IQ=1,imax
                SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)/&
                U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1)
		
  	        MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=((SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)))
		END DO
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! ! 		DO IQ=1,imax
! ! 		      do lq=1,ideg
! ! 		      MATRIX_2(VAR2,lq)=MATRIX_2(VAR2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%STENCILS(1,IQ,lq)
! ! 		      end do
! ! 		END DO
! 		
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(1,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO

                   call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )


		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:NUMBER_OF_DOG,VAR2)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		
		
		ELSE
		
		MATRIX_1=ZERO;MATRIX_2=ZERO
		SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
		/U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		
	       
               DO IQ=1,imax
		  IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
		/U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		  
		  else
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
		  IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1)
		  end if
		  
  	        MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=((SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)))
		END DO
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! ! 		DO IQ=1,imax
! ! 		      do lq=1,ideg
! ! 		      MATRIX_2(VAR2,lq)=MATRIX_2(VAR2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%STENCILS(1,IQ,lq)
! ! 		      end do
! ! 		END DO
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(1,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO
                   call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )


		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		ILOCAL_RECON5(1)%GRADIENTSTurb(1,1:NUMBER_OF_DOG,VAR2)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		
		
		END IF
		
		



END SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS





 SUBROUTINE ALLGRADS_INNER2d(N,ICONSIDERED)
 !> @brief
!> This subroutine calls the respective gradient approximation for every interior cell depending on its individual setting in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I
I=iconsidered
NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
imax=NUMBER_OF_NEI-1



IF (FASTEST.ne.1)THEN !GREEN GAUSS EVERYTHING
SELECT CASE(IELEM(N,I)%GGS)

    CASE(0)	!LEAST SQUARES EVERYTHING
    SELECT case (ITESTCASE)
 
 
    CASE(1,2,3)
      
    CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
    CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    IF (TURBULENCE.EQ.1)THEN
    CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    CALL COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    CALL COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    
    END IF
   
    CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
   
   
   END SELECT
   
   
 CASE(1)
      SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      END SELECT
      
 CASE(2)
 
 
 
 
 
     SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END SELECT
      
      
      
 END SELECT
 else
 
 
 
 SELECT case (ITESTCASE)
 
  CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)

      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      END IF
 
 
 END SELECT
 
 
 END IF



END SUBROUTINE ALLGRADS_INNER2d

SUBROUTINE ALLGRADS_MIX2d(N,ICONSIDERED)
 !> @brief
!> This subroutine calls the respective gradient approximation for every non-interior cell depending on its individual setting in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I
I=iconsidered
NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE
NUMBER_OF_NEI=IELEM(N,I)%inumneighbours
imax=NUMBER_OF_NEI-1



IF (FASTEST.ne.1)THEN !GREEN GAUSS EVERYTHING

SELECT CASE(IELEM(N,I)%GGS)

    CASE(0)	!LEAST SQUARES EVERYTHING
    SELECT case (ITESTCASE)
 
 
    CASE(1,2,3)
      
    CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
    CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    IF (ielem(n,i)%walls.NE.1)THEN      
                    CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    ELSE
                    CALL COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    END IF
                    IF (TURBULENCE.EQ.1)THEN
                    CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    CALL COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    
                    
                    IF (ielem(n,i)%walls.NE.1)THEN
                    CALL COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    ELSE
                    CALL COMPUTE_GRADIENTS_wall_turb_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    
                    
                    END IF
                    
                    Else
                    CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
                    
                    
                    end if
                    
    
    
   END SELECT
   
   
 CASE(1)
      SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      END SELECT
      
 CASE(2)
 
 
 
 
 
     SELECT case (ITESTCASE)
 
 
       CASE(1,2,3)
      
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      CALL COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      IF (ielem(n,i)%walls.NE.1)THEN      
    CALL COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    ELSE
    CALL COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
    END IF
      
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      END IF
      
      END SELECT
      
      
      
 END SELECT
 
 
 else
 
 SELECT case (ITESTCASE)
 
  CASE(1,2,3)
      
!       
	CALL    COMPUTE_GRADIENTS_MIX_MEAN_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      
      CASE(4)
      
      CALL    COMPUTE_GRADIENTS_MIX_MEAN_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL    COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      IF (TURBULENCE.EQ.1)THEN
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      CALL COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)
      
      END IF
 
 
 END SELECT
 
 
 END IF



END SUBROUTINE ALLGRADS_MIX2d


SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all 
!> @brief
!> This subroutine computes the gradients of the conserved variables of each interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=ANGLE2
				
				
			SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			
			
			DO K=1,2
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,2
			ILOCAL_RECON5(1)%GRADIENTS(1,K,1:nof_variables)=sOLS_F(1:nof_variables,K)
			END DO
			
end subroutine COMPUTE_GRADIENTS_INNER_MEAN_GGS2d
 
SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the primitive variables of each interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME


	    leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	    call CONS2PRIM2d(n)
	  SOLS1(1:nof_variables)=leftv(1:nof_variables)
	  sols1(4)=leftv(4)/leftv(1)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=ANGLE2
				
			leftv(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			call CONS2PRIM2d(n)
			SOLS2(1:nof_variables)=leftv(1:nof_variables)
			sols2(4)=leftv(4)/leftv(1)
			DO K=1,2
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,2
			ILOCAL_RECON3(I)%GRADs(1:2,k)=sOLS_F(2:3,K)
			ILOCAL_RECON3(I)%GRADs(3,k)=sOLS_F(4,K)
			END DO
			
			
			
			
			
			
			
end subroutine COMPUTE_GRADIENTS_INNER_MEAN_GGS_VISCOUS2d
 
SUBROUTINE COMPUTE_GRADIENTS_MIX_MEAN_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the conserved variables of each non-interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL,TEMP_VERT
REAL::OOV2
INTEGER::I,J,K,L


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	  leftv(1:nof_variables)=SOLS1(1:nof_variables)

DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=ANGLE2
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner2d(N,ICONSIDERED,FACEX)
				  CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2)
				  
				   
				  
				  
				  
				  
				  LEFTV(1:nof_variables)=SOLS1(1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS2d(N,B_CODE,iconsidered)
				  
				  SOLS2(1:nof_variables)=RIGHTV(1:nof_variables)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				    
			     END IF
			END IF
			
			DO K=1,2
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO


			DO K=1,2
			ILOCAL_RECON5(1)%GRADIENTS(1,K,1:nof_variables)=sOLS_F(1:nof_variables,K)
			END DO
			
			
end subroutine COMPUTE_GRADIENTS_MIX_MEAN_GGS2d

SUBROUTINE COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the primitive variables of each non-interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(nof_variables,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L

I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  
	  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	    call CONS2PRIM2d(n)
	  SOLS1(1:nof_variables)=leftv(1:nof_variables)
	  sols1(4)=leftv(4)/leftv(1)
	  
	  
	  leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
	  
	  
	  
	  
DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=angle2
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner2d(N,ICONSIDERED,FACEX)
				  CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2)
				  
				  
				  
				  
				  
				  
				  
				  
				  LEFTV(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS2d(N,B_CODE,iconsidered)
				  
				  SOLS2(1:nof_variables)=RIGHTV(1:nof_variables)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1:nof_variables)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1:nof_variables)
				      ELSE
					SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1:nof_variables)
				      END IF
				    
			     END IF
			END IF
			
			  leftv(1:nof_variables)=sols2(1:nof_variables)
			call CONS2PRIM2d(n)
			SOLS2(1:nof_variables)=leftv(1:nof_variables)
			sols2(4)=leftv(4)/leftv(1)
			
			
			DO K=1,2
			SOLS_F(1:nof_variables,K)=SOLS_F(1:nof_variables,K)+((OO2*(SOLS2(1:nof_variables)+SOLS1(1:nof_variables)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO


			
			DO K=1,2
			ILOCAL_RECON3(I)%GRADs(1:2,k)=sOLS_F(2:3,K)
			ILOCAL_RECON3(I)%GRADs(3,k)=sOLS_F(4,K)
			END DO
			
			
end subroutine COMPUTE_GRADIENTS_MIX_MEAN_GGS_VISCOUS2d
 
SUBROUTINE COMPUTE_GRADIENTS_INNER_turb_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all 
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each non-interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=angle2
				
			SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)
			
			
			DO K=1,2
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			DO K=1,2
			ILOCAL_RECON5(1)%GRADIENTS2(1,K,1:turbulenceequations+passivescalar)=sOLS_F(1:turbulenceequations+passivescalar,K)
			END DO
			
			
			IF (ICOUPLETURB.EQ.0)THEN
			DO VAR2=1,turbulenceequations+passivescalar
			ILOCAL_RECON3(I)%ULEFTTURB(VAR2,:,:)=U_Ct(I)%VAL(1,VAR2)
			END DO
			END IF
			
			
end subroutine COMPUTE_GRADIENTS_INNER_turb_GGS2d
 
SUBROUTINE COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)/U_C(I)%VAL(1,1)

DO J=1,IELEM(N,I)%IFCA
			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=angle2
				
			SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
			
			
			DO K=1,2
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

			  DO VAR2=1,turbulenceequations+passivescalar
			    ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:2,VAR2)=SOLs_f(var2,1:2)
			    ILOCAL_RECON3(I)%GRADs(3+var2,1:2)=SOLs_f(var2,1:2)
			 END DO

			

			
			
			
end subroutine COMPUTE_GRADIENTS_INNER_turb_GGS_VISCOUS2d
 
SUBROUTINE COMPUTE_GRADIENTS_MIX_turb_GGS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each non-interior cell using the Green-Gauss technique in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL,TEMP_VERT
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
	  LEFTV(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)

DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=angle2
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)
				  
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner2d(N,ICONSIDERED,FACEX)
				  CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2)
				  
				   
				  
				  
				  
				   LEFTV(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  cturbl(1:turbulenceequations+passivescalar)=SOLS1(1:turbulenceequations+passivescalar)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS2d(N,B_CODE,iconsidered)
				  
				  SOLS2(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END IF
				    
			     END IF
			END IF
			
			DO K=1,2
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO


			 DO K=1,2
			ILOCAL_RECON5(1)%GRADIENTS2(1,K,1:turbulenceequations+passivescalar)=sOLS_F(1:turbulenceequations+passivescalar,K)
			END DO
			
			
			
			
			IF (ICOUPLETURB.EQ.0)THEN
			DO VAR2=1,turbulenceequations+passivescalar
			ILOCAL_RECON3(I)%ULEFTTURB(VAR2,:,:)=U_Ct(I)%VAL(1,VAR2)
			END DO
			END IF
			
end subroutine COMPUTE_GRADIENTS_MIX_turb_GGS2d

SUBROUTINE COMPUTE_GRADIENTS_MIX_turb_GGS_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables of each non-interior cell using the Green-Gauss technique
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(turbulenceequations+passivescalar,2)::SOLS_F
REAL,DIMENSION(2)::NORMAL_ALL,TEMP_VERT
REAL::OOV2
INTEGER::I,J,K,L,var2


I=ICONSIDERED
SOLS_F=zero
OOV2=1.0D0/IELEM(N,I)%TOTVOLUME

	  SOLS1(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)/U_C(I)%VAL(1,1)
	  

DO J=1,IELEM(N,I)%IFCA
			 FACEX=J
			 

			ANGLE1=IELEM(N,I)%FACEANGLEX(J)
			ANGLE2=IELEM(N,I)%FACEANGLEY(J)
				NORMAL_ALL(1)=angle1
				NORMAL_ALL(2)=angle2
				nx=NORMAL_ALL(1);ny=NORMAL_ALL(2)
				
			
			IF (IELEM(N,I)%INEIGHB(J).EQ.N)THEN	!MY CPU ONLY
			    IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN MY CPU
				  SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/&
				  U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
				  ELSE
				  !NOT PERIODIC ONES IN MY CPU
				  
				  CALL coordinates_face_inner2d(N,I,FACEX)
				  CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
				  Pox(1)=CORDS(1);Poy(1)=CORDS(2);
				  
				    LEFTV(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  
				  
				  
				  
				  cturbl(1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(j))%icode
				  CALL BOUNDARYS2d(N,B_CODE,i)
				  
				  SOLS2(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)/RIGHTV(1)
				  				  
				  END IF
			    ELSE
				    SOLS2(1:turbulenceequations+passivescalar)=U_Ct(IELEM(N,I)%INEIGH(J))%VAL(1,1:turbulenceequations+passivescalar)/&
				  U_C(IELEM(N,I)%INEIGH(J))%VAL(1,1)
			    
			    
			    
			    
			    END IF
			ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
			      IF (IELEM(N,I)%IBOUNDS(J).GT.0)THEN	!CHECK FOR BOUNDARIES
				  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.5)then	!PERIODIC IN OTHER CPU
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)&
					/SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1)
				      END IF
				  END IF
			      ELSE
			      
				      IF (FASTEST.EQ.1)THEN
					SOLS2(1:turbulenceequations+passivescalar)=SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)&
					/SOLCHANGER(IELEM(N,I)%INEIGHN(J))%SOL(IELEM(N,i)%Q_FACE(j)%Q_MAPL(1),1)
				      ELSE
					SOLS2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
					IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(J)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(J)),1)
				      END IF
				    
			     END IF
			END IF
			
			DO K=1,2
			SOLS_F(1:turbulenceequations+passivescalar,K)=SOLS_F(1:turbulenceequations+passivescalar,K)+((OO2*(SOLS2(1:turbulenceequations+passivescalar)+SOLS1(1:turbulenceequations+passivescalar)))*NORMAL_ALL(K)*IELEM(N,I)%SURF(J)*OOV2)
			
			END DO
END DO

					 DO VAR2=1,turbulenceequations+passivescalar
			    ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:2,VAR2)=SOLs_f(var2,1:2)
			    ILOCAL_RECON3(I)%GRADs(3+var2,1:2)=SOLs_f(var2,1:2)
			 END DO
			 
			 
			
			
end subroutine COMPUTE_GRADIENTS_MIX_turb_GGS_viscous2d
 
SUBROUTINE COMPUTE_GRADIENTS_MEAN_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check all 
!> @brief
!> This subroutine computes the gradients of the conserved variables of each cell using the least-squares 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
   REAL,DIMENSION(imax,nof_variables)::MATRIX_1
   REAL,DIMENSION(NUMBER_OF_DOG,NOF_VARIABLES)::MATRIX_2
   REAL,DIMENSION(NUMBER_OF_DOG,NOF_VARIABLES)::SOL_M
INTEGER::I,VAR2,iq,ll

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO




	    I=ICONSIDERED
   SOLS1=ZERO;
   SOLS2=ZERO

   SOLS1(1:nof_variables)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:nof_variables)
   if (ILOCAL_RECON3(I)%LOCAL.eq.1)then

      DO LL=1,IELEM(N,I)%ADMIS;
         MATRIX_1=ZERO;MATRIX_2=ZERO
                    if ((ees.ne.5).or.(ll.eq.1))then
                    DO IQ=1,imax
                        SOLS2(1:nof_variables)=U_C(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:nof_variables)
                        MATRIX_1(IQ,1:nof_variables)=(SOLS2(1:nof_variables)-SOLS1(1:nof_variables))

                    END DO
                        ELSE
                        DO IQ=1,numneighbours2-1
                        SOLS2(1:nof_variables)=U_C(ILOCAL_RECON3(I)%IHEXLC(LL,IQ+1))%VAL(1,1:nof_variables)
                        MATRIX_1(IQ,1:nof_variables)=(SOLS2(1:nof_variables)-SOLS1(1:nof_variables))
                        

                        
                    END DO
                        
                        
                        
                        END IF
        if ((ees.ne.5).or.(ll.eq.1))then
        call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )
         

         
         
         ILOCAL_RECON5(1)%GRADIENTS(LL,1:NUMBER_OF_DOG,1:nof_variables)=SOL_M(1:NUMBER_OF_DOG,1:nof_variables)
         ELSE
         call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(1:IDEGFREE2,1:numneighbours2-1,LL),                &
            MATRIX_1(1:numneighbours2-1,1:nof_variables),                                                &
            SOL_M(1:IDEGFREE2,1:nof_variables)                                                    &
         )
        ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,1:nof_variables)=SOL_M(1:IDEGFREE2,1:nof_variables)

         
         END IF
         
         
         
      end do

   else
     

      DO LL=1,IELEM(N,I)%ADMIS;
        if ((ees.ne.5).or.(ll.eq.1))then
         DO IQ=1,imax
            IF (ILOCAL_RECON3(I)%IHEXB(LL,IQ+1).EQ.N)THEN
               SOLS2(1:nof_variables)=U_C(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:nof_variables)
            else
               SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1),1:nof_variables)
            end if

            MATRIX_1(IQ,1:nof_variables)=(SOLS2(1:nof_variables)-SOLS1(1:nof_variables))
            
            

            
            
         END DO
        ELSE
            DO IQ=1,numneighbours2-1
            IF (ILOCAL_RECON3(I)%IHEXBC(LL,IQ+1).EQ.N)THEN
               SOLS2(1:nof_variables)=U_C(ILOCAL_RECON3(I)%IHEXLC(LL,IQ+1))%VAL(1,1:nof_variables)
            else
               SOLS2(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXNC(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXLC(LL,IQ+1),1:nof_variables)
            end if

            MATRIX_1(IQ,1:nof_variables)=(SOLS2(1:nof_variables)-SOLS1(1:nof_variables))

         END DO
        
        
        
        
        END IF
         
         if ((ees.ne.5).or.(ll.eq.1))then
         
         call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )

         
         
         ILOCAL_RECON5(1)%GRADIENTS(LL,1:NUMBER_OF_DOG,1:nof_variables)=SOL_M(1:NUMBER_OF_DOG,1:nof_variables)
            ELSE
            call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(1:IDEGFREE2,1:numneighbours2-1,LL),                &
            MATRIX_1(1:numneighbours2-1,1:nof_variables),                                                &
            SOL_M(1:IDEGFREE2,1:nof_variables)                                                    &
         )
            ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,1:nof_variables)=SOL_M(1:IDEGFREE2,1:nof_variables)
            

            END IF
      END DO

   END IF
		
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_MEAN_LSQ2d
! 
SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the primitive variables for interior cells the least-squares 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(imax,nof_variables)::MATRIX_1
REAL,DIMENSION(NOF_VARIABLES,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(NUMBER_OF_DOG,NOF_VARIABLES)::SOL_M
INTEGER::I,VAR2,iq,lq,ll

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO


		    ll=1
	      
		
		if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
		MATRIX_1=ZERO;MATRIX_2=ZERO
		LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:4)
		CALL CONS2PRIM2d(N)
		
	       SOLS1(2:3)=LEFTV(2:3)
	       SOLS1(1)=LEFTV(4)/LEFTV(1)
	       
               DO IQ=1,imax
                LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:4)
		CALL CONS2PRIM2d(N)
	       SOLS2(2:3)=LEFTV(2:3)
	       SOLS2(1)=LEFTV(4)/LEFTV(1)
  	        MATRIX_1(iq,1:3)=((SOLS2(1:3)-SOLS1(1:3)))
  	        
		END DO
		
		     call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )



		DO VAR2=2,3
		ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		ILOCAL_RECON5(1)%GRADIENTSTEMP(1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,1)
		
		ELSE
		
		MATRIX_1=ZERO;MATRIX_2=ZERO
		LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:4)
		CALL CONS2PRIM2d(N)
		
	       SOLS1(2:3)=LEFTV(2:3)
	       SOLS1(1)=LEFTV(4)/LEFTV(1)
	       
               DO IQ=1,imax
		  IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		  LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:4)
		  
		  else
		  LEFTV(1:4)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1:4)
		  end if
		  CALL CONS2PRIM2d(N)
	       SOLS2(2:3)=LEFTV(2:3)
	       SOLS2(1)=LEFTV(4)/LEFTV(1)
  	        MATRIX_1(iq,1:3)=((SOLS2(1:3)-SOLS1(1:3)))
		END DO
		
		
		
! 		
                   call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )


		DO VAR2=2,3
		ILOCAL_RECON5(1)%VELOCITYDOF(VAR2-1,1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		ILOCAL_RECON5(1)%GRADIENTSTEMP(1:NUMBER_OF_DOG)=SOL_M(1:NUMBER_OF_DOG,1)
		
		
		END IF
		
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_INNER_MEAN_LSQ_VISCOUS2d

SUBROUTINE COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check all
!> @brief
!> This subroutine computes the gradients of the primitive variables for non-interior cells using the least-squares 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(nof_variables)::SOLS1,SOLS2
REAL,DIMENSION(3,imax)::MATRIX_1
REAL,DIMENSION(3,imax+NUMBEROFPOINTS2*5)::MATRIX_2
REAL,DIMENSION(3)::MATRIX_3
REAL,DIMENSION(NUMBER_OF_DOG,3)::SOL_M
INTEGER::I,VAR2,ii,k0,g0,ttk,ivvm,iq,lq
real::attt
integer::ll
ll=1

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO
ll=1
	    
	    k0=ilocal_recon3(i)%k0
	    g0=ilocal_recon3(i)%g0
	    
	    
	     MATRIX_1=ZERO;MATRIX_2=ZERO;sol_m=zero;
		LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:4)
		CALL CONS2PRIM2d(N)
		
	       SOLS1(2:4)=LEFTV(2:4)
	       SOLS1(1)=LEFTV(4)/LEFTV(1)
	    
	   
	   
	   
	   
	      DO IQ=1,imax
	      if (ilocal_Recon3(i)%local.eq.1)then
	       LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:4)
	      else
		IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		LEFTV(1:4)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:4)
	    else
		LEFTV(1:4)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1:4)
	    END IF
	      end if
   
		CALL CONS2PRIM2d(N)
	       SOLS2(2:3)=LEFTV(2:3)
	       SOLS2(1)=LEFTV(4)/LEFTV(1)
  	        MATRIX_1(1:3,IQ)=(ILOCAL_RECON3(I)%VOLUME(1,IQ+1)*(SOLS2(1:3)-SOLS1(1:3)))
  	        MATRIX_1(2:3,IQ)=MATRIX_1(2:3,IQ)+((SOLs1(2:3)*ILOCAL_RECON3(I)%STENCILS(LL,IQ,K0))/ILOCAL_RECON3(I)%WALLCOEFF(K0))
		END DO
		matrix_3(1:3)=-sols1(1:3)
		matrix_3(1)=zero
		
		DO VAR2=1,nof_variables-1
		  MATRIX_2=ZERO
		  IF (VAR2.gt.1)THEN
		  DO IQ=1,imax
		     
		      do lq=1,NUMBER_OF_DOG-1
		      MATRIX_2(VAR2,lq)=matrix_2(var2,lq)+MATRIX_1(var2,IQ)*ILOCAL_RECON3(I)%VELLSQ(IQ,lq)
		      end do
		      
		  END DO
		  ELSE
		  DO IQ=1,imax
		      do lq=1,NUMBER_OF_DOG-1
		      MATRIX_2(VAR2,lq)=matrix_2(var2,lq)+MATRIX_1(var2,IQ)*ILOCAL_RECON3(I)%tempsq(IQ,lq)
		      end do
		  END DO
		  end if	  
		if (var2.eq.1)then
		SOL_M(1:NUMBER_OF_DOG-1,VAR2)=MATMUL(ILOCAL_RECON3(I)%TEMPSQMAT(1:NUMBER_OF_DOG-1,1:NUMBER_OF_DOG-1),MATRIX_2(VAR2,1:NUMBER_OF_DOG-1))
		else
		SOL_M(1:NUMBER_OF_DOG-1,VAR2)=MATMUL(ILOCAL_RECON3(I)%VELINVLSQMAT(1:NUMBER_OF_DOG-1,1:NUMBER_OF_DOG-1),MATRIX_2(VAR2,1:NUMBER_OF_DOG-1))
		
		end if
		
	     END DO
		DO VAR2=2,3
			
		
		ILOCAL_RECON5(1)%VELOCITYDOF(var2-1,1:IDEGFREE)=-tolbig
		    IVVM=0
		    DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.EQ.K0) CYCLE
					  IVVM=IVVM+1
					    ILOCAL_RECON5(1)%VELOCITYDOF(var2-1,TTK)=SOL_M(ivvm,VAR2)
		  END DO
		  ATTT=zero
		  ATTT=-SOLs1(var2)
			  DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.NE.K0) &
				  ATTT=ATTT-ILOCAL_RECON5(1)%VELOCITYDOF(var2-1,TTK)*&
						    ILOCAL_RECON3(I)%WALLCOEFF(TTK)
			  END DO
			    ATTT=ATTT/ILOCAL_RECON3(I)%WALLCOEFF(K0)
			    ILOCAL_RECON5(1)%VELOCITYDOF(var2-1,K0)=ATTT
 
		
		
		
		
		
		
		
		
		
		END DO
		
		ILOCAL_RECON5(1)%GRADIENTSTEMP(1:NUMBER_OF_DOG)=-tolbig
	  
		    IVVM=0
		    DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.EQ.G0) CYCLE
					  IVVM=IVVM+1
					    ILOCAL_RECON5(1)%GRADIENTSTEMP(TTK)=SOL_M(ivvm,1)
		    END DO
		    ATTT=zero
			  DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.NE.G0) &
				  ATTT=ATTT-ILOCAL_RECON5(1)%GRADIENTSTEMP(TTK)*&
						    ILOCAL_RECON3(I)%WALLCOEFG(TTK)
			  END DO
			    ATTT=ATTT/ILOCAL_RECON3(I)%WALLCOEFG(G0)
			    ILOCAL_RECON5(1)%GRADIENTSTEMP(G0)=ATTT
	   
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_wall_mean_LSQ_VISCOUS2d

SUBROUTINE COMPUTE_GRADIENTS_wall_turb_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables for non-interior cells using the least-squares 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(1:turbulenceequations+passivescalar)::SOLS1,SOLS2
REAL,DIMENSION(1:turbulenceequations+passivescalar,imax)::MATRIX_1
REAL,DIMENSION(1:turbulenceequations+passivescalar,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(1:turbulenceequations+passivescalar)::MATRIX_3
REAL,DIMENSION(NUMBER_OF_DOG,1:turbulenceequations+passivescalar)::SOL_M
INTEGER::I,VAR2,ii,k0,g0,ttk,ivvm,iq,lq
real::attt
integer::ll
ll=1

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO

	    ll=1
	    
	    k0=ilocal_recon3(i)%k0
	    
	     MATRIX_1=ZERO;MATRIX_2=ZERO;sol_m=zero;
		
		
		
		sols1(1:turbulenceequations+passivescalar)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:turbulenceequations+passivescalar)/&
		U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		
		
	       
	    
	   
	   
	   
	   
	      DO IQ=1,imax
	      if (ilocal_Recon3(i)%local.eq.1)then
	      
	      sols2(1:turbulenceequations+passivescalar)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:turbulenceequations+passivescalar)/&
	      U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1)
	      else
	      
		 IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		sols2(1:turbulenceequations+passivescalar)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:turbulenceequations+passivescalar)/&
		U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1)
	    else
		sols2(1:turbulenceequations+passivescalar)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),5:4+turbulenceequations+passivescalar)/&
		IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1)
	    END IF
	      end if
   
		
  	        MATRIX_1(1:turbulenceequations+passivescalar,IQ)=(ILOCAL_RECON3(I)%VOLUME(1,IQ+1)*(SOLS2(1:turbulenceequations+passivescalar)-SOLS1(1:turbulenceequations+passivescalar)))
  	        MATRIX_1(1:turbulenceequations+passivescalar,IQ)=MATRIX_1(1:turbulenceequations+passivescalar,IQ)+((SOLs1(1:turbulenceequations+passivescalar)*ILOCAL_RECON3(I)%STENCILS(LL,IQ,K0))/ILOCAL_RECON3(I)%WALLCOEFF(K0))
		END DO
		matrix_3(1:turbulenceequations+passivescalar)=-sols1(1:turbulenceequations+passivescalar)
		if (turbulencemodel.eq.2)then
		matrix_3(2)=60.0D0*VISC/(BETA_I1*(IELEM(N,ICONSIDERED)%WallDist**2))
		end if
		
		DO VAR2=1,turbulenceequations+passivescalar
		  MATRIX_2=ZERO
		  DO IQ=1,imax
		   do lq=1,NUMBER_OF_DOG-1
		      MATRIX_2(VAR2,lq)=matrix_2(var2,lq)+MATRIX_1(var2,IQ)*ILOCAL_RECON3(I)%VELLSQ(IQ,lq)
		      end do
		 end do 	  
		
		SOL_M(1:NUMBER_OF_DOG-1,VAR2)=MATMUL(ILOCAL_RECON3(I)%VELINVLSQMAT(1:NUMBER_OF_DOG-1,1:NUMBER_OF_DOG-1),MATRIX_2(VAR2,1:NUMBER_OF_DOG-1))
		
		
		
	     END DO
				
	   
		
		DO VAR2=1,turbulenceequations+passivescalar
			
		ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:NUMBER_OF_DOG,VAR2)=-tolbig
		
		    IVVM=0
		    DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.EQ.K0) CYCLE
					  IVVM=IVVM+1
					    ILOCAL_RECON5(1)%GRADIENTSTURB(1,TTK,var2)=SOL_M(ivvm,VAR2)
		  END DO
		  ATTT=zero
		  ATTT=-SOLs1(var2)
		  if ((turbulencemodel.eq.2).and.(var2.eq.2))then
		attt=60.0D0*VISC/(BETA_I1*(IELEM(N,ICONSIDERED)%WallDist**2))-SOLs1(var2)
		end if
			  DO TTK=1,NUMBER_OF_DOG
				    IF (TTK.NE.K0) &
				  ATTT=ATTT-ILOCAL_RECON5(1)%GRADIENTSTURB(1,TTK,var2)*&
						    ILOCAL_RECON3(I)%WALLCOEFF(TTK)
			  END DO
			    ATTT=ATTT/ILOCAL_RECON3(I)%WALLCOEFF(K0)
			    ILOCAL_RECON5(1)%GRADIENTSTURB(1,k0,var2)=ATTT
  		
		
		
		END DO
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_wall_turb_LSQ_VISCOUS2d

SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI) !check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables for all cells using the least-squares 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR)::SOLS1,SOLS2
REAL,DIMENSION(imax,TURBULENCEEQUATIONS+PASSIVESCALAR)::MATRIX_1
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(NUMBER_OF_DOG,TURBULENCEEQUATIONS+PASSIVESCALAR)::SOL_M
INTEGER::I,VAR2,ll,iq

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO


	      SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
	      if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
	      
	      DO LL=1,IELEM(N,I)%ADMIS;
		MATRIX_1=ZERO;MATRIX_2=ZERO
		if ((ees.ne.5).or.(ll.eq.1))then
		DO IQ=1,imax
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
		else
                DO IQ=1,numneighbours2-1
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXLc(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
            
            
                end if
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! 
! 		
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(LL,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO

                        if ((ees.ne.5).or.(ll.eq.1))then
                        call gemm(                                                  &
                        ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
                        MATRIX_1,                                                &
                        SOL_M                                                    &
                    )
                        else
                        call gemm(                                                  &
                        ILOCAL_RECON3(I)%invmat_stenciltC(:,:,LL),                &
                        MATRIX_1(1:numneighbours2-1,1:TURBULENCEEQUATIONS+PASSIVESCALAR),                                                &
                        SOL_M(1:IDEGFREE2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)                                                    &
                    )

                        end if    

                        if ((ees.ne.5).or.(ll.eq.1))then
                        ILOCAL_RECON5(1)%GRADIENTS2(LL,1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                        
                        else
                        ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                        
                        
                        end if
	        end do
		
	       else
	       
	       DO LL=1,IELEM(N,I)%ADMIS;
		 if ((ees.ne.5).or.(ll.eq.1))then
		DO IQ=1,imax
		
		  IF (ILOCAL_RECON3(I)%IHEXB(LL,IQ+1).EQ.N)THEN
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  
		  else
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(LL,IQ+1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
		  end if
		   
		   MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
		else
		DO IQ=1,numneighbours2-1
		
		  IF (ILOCAL_RECON3(I)%IHEXBc(LL,IQ+1).EQ.N)THEN
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXLc(LL,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		  
		  else
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXNc(LL,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXLc(LL,IQ+1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
		  end if
		   
		   MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=(SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		END DO
		
		
		
		end if
		
		
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! 		DO IQ=1,imax
! 		MATRIX_2(VAR2,1:NUMBER_OF_DOG)=MATRIX_2(VAR2,1:NUMBER_OF_DOG) + MATRIX_1(VAR2,IQ)*ILOCAL_RECON3(I)%STENCILS(LL,IQ,1:NUMBER_OF_DOG)
! 		END DO
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(LL,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO

                if ((ees.ne.5).or.(ll.eq.1))then
                call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )
                else
                call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stenciltC(:,:,LL),                &
            MATRIX_1(1:numneighbours2-1,1:nof_variables),                                                &
            SOL_M(1:IDEGFREE2,1:nof_variables)                                                    &
         )
            end if


		  if ((ees.ne.5).or.(ll.eq.1))then
		ILOCAL_RECON5(1)%GRADIENTS2(LL,1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:NUMBER_OF_DOG,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		else
		ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOL_M(1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		
		end if
	            
	       end do
		
		
		END IF
		
		
		
		
		
		



END SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ2d

SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS2d(N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI)!check_all
!> @brief
!> This subroutine computes the gradients of the turbulence variables for all cells using the least-squares 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED,NUMBER_OF_DOG,NUMBER_OF_NEI
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR)::SOLS1,SOLS2
REAL,DIMENSION(imax,TURBULENCEEQUATIONS+PASSIVESCALAR)::MATRIX_1
REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR,NUMBER_OF_DOG)::MATRIX_2
REAL,DIMENSION(NUMBER_OF_DOG,TURBULENCEEQUATIONS+PASSIVESCALAR)::SOL_M
INTEGER::I,VAR2,iq,lq,ll

I=ICONSIDERED
SOLS1=ZERO;
SOLS2=ZERO
ideg=ielem(n,i)%idegfree
		 ll=1
		
		
		  
	      
		
		
		
		if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
		MATRIX_1=ZERO;MATRIX_2=ZERO
		SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)/&
		U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		
		
	       
	       
               DO IQ=1,imax
                SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)/&
                U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1)
		
  	        MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=((SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)))
		END DO
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! ! 		DO IQ=1,imax
! ! 		      do lq=1,ideg
! ! 		      MATRIX_2(VAR2,lq)=MATRIX_2(VAR2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%STENCILS(1,IQ,lq)
! ! 		      end do
! ! 		END DO
! 		
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(1,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO

                   call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )


		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		ILOCAL_RECON5(1)%GRADIENTSTURB(1,1:NUMBER_OF_DOG,VAR2)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		
		
		ELSE
		
		MATRIX_1=ZERO;MATRIX_2=ZERO
		SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
		/U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		
	       
               DO IQ=1,imax
		  IF (ILOCAL_RECON3(I)%IHEXB(1,IQ+1).EQ.N)THEN
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ+1))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
		/U_C(ILOCAL_RECON3(I)%IHEXL(1,1))%VAL(1,1)
		  
		  else
		  SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)/&
		  IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ+1))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ+1),1)
		  end if
		  
  	        MATRIX_1(iq,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=((SOLS2(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-SOLS1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)))
		END DO
! 		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! ! 		DO IQ=1,imax
! ! 		      do lq=1,ideg
! ! 		      MATRIX_2(VAR2,lq)=MATRIX_2(VAR2,lq)+MATRIX_1(VAR2,iq)*ILOCAL_RECON3(I)%STENCILS(1,IQ,lq)
! ! 		      end do
! ! 		END DO
! 		 matrix_2(var2,1:NUMBER_OF_DOG)=matmul(MATRIX_1(VAR2,1:imax),ILOCAL_RECON3(I)%STENCILS(LL,1:imax,1:NUMBER_OF_DOG))
! 		SOL_M(1:NUMBER_OF_DOG,VAR2)=MATMUL(ILOCAL_RECON3(I)%INVMAT(1,1:NUMBER_OF_DOG,1:NUMBER_OF_DOG),MATRIX_2(VAR2,1:NUMBER_OF_DOG))
! 		END DO
                   call gemm(                                                  &
            ILOCAL_RECON3(I)%invmat_stencilt(:,:,LL),                &
            MATRIX_1,                                                &
            SOL_M                                                    &
         )


		DO VAR2=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		ILOCAL_RECON5(1)%GRADIENTSTurb(1,1:NUMBER_OF_DOG,VAR2)=SOL_M(1:NUMBER_OF_DOG,VAR2)
		END DO
		
		
		END IF
		
		



END SUBROUTINE COMPUTE_GRADIENTS_TURB_LSQ_VISCOUS2d











END MODULE GRADIENTS
