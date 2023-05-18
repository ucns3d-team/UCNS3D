SUBROUTINE  GET_STATES_INTERIOR(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP
L=FACEX
NGP=POINTX
I=ICONSIDERED


                      IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                    ELSE !FV

				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state

				      END IF

					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					END IF

					IF (SRF.EQ.1)THEN

                        SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)
                        CALL ROTATEF(N,TRI,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
					END IF



END SUBROUTINE GET_STATES_INTERIOR



SUBROUTINE  GET_STATES_INTERIOR2D(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP
L=FACEX
NGP=POINTX
I=ICONSIDERED


IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables)= ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES,IELEM(N,I)%INEIGHN(L), NGP)
                        ELSE !FV

                        CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
                        CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state

                        END IF
!
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					END IF


END SUBROUTINE GET_STATES_INTERIOR2D


SUBROUTINE CALCULATE_INTERIOR_VISCOUS(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex
L=FACEX
NGP=POINTX
I=ICONSIDERED



					IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
						ELSE !FV

				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state

				      END IF



				      LCVGRAD(1,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,2,L,NGP);LCVGRAD(2,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,3,L,NGP);
				      LCVGRAD(3,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,4,L,NGP);LCVGRAD(4,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP);
				      RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,IELEM(N,I)%INEIGHN(L),NGP);

					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)


					  do nvar=1,turbulenceequations+passivescalar
					  RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
					  LCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(I)%ULEFTTURBV(1:3,nvar,L,NGP)
					  end do


					END IF



END SUBROUTINE CALCULATE_INTERIOR_VISCOUS


SUBROUTINE CALCULATE_INTERIOR_VISCOUS2D(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex
L=FACEX
NGP=POINTX
I=ICONSIDERED




					  IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
						LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
					  RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP)
						ELSE !FV
							CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
				      RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP)

						end if



					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)


					  do nvar=1,turbulenceequations+passivescalar
					  RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
					  LCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,nvar,L,NGP)
					  end do


					END IF





END SUBROUTINE CALCULATE_INTERIOR_VISCOUS2D


SUBROUTINE CALCULATE_BOUNDED_VISCOUS(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex,kk

L=FACEX
NGP=POINTX
I=ICONSIDERED




						IF (DG == 1) THEN
									CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
									else
								CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
								END IF



!   CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,2,L,NGP);LCVGRAD(2,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,3,L,NGP);
				      LCVGRAD(3,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,4,L,NGP);LCVGRAD(4,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP);

					  IF (SRF.EQ.1)THEN
						SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)
						CALL ROTATEF(N,TRI,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
				END IF

					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
						do nvar=1,turbulenceequations+passivescalar
						LCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(I)%ULEFTTURBV(1:3,nvar,L,NGP)
						end do
					end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if ((ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(L))%icode.eq.50))then

								  IF (DG == 1) THEN
                                    CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                    ELSE
                                    CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                    END IF



								  RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,L,NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,L,NGP);
								  RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,L,NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,L,NGP);
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF

									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do


								    END IF
									IF(PER_ROT.EQ.1)THEN
                                        CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                        DO KK=1,4
                                            RCVGRAD(KK,1:3)=ROTATE_PER_1(RCVGRAD(KK,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                        END DO
                                        RCVGRAD_T(1,1:3)=ROTATE_PER_1(RCVGRAD_T(1,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                    END IF



								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_innerx(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)

								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)


								    RCVGRAD(:,:)=LCVGRAD(:,:)
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				  				    RCVGRAD_T(:,:)=LCVGRAD_T(:,:)
				  				    end if
				  				     if (B_CODE.eq.4)then
				  				    rightv=zero
				  				    rightv(2:4)=LCVGRAD(4,1:3)
				  				    leftv=zero

				  				    CALL ROTATEF(N,TRI,leftv,rightv,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
								    leftv(2)=-leftv(2)
								    CALL ROTATEB(N,INVTRI,rightv,leftv,ANGLE1,ANGLE2)
				  				    RCVGRAD(4,1:3)=rightv(2:4)
				  				  	end if
!
!
								  END IF
							ELSE
							      IF (DG == 1) THEN
                                CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE
							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
							      END IF
							      RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,IELEM(N,I)%INEIGHN(L),NGP);
							      RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,IELEM(N,I)%INEIGHN(L),NGP);
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
								    END IF




							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if ((ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(L))%icode.eq.50))then	!periodic in other
									  IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								ELSE

									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                END IF

									  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE
                                                        IF (IEX.EQ.1)THEN
                                                RCVGRAD(4,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
                                                        ELSE
                                                RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
                                                        END IF

										      END IF


										END DO
									  END DO






								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF



									  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO

									  IF(PER_ROT.EQ.1)THEN
										CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
										  DO KK=1,4
											  RCVGRAD(KK,1:3)=ROTATE_PER_1(RCVGRAD(KK,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
										  END DO
										  RCVGRAD_T(1,1:3)=ROTATE_PER_1(RCVGRAD_T(1,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
										END IF



								END IF
							ELSE

								   IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                ELSE
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  END IF
								  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE



										      IF (IEX.EQ.1)THEN
									  RCVGRAD(4,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      END IF

										      END IF
										END DO
									  END DO



								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF

									 DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO


								    END IF

!
							END IF
					    END IF



END SUBROUTINE CALCULATE_BOUNDED_VISCOUS



SUBROUTINE CALCULATE_BOUNDED_VISCOUS2D(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex

L=FACEX
NGP=POINTX
I=ICONSIDERED

IF (DG == 1) THEN
	CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
	else
CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
END IF


				      LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
						do nvar=1,turbulenceequations+passivescalar
						LCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,nvar,L,NGP)
						end do
					end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  !CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)

									IF (DG == 1) THEN
										CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
										ELSE
										CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
										END IF


								  RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,L,NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,L,NGP);
								  RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,L,NGP);
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF

									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do


								    END IF




								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_inner2dx(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)


								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)

								    RCVGRAD(:,:)=LCVGRAD(:,:)


								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				  				    RCVGRAD_T(:,:)=LCVGRAD_T(:,:)
				  				    end if
! 				  				     if (B_CODE.eq.4)then
!
!
! 				  				    rightv=zero
! 				  				    rightv(2:3)=LCVGRAD(3,1:2)
! 				  				    leftv=zero
!
! 				  				    CALL ROTATEF2d(N,TRI,leftv,rightv,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 								    leftv(2)=-leftv(2)
! 								    CALL ROTATEB2d(N,INVTRI,rightv,leftv,ANGLE1,ANGLE2)
! 				  				    RCVGRAD(3,1:2)=rightv(2:3)
! 				  				  	end if
!
!
								  END IF
							ELSE


								IF (DG == 1) THEN
									CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
									ELSE
									  CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
									  END IF




							      RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
							      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP);
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
								    END IF




							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									 ! CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)

									IF (DG == 1) THEN
										CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									ELSE

										  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									END IF


									  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE



										      IF (IEX.EQ.1)THEN
									  RCVGRAD(3,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      END IF

										      END IF
										END DO
									  END DO






								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF




									  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO

								END IF



								END IF
							ELSE

								  !CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								IF (DG == 1) THEN
									CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									ELSE
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									  END IF


								ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE


										      IF (IEX.EQ.1)THEN
									  RCVGRAD(3,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      END IF

										      END IF
										END DO
									  END DO



								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF

									 DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO


								    END IF

!
							END IF
					    END IF






END SUBROUTINE CALCULATE_BOUNDED_VISCOUS2D


SUBROUTINE  GET_STATES_BOUNDS(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP
L=FACEX
NGP=POINTX
I=ICONSIDERED



						IF (DG == 1) THEN
									CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
									else
								CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
								END IF

								IF (SRF.EQ.1)THEN
									SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)
									CALL ROTATEF(N,TRI,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
							END IF



								IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
										CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
									ELSE
										CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									end if
								end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  !if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
									if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then

								 IF (DG == 1) THEN
                                    CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                    ELSE
                                    CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                    END IF

									IF(PER_ROT.EQ.1)THEN
										CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
									  END IF

								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF




								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_innerx(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)

								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode



								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)


								  END IF
							ELSE


                                   IF (DG == 1) THEN
                                CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE
							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
							      END IF

								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF




							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								!if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then
								 IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								ELSE

									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                END IF
								IF(PER_ROT.EQ.1)THEN
									CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
								  END IF

								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF



								END IF
							ELSE
                                   IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                ELSE
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  END IF

!
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF

!
							END IF
					    END IF






END SUBROUTINE GET_STATES_BOUNDS



SUBROUTINE  GET_STATES_BOUNDS2D(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,IKAS
L=FACEX
NGP=POINTX
I=ICONSIDERED

IF (DG == 1) THEN
                        CLEFT(1:nof_Variables)= ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                    ELSE

				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)

                    END IF

					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
					end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                        IF (DG == 1) THEN
                                            CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                        ELSE !FV
                                            CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                        END IF
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
                                        if (icoupleturb.eq.1)then
                                        CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
                                        (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
                                        ELSE
                                        CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
                                        END IF
								    END IF





								  IKAS=1

								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_inner2dx(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)


								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
!
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
!
				  				  	 IKAS=2

								  END IF
							ELSE
                                    IF (DG == 1) THEN
                                        CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
        !
                                    ELSE !FV
                                        CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                    END IF

								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF


							       IKAS=3

							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                        IF (DG == 1) THEN
                                            CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                        else
                                            CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
 									   END IF
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF




								    END IF



								END IF
							ELSE
                                    IF (DG == 1) THEN
                                        CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                    ELSE
                                        CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                    END IF
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF


								    END IF
								   IKAS=4
!
							END IF
					    END IF







END SUBROUTINE GET_STATES_BOUNDS2D
