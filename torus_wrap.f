      SUBROUTINE GFOSUB(ID, TIME, PAR, NPAR, DFLAG,
     &                  IFLAG, RESULT)

	IMPLICIT NONE
C
C === Type and dimension statements =======================
C
C --- External variable definitions -----------------------
C
      INTEGER                             ID
      DOUBLE PRECISION                    TIME
      DOUBLE PRECISION                    PAR( * )
      INTEGER                             NPAR
      LOGICAL                             DFLAG
      INTEGER                             IFLAG
      DOUBLE PRECISION                    RESULT(6)
C
C   ID           Identifier of calling GFORCE statement
C   TIME         Current time
C   PAR          Array containing passed parameters
C   NPAR         Number of passed parameters
C   DFLAG        Differencing flag
C   IFLAG        Initial pass flag
C   RESULT       Array (dimension 6) of computed GFORCE
C                components returned to ADAMS
C
C --- Local variable and parameter definitions ------------

C---- Passed parameters
	integer I_mar, C_mar, O_mar, splineID
	double precision R,r_, E_inf, e_star, area, l0


C---- Variables required for build in functions
	integer NSTATES
	logical ERRFLG

C-----Intermediate variables used in calculations
	logical penetration 					!used to check penetration.
	double precision I(3),O(3),IO(3),D(3),rA(3),rB(3),a_
	double precision u(3),v(3),w(3), R_mat(3,3)
	double precision l
	integer n
	

C---  Dummy variables used over and over
	double precision A, B, C, uF(3)
	double precision t1, t2,sA,sB,sold,snew

C     Assign readable variable names to passed parameters

	I_mar=nint(par(1)) ! ID of the insertion marker
	C_mar=nint(par(2)) ! ID of the center of the wrapping torus
	O_mar=nint(par(3)) ! ID of the origin marker
	R=par(4)           ! major radius of the torus
	r_=par(5)          !minor radius of the torus 
	E_inf=par(6)	 ! Tangent modulus in the linear region
	e_star=par(7)	 ! Linearity strain
	area=par(8)		 ! Cross sectional area of the ligament
	l0=par(9)		 ! Free length of the ligament
	splineID=nint(par(10)) !ID of the ADAMS spline



C ===Executable code ======================================

C --- Use TDISP to get marker translational displacements
C---- Get the componenents of I,O wrt. C in C coordinates

	CALL SYSARY('TDISP', (/I_mar,C_mar,C_mar/), 3, I, NSTATES, ERRFLG) ! Componenets I with respecto C in C
	CALL SYSARY('TDISP', (/O_mar,C_mar,C_mar/), 3, O, NSTATES, ERRFLG) ! Componenets O with respecto C in C
	IO=O-I

	A=IO(1)**2+IO(2)**2
	B=2*(IO(1)*I(1)+IO(2)*I(2))
	C=I(1)**2+I(2)**2-(R-r_)**2
	t1=(-B+sqrt(B**2-4*A*C))/(2*A)
	C=I(1)**2+I(2)**2-(R+r_)**2;
	t2=(-B+sqrt(B**2-4*A*C))/(2*A);

	a_=sqrt(A)*(t2-t1)/2
	D=(/IO(1)*(t1+t2)/2+I(1),IO(2)*(t1+t2)/2+I(2),0.d0/)

	u=(/IO(1),IO(2),0.d0/)/sqrt(IO(1)**2+IO(2)**2)
	v=(/0.d0,0.d0,1.d0/)
	w=(/u(2),-u(1),0.d0/)
	R_mat(1:3,1)=u
	R_mat(1:3,2)=v
	R_mat(1:3,3)=w

	I=matmul(transpose(R_mat),I-D)
	O=matmul(transpose(R_mat),O-D)
	IO=O-I

	A=r_**2*IO(1)**2+a_**2*IO(2)**2
	B=2*(r_**2*IO(1)*I(1)+a_**2*IO(2)*I(2))
	C=r_**2*I(1)**2+a_**2*I(2)**2-a_**2*r_**2;
	
	penetration=B**2-4*A*C>0

	
	if (penetration) then

	t1=sqrt(a_**2*I(2)**2+r_**2*I(1)**2-a_**2*r_**2) !tA
	t1=(I(2)*a_-t1)/(  r_*(I(1)+a_)  )
	sA=2*atan(t1);
	rA(1)=a_*cos(sA);
	rA(2)=r_*sin(sA);

	t2=sqrt(a_**2*(O(2))**2+r_**2*O(1)**2-a_**2*r_**2)
	t2=(O(2)*a_+t2)/(  r_*(O(1)+a_)  )
	sB=2*atan(t2);
	rB(1)=a_*cos(sB);
	rB(2)=r_*sin(sB)
	rB(3)=0

	l=0.d0
	sold=sA

	do n=1,5
	snew=sold+(sA-sB)/5
	l=l+sqrt(a_**2*(cos(sold)-cos(snew))**2+r_**2*(sin(sold)-sin(snew))**2)

	sold=snew
	
	end do

	l=l+sqrt((rA(1)-I(1))**2+(rA(2)-I(2))**2)
	l=l+sqrt((rB(1)-O(1))**2+(rB(2)-O(2))**2)
	
	uF=rB-O
	uF=uF/sqrt(uF(1)**2+uF(2)**2+uF(3)**2)
	uF=matmul(R_mat,uF)

	A=(l-l0)/l0/e_star ! dimensionless strain, e/e_star
	CALL AKISPL(A,0,splineID,0,B,ERRFLG) !dimensionless stress
	C=E_inf*e_star*area*B !force value

	else

	l=sqrt(IO(1)**2+IO(2)**2+IO(3)**2)
	uF=-matmul(R_mat,IO)
	uF=uF/sqrt(uF(1)**2+uF(2)**2+uF(3)**2)

	A=(l-l0)/l0/e_star ! dimensionless strain, e/e_star
	CALL AKISPL(A,0,splineID,0,B,ERRFLG) !dimensionless stress
	C=E_inf*e_star*area*B !force value

	end if





	RESULT(1) = C*uF(1)
	RESULT(2) = C*uF(2)
	RESULT(3) = C*uF(3)
	RESULT(4) = 0
	RESULT(5) = 0
	RESULT(6) = 0



      RETURN
      END
