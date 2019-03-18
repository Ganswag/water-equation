!Se modificaron las subrutinas de Unidades, para dejar todas las unidades en 
!las originales, en éste caso presión en MPA y temperatura en K.

!Basado en el programa original de Lester Haar, Jhon S. Gallagher y George 
!S.Kell


!-------------------------PROGRAMA PRINCIPAL
! Roberto Juárez Yáñez


	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 NT,ND,NP,NH
	COMMON /UNITS/ IT,ID,IP,IH,NT,ND,NP,NH,FT,FD,FP,FH
	COMMON /QQQQ/ Q0,Q5
	COMMON /FCTS/ AD,GD,SD,UD,HD,CVD,CPD,DPDT,DVDT,CJIT,CJTH
	COMMON /ACONST/ WM,GASCON,TZ,AA,Z,DZ,Y,UREF,SREF
	COMMON /NCONST/ G(40),II(40),JJ(40),NC
	!DATA   NSS1,NSS2/2HM,2Hft/

	!CALL UNIT
	!NS=NSS1
	IF(ID .EQ. 4) NS=NSS2
15    WRITE(6,11)
100   READ(5,*,END=9) IOPT,X,T
!	T=TTT(TT)
	RT=GASCON*T
	CALL BB(T)
	IF(IOPT .LE. 0) GOTO 9
	GOTO (101,102,100,100,100), IOPT
101   DD=X
	D=DD*FD
	CALL QQ(T,D)
	ZDUM = BASE(D,T)
	PRES = FP*(RT*D*Z+Q0)
	DQ=RT*(Z+Y*DZ)+Q5
	GOTO 111
102   PRES=X
	P=PRES/FP
	DGSS=P/T/.4D0
	PSAT=2.D4
	DLL=0.D0
	DVV=0.D0
	IF(T .LT. TZ) CALL PCORR(T,PSAT,DLL,DVV)
	IF(P .GT. PSAT) DGSS=DLL
	CALL DFIND(D,P,DGSS,T,DQ)
	DD=D/FD
111    CALL THERM(D,T)
	U = UD*RT*FH
	C=DSQRT(DABS(CPD*DQ*1.D3/CVD))
	IF(ID.EQ.4) C=C*3.280833D0
	H = HD*RT*FH
	S = SD*GASCON*FH*FT
	CP=CPD*GASCON*FH*FT
	CV=CVD*GASCON*FH*FT
	VL=1.D0/D
	DPDD = DQ*FD*FP
	DPDT1=DPDT*FP*FT
	WRITE(6,20) IT,NT,PRES,NP,DD,ND
	WRITE(6,21) DPDT1,DPDD,CV,NH,NT,CP,S,H,NH,U,C,NS,CJTT,CJTH,DVDT
20    FORMAT(' T =',F12.4,' DEG ',A1,5X,'P =',F13.6,1X,A6,5X,'D =', F14.10,1X,A8)
21    FORMAT('DP/DT =',F16.9,6X, 'DP/DD =',F16.5/, &
       ' CV =',F12.6,1X,A6,A1,5X,'CP =',F12.6,6X, 's =',F12.6/, &
       ' H =',F14.6,1X,A6,5X, 'u =',F14.6,6X, 'VEL SND =',F14.6,A2,/ '/SEC',&
       ' JT(T) =',F11.5,5X, 'JT(H) =',F11.5,5X, 'DV/DT =',F12.6/)
11    FORMAT('Ingrese la Densidad en MPA, y la temperatura en K; (Ingrese 0 para salir)')
	GO TO 100
9	STOP
	END


!--------------------------CONTENIDOS

	BLOCK DATA
! 	Aquí se suplen la mayoría de los datos que se necesitan para realizar los
! 	cálculos necesarios en las subrutinas empleadas.
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /ACONST/ WM,GASCON,TZ,AA,ZB,DZB,YB,UREF,SREF
	COMMON /NCONST/ G(40),II(40),JJ(40),NC
	COMMON /ELLCON/ GI,G2,GF,BI,B2,BIT,B2T,BITT,B2TT
	COMMON /BCONST/ BP(10),BQ(10)
	COMMON /ADDCON/ ATZ(4),ADZ(4),AAT(4),AAD(4)


	DATA ATZ/2*64.D1,641.6D0,27.D1/, ADZ/3*.319D0,1.55D0/,AAT/2*2.D4,    &
 	4.D4,25.D0/,AAD/34.D0,4.D1,3.D1,1.05D3/
	DATA WM/18.0152D0/,GASCON/.461522D0/,TZ/647.073D0/,AA/1.D0/,NC/36/
	DATA UREF,SREF/-4328.455039D0,7.6180802D0/
	DATA GI,G2,GF/11.D0,44.333333333333D0,3.5D0/
	DATA BP/.7478629D0,-.3540782D0,2*0.D0,.7159876D-2,0.D0,-.3528426D-2, &
      3*0.D0/,BQ/1.1278334D0,0.D0,-.5944001D0,-5.010996D0,0.D0,           &
      .63684256D0,4*0.D0/
	DATA G/-.53062968529023D3,.22744901424408D4,.78779333020687D3,  &
       -.69830527374994D2,.17863832875422D5,-.39514731563338D5,  &
        .33803884280753D5,-.13855050202703D5,-.25637436613260D6,  &
       .48212575981415D6,-.34183016969660D6, .12223156417448D6,  &
       .11797433655832D7,-.21734810110373D7, .10829952168620D7,  &
       -.25441998064049D6,-.31377774947767D7, .52911910757704D7,  &
       -.13802577177877D7,-.25109914369001D6, .46561826115608D7,  &
       -.72752773275387D7,.41774246148294D6, .14016358244614D7,  &
       -.31555231392127D7,.47929666384584D7,.40912664781209D6,  &
       -.13626369388386D7, .69625220862664D6,-.10834900096447D7,  &
       -.22722827401688D6, .38365486000660D6, .68833257944332D4,  &
        .21757245522644D5,-.26627944829770D4,-.70730418082074D5,  &
       -.225D0,-1.68D0, .055D0,-93.0D0/
	DATA II/4*0,4*1,4*2,4*3,4*4,4*5,4*6,4*8,2*2,0,4,3*2,4/
	DATA JJ/2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7, &
       2,3,5,7,1,3*4,0,2,0,0/
	END

!-----------------------------------------------

	SUBROUTINE BB(T)
!    CALCULA LAS B'S DE LAS EQS. 3,4 USANDO LOS COEFICIENTES DESDE
!    EL BLOCKDATA, CALCULA TAMBIÉN LAS PRIMERAS Y SEGUNDAS DERIVADAS
!    W.R. PARA LAS TEMPERATURAS. LAS B'S CALCULADAS EN CM3/G.

	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /ELLCON/ GI,G2,GF,BI,B2,BIT,B2T,BITT,B2TT
	COMMON /ACONST/ WM,GASCON,TZ,AA,Z,DZ,Y,UREF,SREF
	COMMON /BCONST/ BP(10),BQ(10)
	DIMENSION V(10)

      V(1)=1.
      DO 2 I=2,10
2     V(I)= V(I-1)*TZ/T
      B1  = BP(1)+BP(2)*DLOG(1./V(2))
      B2  = BQ(1)
      B1T = BP(2)*V(2)/TZ
      B2T = 0.
      B1TT= 0.
      B2TT= 0.
      DO 4 I=3,10
      B1  = B1+BP(I)*V(I-1)
      B2  = B2+BQ(I)*V(I-1)
      B1T = B1T-(I-2)*BP(I)*V(I-1)/T
      B2T = B2T-(I-2)*BQ(I)*V(I-1)/T
      B1TT= B1TT+BP(I)*(I-2)**2*V(I-1)/T/T
4     B2TT= B2TT+BQ(I)*(I-2)**2*V(I-1)/T/T
      B1TT= B1TT-B1T/T
      B2TT= B2TT-B2T/T
      RETURN
      END

!------------------------------------------------------------------

	FUNCTION BASE(D,T)
! Calcula Z (=Pbase/(DRT) de la eq. q) (llamada BASE),
! y también Abase,Gbase,Sbase,Ubase,Hbase,CVbase, 
!  l/(DRT) * DP/DT para la base fct (llamada DPDTB).
! Las AB,GB,SB,UB,HB and CVB son adimencionales:
! AB/RT, GB/RT, SB/R, etc...

	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON /ELLCON/ G1,G2,GF,B1,B2,B1T,B2T,B1TT,B2TT
! G1,G2 y GF son los ALPHA, BETA y GAMMA de la EQ 2, incluidas en el 
! blockdata. B1 AND B2 son los "EXCLUDED VOLUME" y "2ND VIRIAL" (EQS 3 y 4) 
! Provistos por la subrutina BB(T), junto sus primeras y seegundas 
! derivadas respecto al tiempo (BIT,B2T,BITT,B2TT).

	COMMON /BASEF/ AB,GB,SB,UB,HB,CVB,DPDTB
	COMMON /ACONST/ WM,GASCON,TZ,A,Z,DZ,Y,UREF,SREF
	Y  = .25D0*B1*D
	X  = 1.D0-Y
	Z0 = (1.D0+G1*Y+G2*Y*Y)/X**3
	Z  = Z0+4.D0*Y*(B2/B1-GF)
	DZ0= (G1+2.D0*G2*Y)/X**3 + 3.D0*(1.D0+G1*Y+G2*Y*Y)/X**4
	DZ = DZ0+4.D0*(B2/B1-GF)
	AB = -DLOG(X)-(G2-1.D0)/X+28.16666667D0/X/X+4.D0*Y*(B2/B1-GF) 	+15.166666667D0 + DLOG(D*T*4.55483D0)
	GB = AB + Z
	BASE = Z
	BB2TT=T*T*B2TT
	UB = -T*B1T*(Z-1.D0-D*B2)/B1-D*T*B2T
	HB = Z+UB
	CVB=2.D0*UB+(Z0-1.D0)*((T*B1T/B1)**2-T*T*B1TT/B1) - D*(BB2TT - GF*B1TT*T*T) -(T*B1T/B1)**2*Y*DZ0
	DPDTB=BASE/T + BASE*D/Z*(DZ*B1T/4.D0+B2T-B2/B1*B1T)
	SB = UB - AB
	RETURN
	END

!------------------------------------

	SUBROUTINE QQ(T,D)
! Para unas T(K) y D(G/CM3) dadas, calcula las contrubuciones residuales
! para: Presión (Q), Fact. Helmholtz (AR), DP/DRHO (Q5),
! y también p/ la función de Gibbs: Entropía, Energía interna, Entalpía,
! Capacidad térmica isocórica y DPDT. (EQ 5)
! Terms 37 al 39 son los terminos adicionales que afectan solo
! a vecindades inmediatas al punto crítico y el term 40 es el término
! adicional cuando estamos a bajas T,y altas P.
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /RESF/ AR,GR,SR,UR,HR,CVR,DPDTR
	COMMON /QQQQ/ Q,Q5
	DIMENSION QR(11),QT(10),QZR(9),QZT(9)
	EQUIVALENCE (QR(3),QZR(1)),(QT(2),QZT(1))
	COMMON /NCONST/ G(40),II(40),JJ(40),N
	COMMON /ACONST/ WM,GASCON,TZ,AA,Z,DZ,Y,UREF,SREF
	COMMON /ADDCON/ ATZ(4),ADZ(4),AAT(4),AAD(4)
	RT=GASCON*T
	QR(1)=0.D0
	Q5 = 0.D0
	Q  = 0.D0
	AR = 0.D0
	DADT=0.D0
	CVR=0.D0
	DPDTR=0.D0
	E  = DEXP(-AA*D)
	Q10= D*D*E
	Q20= 1.D0-E
	QR(2)=Q10
	V  = TZ/T
	QT(1)=T/TZ
	DO 4 I=2,10
      QR(I+1)=QR(I)*Q20
4     QT(I)=QT(I-1)*V
      DO 10 I=1,N
      K=II(I)+1
      L=JJ(I)
      ZZ=K
      QP=G(I)*AA*QZR(K-1)*QZT(L)
      Q=Q+QP
      Q5 = Q5 + AA*(2./D-AA*(1.-E*(K-I)/Q20))*QP
      AR=AR+G(I)*QZR(K)*QZT(L)/Q10/ZZ/RT
      DFDT=Q20**K*(I-L)*QZT(L+I)/TZ/K
      D2F=L*DFDT
      DPT =DFDT*Q10*AA*K/Q20
      DADT=DADT+G(I)*DFDT
      DPDTR=DPDTR+G(I)*DPT
10    CVR=CVR+G(I)*D2F/GASCON
	QP=0.D0
	Q2A=0.D0
      DO 20 J =37,40
	IF(G(J) .EQ. 0.D0) GO TO 20
      K  = II(J)
      KM = JJ(J)
      DDZ= ADZ(J-36)
      DEL= D/DDZ-1.D0
      IF(DABS(DEL).LT. 1.D-10) DEL=1.D-10
      DD  = DEL *DEL
      EX1 = -AAD(J-36)*DEL**K
      DEX = DEXP(EX1)*DEL**KM
      ATT = AAT(J-36)
      TX  = ATZ(J-36)
      TAU = T/TX-1.D0
      EX2 = -ATT*TAU*TAU
      TEX = -DEXP(EX2)
      Q10 = DEX*TEX
      QM  = KM/DEL - K*AAD(J-36)*DEL**(K-1)
      FCT = QM*D**2*Q10/DDZ
      Q5T = FCT*(2./D+QM/DDZ)-(D/DDZ)**2*Q10*(KM/DEL/DEL+K*(K-1)*AAD(J-36)*DEL**(K-2))
      Q5  = Q5+Q5T*G(J)
      QP  = QP+G(J)*FCT
      DADT= DADT-2.D0*G(J)*ATT*TAU*Q10/TX
      DPDTR=DPDTR-2.D0*G(J)*ATT*TAU*FCT/TX
      Q2A = Q2A+T*G(J)*(4.D0*ATT*EX2+2.D0+ATT)*Q10/TX/TX
      AR  = AR+Q10*G(J)/RT
20    CONTINUE
      SR  = -DADT/GASCON
	UR  = AR+SR
	CVR=CVR+Q2A/GASCON
	Q=Q+QP
	RETURN
	END

!--------------------------------------

	SUBROUTINE DFIND(DOUT,P,D,T,DPD)
! Subrutina para calcular la densidad con entrda presión P(MPA), y
! temperatura T(K), calculada inicialmente en D(G/CM3). La salida es
! en G/CM3, también, como un biproducto, DP/DRHO es calculado
! ("DPD", MPA CM3/G)
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /QQQQ/ Q0,Q5
	COMMON /ACONST/ WM,GASCON,TZ,AA,Z,DZ,Y,UREF,SREF
	DD=D
	RT=GASCON*T
	IF (DD .LE. 0.D0) DD=1.D-8
	IF (DD .GT. 1.9D0) DD=1.9D0
	L=0
9     L=L+1
11    IF(DD .LE. 0.D0) DD=1.D-8
      IF(DD .GT. 1.9D0) DD=1.9D0
      CALL QQ(T,DD)
      PP = RT*DD*BASE(DD,T)+Q0
      DPD=RT*(Z+Y*DZ)+Q5
! las siguientes 3 lineas es para revisar la densidad negativa DP/DRHO, para
! asumir que estamos en la región de dos fases, para ajustarlo adecuadamente.
 	IF (DPD .GT. 0.D0) GO TO 13
 	IF (D .GE. .2967D0) DD=DD*1.02D0
	IF (D .LT. .2967D0) DD=DD*.98D0
	IF (L .LE. 10) GO TO 9
13    DPDX=DPD*1.1D0
      IF(DPDX .LT. .1D0) DPDX=.1D0
      DP=DABS(1.D0-PP/P)
      IF(DP .LT. 1.D-8) GO TO 20
      IF(D .GT. .3 .AND. DP .LT. 1.D-7) GO TO 20
      IF(D .GT. .7 .AND. DP .LT. 1.D-6) GO TO 20
      X=(P-PP)/DPDX
      IF(DABS(X) .GT. .1D0) X=X*.1D0/DABS(X)
      DD=DD+X
      IF(DD .LE. 0.) DD=1.D-8
19    IF(L .LE. 30) GO TO 9
20    CONTINUE
	DOUT=DD
	RETURN
	END

!----------------------------------------------------

	SUBROUTINE THERM(D,T)
	IMPLICIT REAL*8(A-H,O-Z)

	COMMON /ACONST/ WM,GASCON,TZ,AA,ZB,DZB,Y,UREF,SREF
	COMMON /QQQQ/ QP,QDP
	COMMON /BASEF/AB,GB,SB,UB,HB,CVB,DPDTB
	COMMON /RESF/AR,GR,SR,UR,HR,CVR,DPDTR
	COMMON /IDF/ AI,GI,SI,UI,HI,CVI,CPI
	COMMON /FCTS/AD,GD,SD,UD,HD,CVD,CPD,DPDT,DVDT,CJTT,CJTH
! Subrutina para calcular las funciones termidinamicas en unidades
! adimencionales (AD=A/RT. GD=G/RT, SD=S/R. UD=U/RT,
! HD=H/RT, CVD=CV/R, AND CPD=CP/R)
	CALL IDEAL(T)
	RT = GASCON*T
	Z  = ZB + QP/RT/D
	DPDD=RT*(ZB+Y*DZB) + QDP
	AD = AB + AR + AI - UREF/T + SREF
	GD = AD + Z
	UD = UB + UR + UI - UREF/T
	DPDT = RT*D*DPDTB + DPDTR
	CVD = CVB + CVR + CVI
	CPD = CVD + T*DPDT**2/(D*D*DPDD*GASCON)
	HD = UD + Z
	SD = SB + SR + SI - SREF
	DVDT=DPDT/DPDD/D/D
	CJTT = 1./D-T*DVDT
	CJTH =-CJTT/CPD/GASCON
	RETURN
	END

!------------------------------------------------------------

	FUNCTION PS(T)
! Función para calcular una aproximación para la presión del vapor. PS,
! como funcion de la temperatura. La presión del vapor s calcula acuerdo 
! con la presion del vapor predicha en la superficie con 0.02% acorde a los
! grados o m¿cerca de la temperatura critica y puede servir como un modo 
! inicial para despues refinarlo imponiendo la condición de Gl=Gv.
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(8)
	DATA A/-7.8889166D0,2.5514255D0,-6.716169D0, &
	33.239495D0,-105.38479D0,174.35319D0,-148.39348D0, &
	48.631602D0/
	IF(T .GT. 314.D0) GO TO 2
      PL=6.3573118D0-8858.843D0/T+607.56335D0*T*(-.6)
      PS=.1D0*DEXP(PL)
      RETURN
2     V=T/647.25D0
      W=DABS(1.D0-V)
      B=0.D0
	DO 4 I=1,8
      Z=I
4     B= B+A(I)*W**((Z+1.D0)/2.D0)
      Q = B/V
	PS= 22.093D0*DEXP(Q)
	RETURN
	END

!---------------------------------------

	SUBROUTINE IDEAL(T)
!   Subrutina para calcular propiedades termodinamicas para el agua
!   en la función del gas ideal de OF H.W. WOOLLEY
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON /IDF/ AI,GI,SI,UI,HI,CVI,CPI
	DIMENSION C(18) 
	DATA C/.19730271018D2, .209662681977D2,-.483429455355D0, &
	.605743189245D1 ,22.56023885D0,-9.87532442D0,-.43135538513D1, &
	.458155781D0,-.47754901883D-1, .41238460633D-2,-.27929052852D-3, &
	.14481695261D-4,-.56473658748D-6, .16200446D-7,-.3303822796D-9, &
	.451916067368D-11,-.370734122708D-13,.137546068238D-15/
	TT=T/1.D2
	TL=DLOG(TT)
	GI=-(C(1)/TT+C(2))*TL
	HI=(C(2)+C(1)*(1.D0-TL)/TT)
	CPI=C(2)-C(1)/TT
	DO 8 I=3,18
	GI=GI-C(I)*TT**(I-6)
	HI=HI+C(I)*(I-6)*TT**(I-6)
8     CPI=CPI+C(I)*(I-6)*(I-5)*TT**(I-6)
	AI=GI-1.D0
	UI=HI-1.D0
	CVI=CPI-1.D0
	SI=UI-AI
	RETURN
	END

!-------------------------------------------------------

	SUBROUTINE CORR(T,P,DL,DV,DELG)
! Para calcular, para t y p entrante cerca de la presión del vapor, las 
! correspondientes densidades del vapor y del liquido. donde 
! DELG = (GL-GV)/RT para usarse en los calculas de correcci'on, la densidad 
! del vapor para DELG = O.
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /QQQQ/ Q00,Q11
	COMMON /ACONST/ WM,GASCON,TZ,AA,ZB,DZB,YB,CREF,SREF
	COMMON /FCTS/ AD,GD,SD,UD,HD,CVD,CPD,DPDT,DVDT,CJTT,CJTH

	DLIQ=DL
	IF(DL .LE. 0.) DLIQ=1.11-.0004*T
	CALL BB(T)
	RT=GASCON*T
	CALL DFIND(DL,P,DLIQ,T,DQ)
	CALL THERM(DL,T)
	GL=GD
	DVAP=DV
	IF(DV .LE. 0.) DVAP=P/GASCON/T
	CALL DFIND(DV,P,DVAP,T,DQ)
	IF(DV.LT. 5.D-7) DV=5.D-7
	CALL THERM(DV, T)
	GV=GD
	DELG = GL-GV
	RETURN
	END

!----------------------------------------------------------

	SUBROUTINE PCORR(T,P,DL,DV)
! Subrutina para calcular la presión del vapor P y del liquido y densidad
! de vapor correspondientes a la entrada T, corrigiendo cuandoT GL-GV=O.
! La funcion PS nos da una razonablemente buena aproximación para la 
! presión del vapor, y es usada como un punto inicial para las iteraciones

	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON /ACONST/ WM,GASCON,TZ,AA,ZB,DZB,YB,UREF,SREF
	P = PS(T)
2     CALL CORR(T,P,DL,DV,DELG)
      DP=DELG*GASCON*T/(1./DV-1./DL)
      P = P+DP
	IF(DABS(DELG) .LT. 1.D-4) RETURN
	GO TO 2
	END

! Subrutinas para cambiar las unidades con que se trabaja
!
!	SUBROUTINE UNIT
!	IMPLICIT REAL*8 (A-H,O-Z)
!	REAL*8 NT,ND,NP,NH,NNT,NND,NNP,NNH
!	COMMON /UNITS/ IT,ID,IP,IH,NT,ND,NP,NH,FT,FD,FP,FH
!	DIMENSION FFD(4),FFP(5),FFH(6),NNT(4),NND(4),NNP(5),NNH(6)
!	DATA FFD/1.D-3,1.D0, .018015200, .016018D0/
!	DATA FFP/1.D0,10.D0,9.869232667D0,145.038D0, 10.1971D0/
!	DATA FFH/2*1.D0, 18.0152D0, .23884590D0,4.30285666D0,.4299226D0/
!	DATA NNT/1HK,1HC,1HR,1HF/
!	DATA NND/6Hkg/m3,6Hg/cm3,6Hmol/L,6Hlb/ft3/
!	DATA NNP/6HMPa ,6HBar,6HAtm,6HPSI,6Hkg/cm2/
!	DATA NNH/6hkJ/kg,6H J/g,6HJ/mol,6Hcal/g,7Hcal/mol,6HBTU/lb/
!	DATA A1,A2,A3,A4/8HTEMPERAT,7HDENSITY,8HPRESSURE,8HENERGY/
! 
!	WRITE(6, 11) A1
!30    WRITE(6,12)
!      READ(5,*,END=99) IT
!     IF(IT .EQ. 0) STOP
!      IF(IT .GT. 4) GOTO 30
!      NT=NNT(IT)
!      WRITE(6,11) A2
!31    WRITE(6,13)
!      READ(5,*,END=99) ID
!      IF(ID .GT. 4 .OR. ID .LT. 1) GOTO 31
!      ND=NND(ID)
!      FD=FFD(ID)
!      WRITE(6, 11) A3
!32    WRITE(6,14)
!      READ(5,*,END=99) IP
!      IF(IP .GT. 5 .OR. IP .LT. 1) GOTO 32
!      NP=NNP(IP)
!      FP=FFP(IP)
!      WRITE(6, 11) A4
!33    WRITE(6,15)
!      READ(5,*,END=99) IH
!	IF(IH .GT. 6 .OR. IH .LT. 1) GOTO 33
!	NH=NNH(IH)
!	FH=FFH(IH)
!	RETURN
!99    STOP
!11    FORMAT(' Elija unidades para',A8)
!12    FORMAT(' Elija para 1=K, 2= °C, 3= °R, 4= °F')
!13    FORMAT(' Elija para 1=KG/M3, 2=G/CM3, 3=MOL/L, 4=LB/FT3')
!14    FORMAT(' Elija para 1=MPA, 2=BAR, 3=ATM, 4=PSIA, 5=KG/CM2')
!15    FORMAT(' Elija para 1=KJ/KG, 2=J/G, 3=J/MOL, 4=CALORIA/G, 5=CALORIA/MOL, 6=BTU/LB')
!	END
!
!	FUNCTION TTT(T)
!	REAL*8 T,TTT,FT,FD,FP,FH,NT,ND,NP,NH
!	COMMON /UNITS/ IT,ID,IP,IH,NT,ND,NP,NH,FT,FD,FP,FH
!1     GO TO (1,2,3,4),IT
!      TTT=T
!	FT=1.
!	RETURN
!2     TTT=T+273.15D0
!	FT=1.
!	RETURN
!3     TTT=T/1.8D0
!	FT=.5555555555556D0
!	RETURN
!4     TTT=(T+459.67D0)/1.8D0
!	FT=.5555555555556D0
!	RETURN
!	END
!
!	FUNCTION TTI (T)
!	REAL*8 T,TTI,FT,FD,FP,FH,NT,ND,NP,NH
!	COMMON /UNITS/ IT,ID,IP,IH,NT,ND,NP,NH,FT,FD,FP,FH
!	GO TO (5,6,7,8),IT
!5     TTI=T
!	RETURN
!6     TTI=T-273.15D0
!	RETURN
!7     TTI=T*1.8D0
!	RETURN
!8     TTI=T*1.8D0-459.67D0
!	RETURN
!	END
