! /***************************************************************************
!  *   Copyright (C) 2008 by Diego Martínez Gutiérrez                        *
!  *   diegom@yoda                                                           *
!  *                                                                         *
!  *   This program is free software; you can redistribute it and/or modify  *
!  *   it under the terms of the GNU General Public License as published by  *
!  *   the Free Software Foundation; either version 2 of the License, or     *
!  *   (at your option) any later version.                                   *
!  *                                                                         *
!  *   This program is distributed in the hope that it will be useful,       *
!  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
!  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
!  *   GNU General Public License for more details.                          *
!  *                                                                         *
!  *   You should have received a copy of the GNU General Public License     *
!  *   along with this program; if not, write to the                         *
!  *   Free Software Foundation, Inc.,                                       *
!  *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
!  ***************************************************************************/
      PROGRAM infwire
c     Programa para calcular modos normales acusticos de hilos infinitos 
c     22-Sep-2008
c     Notas:Las unidades de los valores propios (frecuencias propias) estan en S.I.
      implicit none
      integer N
      parameter (N=12)        							!antes N=12
      integer NDIM
      parameter (NDIM=((((N+1)*(N+2))*3)/2)) 					!NDIM ha de ser igual al valor maximo de INDICE_GENERAL
      integer TIPO           						!TIPO=1 cuadrado,TIPO=2 circulo,TIPO=3 hexagono
      integer i,j
      integer jx
      real*8 Rho,dim_A,dim_B,C(6,6),F,PI
      integer R(3,3) 
      integer m,ne
      real*8 q
      complex*16 UI,in_Gamma(NDIM,NDIM),E(NDIM,NDIM)	
      complex*16 omega
      complex*16 OUT_A(NDIM),OUT_B(NDIM)
     &     ,Vect_L_out(NDIM,NDIM),Vect_R_out(NDIM,NDIM)
     &     ,work(2*NDIM),rwork(8*NDIM)
      integer info
      data R/1,6,5,6,2,4,5,4,3/
      integer indice_espacial(NDIM)
      integer INDICE_GENERAL,suma,etiqueta_m(NDIM),etiqueta_n(NDIM)
     &     ,ii,jj
      integer indice_i,indice_j

      integer omega2i(NDIM)
      complex*16 omega2c(NDIM)

      common /tip/ TIPO
      UI=(0.d0,1.d0)								!Numero 'i'
      PI=dacos(-1.d0)
c     Definimos el hilo 
      open(unit=1,file='hilo_inf.dat')
c     Rho  (kg/m**3) ; dim_A (m) ; dim_B (m) ; TIPO (sin unidades:1=cuadrado,2=circulo,3=hexagono)
      read(1,*) Rho,dim_A,dim_B,TIPO
      close(1)
c     Definimos la matriz C(i,j) 
      open(unit=2,file='c.dat')
      do i=1,6
         do j=1,6
            read(2,*) ii,jj,C(i,j)
            if(ii.ne.i) pause 'No coincide i'
            if(jj.ne.j) pause 'No coincide j'
         end do
      end do
      close(2)

c     e^iqz 
c      do j=0,100
      do jx=50,50
         q=dfloat(jx)*1.4d5

         INDICE_GENERAL=0

         do i=1,3								!indice x,y,z
            do m=0,N
               do ne=0,N
                  suma=m+ne
                  if(suma.le.N) then 						!mirar si es m+n menor_igual N
                     INDICE_GENERAL=INDICE_GENERAL+1 				!entonces hay mas terminos en la matriz
                     indice_espacial(INDICE_GENERAL)=i
                     Etiqueta_m(INDICE_GENERAL)=m
                     Etiqueta_n(INDICE_GENERAL)=ne
                  endif
               enddo
            enddo
         enddo

c--------------------------------------------------------------????      
c     apanyo temporal para verificar el tamanyo de INDICE_GENERAL 
         print *, 'INDICE_GENERAL maximo=',INDICE_GENERAL,'(',NDIM,')'
c      stop
c     eliminar luego-----------------------------------------------????
c     calculamos Gamma 
         open(unit=4,file='Gamma.dat')
         do ii=1,INDICE_GENERAL
            do jj=ii,INDICE_GENERAL
               indice_i=indice_espacial(ii)
               indice_j=indice_espacial(jj)

               in_Gamma(ii,jj)=
     &              C(R(indice_i,1),R(indice_j,1)) 				!R matriz relaciona c(i,j,k,l) con C(a,b)
     &              *Etiqueta_m(ii)*Etiqueta_m(jj)
     &              *1.d0/(dim_A*dim_A)
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)-2
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj))
     &              +
     &              C(R(indice_i,2),R(indice_j,2))
     &              *Etiqueta_n(ii)*Etiqueta_n(jj)
     &              *1.d0/(dim_B*dim_B)
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj)-2)
     &              +
     &              C(R(indice_i,3),R(indice_j,3))
     &              *Q*Q
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj))
     &              +
     &              C(R(indice_i,1),R(indice_j,2))
     &              *Etiqueta_m(ii)*Etiqueta_n(jj)
     &              *1.d0/(dim_A*dim_B)
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)-1
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj)-1)
     &              +
     &              C(R(indice_i,1),R(indice_j,3))
     &              *Etiqueta_m(ii)*UI*Q
     &              *1.d0/dim_A
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)-1
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj))
     &              +
     &              C(R(indice_i,2),R(indice_j,1))
     &              *Etiqueta_n(ii)*Etiqueta_m(jj)
     &              *1.d0/(dim_A*dim_B)
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)-1
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj)-1)
     &              +
     &              C(R(indice_i,2),R(indice_j,3))
     &              *Etiqueta_n(ii)*UI*Q
     &              *1.d0/dim_B
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj)-1)
     &              +
     &              C(R(indice_i,3),R(indice_j,1))
     &              *(-1)*UI*Q*Etiqueta_m(jj)
     &              *1.d0/dim_A
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)-1
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj))
     &              +
     &              C(R(indice_i,3),R(indice_j,2))
     &              *(-1)*UI*Q*Etiqueta_n(jj)
     &              *1.d0/dim_B
     &              *F(Etiqueta_m(ii)+Etiqueta_m(jj)
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj)-1)

                in_Gamma(jj,ii)=
     &	dcmplx(dble(in_Gamma(ii,jj)),(-1.d0)*dimag(in_Gamma(ii,jj)))
                write (4,*) ii,jj,in_Gamma(ii,jj),Q

c     calculamos E 
c     E(ii,jj)=Rho*F(,)

               E(ii,jj)=0
               if(indice_i.eq.indice_j) 
     &              E(ii,jj)=Rho*F(Etiqueta_m(ii)+Etiqueta_m(jj)
     &              ,Etiqueta_n(ii)+Etiqueta_n(jj))
               E(jj,ii)=E(ii,jj)
            enddo
         enddo
         close(4)
c     Se resuelve la ecuacion de valores propios 


c------------------------------------------------------------------------
c      SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
c     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
c      CHARACTER          JOBVL, JOBVR
c      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
c      DOUBLE PRECISION   RWORK( * )
c      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
c     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
c     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices
*  (A,B), the generalized eigenvalues, and optionally, the left and/or
*  right generalized eigenvectors.
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*  singular. It is usually represented as the pair (alpha,beta), as
*  there is a reasonable interpretation for beta=0, and even for both
*  being zero.
*
*  The right generalized eigenvector v(j) corresponding to the
*  generalized eigenvalue lambda(j) of (A,B) satisfies
*
*               A * v(j) = lambda(j) * B * v(j).
*
*  The left generalized eigenvector u(j) corresponding to the
*  generalized eigenvalues lambda(j) of (A,B) satisfies
*
*               u(j)**H * A = lambda(j) * u(j)**H * B
*
*  where u(j)**H is the conjugate-transpose of u(j).
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N':  do not compute the left generalized eigenvectors;
*          = 'V':  compute the left generalized eigenvectors.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N':  do not compute the right generalized eigenvectors;
*          = 'V':  compute the right generalized eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VL, and VR.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the matrix A in the pair (A,B).
*          On exit, A has been overwritten.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the matrix B in the pair (A,B).
*          On exit, B has been overwritten.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
*          generalized eigenvalues.
*
*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio alpha/beta.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
*          If JOBVL = 'V', the left generalized eigenvectors u(j) are
*          stored one after another in the columns of VL, in the same
*          order as their eigenvalues.
*          Each eigenvector is scaled so the largest component has
*          abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVL = 'N'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL >= 1, and
*          if JOBVL = 'V', LDVL >= N.
*
*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
*          If JOBVR = 'V', the right generalized eigenvectors v(j) are
*          stored one after another in the columns of VR, in the same
*          order as their eigenvalues.
*          Each eigenvector is scaled so the largest component has
*          abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVR = 'N'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR >= 1, and
*          if JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,2*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          =1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHA(j) and BETA(j) should be
*                correct for j=INFO+1,...,N.
*          > N:  =N+1: other then QZ iteration failed in DHGEQZ,
*                =N+2: error return from DTGEVC.
*
*  =====================================================================
*     ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
*     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
      call ZGGEV('N','V',NDIM,in_Gamma,NDIM,E,NDIM,OUT_A,OUT_B,
     $    Vect_L_out,NDIM,Vect_R_out,NDIM, WORK,2*NDIM, RWORK, INFO )
c------------------------------------------------------------------------
c     call dgvlrg(NDIM,in_Gamma,NDIM,E,NDIM,out_a,out_b)
CC         print *, 'error=',info,' ::lwork',work(1) 
c     Se imprimen los valores propios en un archivo
         open(unit=3,file='salida_frecuencias.dat',status='Unknown'
     &        ,position='Append')
         write (3,*) '------------------------------------------------'
         write (3,*) 'Q=',q
         write (3,*) '------------------------------------------------'
!          print *, INDICE_GENERAL,'=',NDIM,' es correcto?'
!          print *, tipo,'=3?'
c         pause 'es todo correcto?'
         do ii=1,INDICE_GENERAL,1
            if(OUT_B(ii).ne.0.d0) omega=OUT_A(ii)/OUT_B(ii)
            write (3,*) ii,cdsqrt(omega)/(2.d0*PI)
            omega2i(ii)=ii                                 			!integer como indice_general
            omega2c(ii)=cdsqrt(omega)/(2.d0*PI)            			!complex*16 como omega
         enddo
         close(3)
c----------------------------------------------------------------------
      call reordena_R(NDIM,INDICE_GENERAL,omega2i,omega2c
     & ,Vect_R_out)
      call salida_frecuencias(NDIM,INDICE_GENERAL,q,omega2i
     & ,omega2c)
      call salida_posiciones(NDIM,INDICE_GENERAL,omega2i
     & ,Vect_R_out,Etiqueta_m,Etiqueta_n,q,TIPO)
c----------------------------------------------------------------------
      enddo                     						!del bucle do j=0,0
c      pause 'todo bien?'
      end
c-------------------------------------------------------------------
c      integer function BIFACT(mm)
      real*8 function BIFACT(mm)
c     Con N=12 BIFACT puede ser entero (los valores caen dentro del rango)
c     El doble factorial de un numero negativo par no esta definido.
c     El doble factorial de un numero negativo impar (-1)!!=1;(-3)!!=-1;(-5)!!=1/3;etc. 
      implicit none
      integer i,mm

      if(mm.eq.-3) then
c      pause 'mm=-3'
      BIFACT=-1
      return
      endif
      if(mm.eq.-2) then						!no está definido el doble factorial de números pares negativos
      BIFACT=0							!pause 'mm=-2'
      return
      endif
      BIFACT=1
      if(mm.gt.1) then
         do i=mm,1,-2
            BIFACT=BIFACT*i
         enddo
      endif
      return
      end
c     Definiciones de la funcion F() de integral (sistema infinito). 
      real*8 function F(a,b)
      implicit none
      integer a,b      !,BIFACT
      integer TIPO
      real*8 PI,alpha,BIFACT,ga2,gb2,gab3,B1,Be1	!,Be2,FACT,ga1

      common /tip/ TIPO
      PI=dacos(-1.d0)
      alpha=0.7d0
c      alpha=0.3333333333333334d0

      select case(TIPO)
      case(1)
         if(((a+1)*(b+1)).ne.0) then
      F=(1.d0+(-1.d0)**a)/(2.d0*(a+1))
     & *(1.d0+(-1.d0)**b)/(2.d0*(b+1)) 						!cuadrada
     & *(1.d0-alpha**(a+b+2))		!arreglo temporal (solo valida para cuadrados) en realidad es alpha**(a+1)*beta**(b+1) para un rectangulo en general
         else
            F=0
         endif
      case(2)
         F=(1.d0-alpha**(a+b+2))*(4.d0*PI*BIFACT(a-1)*BIFACT(b-1) 		!cilindrica
     &   /BIFACT(a+b+2))*
     &((1.d0+(-1.d0)**a)*(1.d0+(-1.d0)**b))/4.d0
      case(3)
	  call GAMMA(dfloat(a+2),ga2)
	  call GAMMA(dfloat(b+2),gb2)
	  call GAMMA(dfloat(a+b+3),gab3)
	  call HYGFX(dfloat(a+1),dfloat(0-b-1),dfloat(a+2)
     &	  ,0.5d0,Be1)
        B1=Be1/2.0d0**(a+1) 	  
c-----------------mediante hipergeometrica-------------------
	if(((a+1)*(b+1)).gt.0) then
         F=2.0d0**(b+1)/3.0d0
     & *(dfloat(1-(-1)**(a+1))/dfloat(a+1))
     & *(dfloat(1-(-1)**(b+1))/dfloat(b+1))
     & *(1.0d0-alpha**(a+b+2))
     & *(
     & 1.0d0/2.0d0**(a+b+2)
     & +(Ga2*Gb2)/Gab3
     & -B1
     &   )
		return
	else
	F=0.0d0
	endif
c------------------------------------------------------------
      case default
         pause 'error de tipo en F(...)'
      end select
      return
      end
c--------------------------------------------------------------------
      subroutine reordena_R(NDIM,INDICE_GENERAL,omega2i,omega2c
     & ,Vect_R_out)
      implicit none
      integer ii,jj,INDICE_GENERAL,mmm,kkk
      integer NDIM
      integer omega2i(NDIM),omega_tempi
      complex*16 omega2c(NDIM),omega_tempc,Vect_R_out_temp(NDIM)
     & ,Vect_R_out(NDIM,NDIM)

c     Reordena los valores propios (y los vectores propios -26-11-2008)
         do ii=1,(INDICE_GENERAL-1),1
            mmm=ii
            do jj=ii+1,INDICE_GENERAL,1
              if(abs(omega2c(jj)).lt.abs(omega2c(mmm))) mmm=jj
c               if(cabs(omega2c(jj)).lt.cabs(omega2c(mmm))) mmm=jj
            enddo
            omega_tempi=omega2i(ii)
            omega_tempc=omega2c(ii)
            omega2i(ii)=omega2i(mmm)
            omega2c(ii)=omega2c(mmm)
            omega2i(mmm)=omega_tempi
            omega2c(mmm)=omega_tempc
		do kkk=1,NDIM,1
	    Vect_R_out_temp(kkk)=Vect_R_out(kkk,ii)
	    Vect_R_out(kkk,ii)=Vect_R_out(kkk,mmm)
	    Vect_R_out(kkk,mmm)=Vect_R_out_temp(kkk)
		enddo
         enddo
      return
      end

      subroutine salida_frecuencias(NDIM,INDICE_GENERAL,q,omega2i
     & ,omega2c)
      implicit none
      integer NDIM
      integer ii,j,INDICE_GENERAL
      integer omega2i(NDIM)
      complex*16 omega2c(NDIM)
      real*8 q
      character(4) pj

      do j=1,20		!antes 100
      write(pj,'(I4)') j
         open(unit=4,file=('salida_'//pj//
     & '_frecuencia_ordenada.dat')
     & ,status='Unknown',position='Append')
         write (4,*) q,
     &	 ' : ',omega2i(j),' = '
     &	 ,real(omega2c(j))
	 close(4)
      enddo

         open(unit=3,file='salida_frecuencias_ordenada.dat'
     & ,status='Unknown',position='Append')
         write (3,*) ' '
         write (3,*) ' '
         write (3,*) ' '
         write (3,*) '#------------------------------------------------'
         write (3,*) '#Q=',q
         write (3,*) '#------------------------------------------------'
         do ii=1,INDICE_GENERAL,1
            write (3,'(I3,A2,E11.4,A2,I3,A2,E11.4)') 
     &	    ii,':',q,':',omega2i(ii),'=',real(omega2c(ii))
         enddo
         close(3)
      return
      end

      subroutine salida_posiciones(NDIM,INDICE_GENERAL,omega2i
     & ,Vect_R_out,Etiqueta_m,Etiqueta_n,q,TIPO)
      implicit none
      integer NDIM,x,y,z,TIPO
      integer ii,jj,j,INDICE_GENERAL,etiqueta_m(NDIM),etiqueta_n(NDIM)
      integer omega2i(NDIM)
      complex*16 Vect_R_out(NDIM,NDIM),Uve_1,Uve_2,Uve_3,u
      real*8 q,Alpha		!,Uve_1,Uve_2,Uve_3
      character(4) pj
      u=(0.d0,1.d0)
      Alpha=0.7d0
      do j=1,20

        do ii=1,NDIM,1

      write(pj,'(I4)') j
         open(unit=4,file=('salida_'//pj//
     & '_posiciones_ordenada.dat')
     & ,status='Unknown',position='Append')
         write (4,'(I3,I3,E11.4,A2,E11.4,A4,E11.4)') 
     &   Etiqueta_m(ii),Etiqueta_n(ii),q,
     &	 ' : ',omega2i(j),' = ',real(Vect_R_out(ii,j))
	 close(4)
      enddo
      enddo

         open(unit=3,file='salida_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')
         write (3,*) ' '
         write (3,*) ' '
         write (3,*) ' '
         write (3,*) '#------------------------------------------------'
         write (3,*) '#Q=',q
         write (3,*) '#------------------------------------------------'
         do j=1,NDIM,1
         do ii=1,INDICE_GENERAL,1
            write (3,'(I3,A2,I3,I3,A4,E11.4,A2,I3,A2,E11.4)') 
     &	    ii,'(',Etiqueta_m(j),Etiqueta_n(j),'):',q,':'
     &      ,omega2i(ii),'=',real(Vect_R_out(j,ii))
         enddo
         enddo
         close(3)
      do j=1,10,1
      do x=-40,40,1
      do y=-40,40,1
      do z=1,1,1
      Uve_1=(0.0d0,0.d0)
      Uve_2=(0.0d0,0.d0)
      Uve_3=(0.0d0,0.d0)
      do jj=1,NDIM/3,1
      Uve_1=Uve_1
     & +Vect_R_out(jj+(0*NDIM/3),j)
     & *dcmplx((float(x)/40.0)**Etiqueta_m(jj+(0*NDIM/3)))
     & *dcmplx((float(y)/40.0)**Etiqueta_n(jj+(0*NDIM/3)))
     & *exp(u*dcmplx(q)*dcmplx(z))
      Uve_2=Uve_2
     & +Vect_R_out(jj+(1*NDIM/3),j)
     & *dcmplx((float(x)/40.0)**Etiqueta_m(jj+(1*NDIM/3)))
     & *dcmplx((float(y)/40.0)**Etiqueta_n(jj+(1*NDIM/3)))
     & *exp(u*dcmplx(q)*dcmplx(z))
      Uve_3=Uve_3
     & +Vect_R_out(jj+(2*NDIM/3),j)
     & *dcmplx((float(x)/40.0)**Etiqueta_m(jj+(2*NDIM/3)))
     & *dcmplx((float(y)/40.0)**Etiqueta_n(jj+(2*NDIM/3)))
     & *exp(u*dcmplx(q)*dcmplx(z))
      enddo

      select case (TIPO)
c---------------------------------------------------------------------
      case(1)		!Seccion Cuadrada
      write(pj,'(I4)') j
      open(unit=5,file='salida_modulo'//pj//
     &  '_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')
      if((abs(x).gt.Alpha*40).or.(abs(y).gt.Alpha*40)) then	!se hace un agujero
	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
      else
	    write(5,*) x,y,0.0
      endif
      close(5)
      open(unit=5,file='salida_xyz'//pj//
     & '_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')
      if((abs(x).gt.Alpha*40).or.(abs(y).gt.Alpha*40)) then
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
      else
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
      endif
      close(5)
c---------------------------------------------------------------------
      case(2)		!Seccion circular
      write(pj,'(I4)') j
      open(unit=5,file='salida_modulo'//pj//
     &  '_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')
      if(y.gt.Alpha*40) then
         if((x.lt. sqrt(40.0**2-y**2)).and.
     & (x.gt.-sqrt(40.0**2-y**2))    ) then
	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
         else
	    write(5,*) x,y,0.0
         endif
      elseif((y.gt.0).and.(y.le.Alpha*40)) then
         if(
     & (	  
     &	     (x.lt. sqrt(40.0**2-y**2))
     &  .and.(x.gt. sqrt((Alpha*40)**2-y**2))
     & ) 
     & .or.
     & (
     &       (x.gt.-sqrt(40.0**2-y**2))
     &  .and.(x.lt.-sqrt((Alpha*40)**2-y**2))
     & )
     &     ) then
	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
         else
	    write(5,*) x,y,0.0
         endif
      elseif(y.lt.-Alpha*40) then
         if((x.lt. sqrt(40.0**2-y**2)).and.
     & (x.gt.-sqrt(40.0**2-y**2))    ) then
	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
         else
	    write(5,*) x,y,0.0
         endif
      elseif((y.le.0).and.(y.ge.-Alpha*40)) then
         if(
     & (	  
     &	     (x.lt. sqrt(40.0**2-y**2))
     &  .and.(x.gt. sqrt((Alpha*40)**2-y**2))
     & ) 
     & .or.
     & (
     &       (x.gt.-sqrt(40.0**2-y**2))
     &  .and.(x.lt.-sqrt((Alpha*40)**2-y**2))
     & )
     &     ) then
	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
         else
	    write(5,*) x,y,0.0
         endif
       else
	    write(5,*) x,y,0.0
       endif
      close(5)
      open(unit=5,file='salida_xyz'//pj//
     & '_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')
      if(y.gt.Alpha*40) then
         if((x.lt. sqrt(40.0**2-y**2)).and.
     & (x.gt.-sqrt(40.0**2-y**2))    ) then
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
         else
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
         endif
      elseif((y.gt.0).and.(y.le.Alpha*40)) then
         if(
     & (	  
     &	     (x.lt. sqrt(40.0**2-y**2))
     &  .and.(x.gt. sqrt((Alpha*40)**2-y**2))
     & ) 
     & .or.
     & (
     &       (x.gt.-sqrt(40.0**2-y**2))
     &  .and.(x.lt.-sqrt((Alpha*40)**2-y**2))
     & )
     &     ) then
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
         else
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
         endif
      elseif(y.lt.-Alpha*40) then
         if((x.lt. sqrt(40.0**2-y**2)).and.
     & (x.gt.-sqrt(40.0**2-y**2))    ) then
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
         else
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
         endif
      elseif((y.le.0).and.(y.ge.-Alpha*40)) then
         if(
     & (	  
     &	     (x.lt. sqrt(40.0**2-y**2))
     &  .and.(x.gt. sqrt((Alpha*40)**2-y**2))
     & ) 
     & .or.
     & (
     &       (x.gt.-sqrt(40.0**2-y**2))
     &  .and.(x.lt.-sqrt((Alpha*40)**2-y**2))
     & )
     &     ) then
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
         else
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
         endif
       else
	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
       endif
      close(5)
c---------------------------------------------------------------------
      case(3)		!Seccion hexagonal
      write(pj,'(I4)') j
      open(unit=5,file='salida_modulo'//pj//
     &  '_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')

      if(x.lt.-20) then
        if((y.lt.2*(40+x)).and.(y.gt.(-2)*(40+x))) then
! -------------------------------------------------------------------
         if(x.lt.-20*Alpha) then
          if((y.gt.2*(Alpha*40+x)).or.(y.lt.(-2)*(Alpha*40+x))) then
      	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
          else
      	    write(5,*) x,y,0.0
          endif
         endif
! -------------------------------------------------------------------
      	else
      	    write(5,*) x,y,0.0
      	endif
      elseif((x.ge.-20).and.(x.le.20)) then
! -------------------------------------------------------------------
       if(x.lt.-20*Alpha) then
        if((y.gt.2*(Alpha*40+x)).or.(y.lt.(-2)*(Alpha*40+x))) then
      	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
      	else
      	    write(5,*) x,y,0.0
      	endif
       elseif(x.gt.20*Alpha) then
        if((y.gt.2*(Alpha*40-x)).or.(y.lt.(-2)*(Alpha*40-x))) then
      	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
      	else
      	    write(5,*) x,y,0.0
      	endif
       else
        if(abs(y).gt.Alpha*40) then
        write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
      	else
      	    write(5,*) x,y,0.0
      	endif
       endif
! -------------------------------------------------------------------
      else
        if((y.lt.2*(40-x)).and.(y.gt.(-2)*(40-x))) then
! -------------------------------------------------------------------
         if(x.gt.20*Alpha) then
          if((y.gt.2*(Alpha*40-x)).or.(y.lt.(-2)*(Alpha*40-x))) then
      	    write(5,*) x,y,dsqrt(dreal(Uve_1)**2+dimag(Uve_1)**2
     &                         + dreal(Uve_2)**2+dimag(Uve_2)**2
     &                         + dreal(Uve_3)**2+dimag(Uve_3)**2)
          else
      	    write(5,*) x,y,0.0
          endif
         endif
! -------------------------------------------------------------------
      	else
      	    write(5,*) x,y,0.0
      	endif
      endif
      close(5)
      
      open(unit=5,file='salida_xyz'//pj//
     & '_posiciones_ordenada.dat'
     & ,status='Unknown',position='Append')
      if(x.lt.-20) then
        if((y.lt.2*(40+x)).and.(y.gt.(-2)*(40+x))) then
! -------------------------------------------------------------------
         if(x.lt.-20*Alpha) then
          if((y.gt.2*(Alpha*40+x)).or.(y.lt.(-2)*(Alpha*40+x))) then
 	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
          else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
          endif
         endif
! -------------------------------------------------------------------
      	else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
      	endif
      elseif((x.ge.-20).and.(x.le.20)) then
! -------------------------------------------------------------------
       if(x.lt.-20*Alpha) then
        if((y.gt.2*(Alpha*40+x)).or.(y.lt.(-2)*(Alpha*40+x))) then
 	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
      	else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
      	endif
       elseif(x.gt.20*Alpha) then
        if((y.gt.2*(Alpha*40-x)).or.(y.lt.(-2)*(Alpha*40-x))) then
 	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
      	else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
      	endif
       else
        if(abs(y).gt.Alpha*40) then
 	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
      	else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
      	endif
       endif
! -------------------------------------------------------------------
      else
        if((y.lt.2*(40-x)).and.(y.gt.(-2)*(40-x))) then
! -------------------------------------------------------------------
         if(x.gt.20*Alpha) then
          if((y.gt.2*(Alpha*40-x)).or.(y.lt.(-2)*(Alpha*40-x))) then
 	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y
     &	    ,real(Uve_1),real(Uve_2),real(Uve_3)
          else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
          endif
         endif
! -------------------------------------------------------------------
      	else
      	    write(5,'(I3,I3,E11.4,E11.4,E11.4)') x,y,0.0,0.0,0.0
      	endif
      endif

      close(5)

      case default
c	No hace nada
      end select

      enddo
      enddo
      enddo
      enddo
      return
      end


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function �(x)
C       Input :  x  --- Argument of �(x)
C                       ( x is not equal to 0,-1,-2,���)
C       Output:  GA --- �(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END


        SUBROUTINE PSI(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
c-----------------------------------------------------------------------
