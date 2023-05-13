         subroutine NotAtALLNeeded ! just has the ARCCOM2e2.CH
c       ARCCOM2e2.CH =filename                                                    
        IMPLICIT REAL*8 (A-H,M,O-Z)                                    
        PARAMETER (NMX=200,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX, 
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2) 
         COMMON/DataForRoutines1/X(NMX3),V(NMX3),WTTL,M(NMX), 
     &   XC(NMX3),WC(NMX3),MC(NMX)  
     &  ,XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),N 
      COMMON/DataForChainRoutinesTwo/MMIJ,CMX(3),CMV(3)
     & ,ENERGY,Energr,CHTIME 
      common/softening/ee,cmethod(3),Clight,NofBH 
      common/TIMECOMMON/Taika,timecomparison              
      common/spincommon/spin(3)! the relative spin of M(1) !Spin=spin*G*M^2/c
      common/tolerancecommon/tolerance
            end

       function cdot(a,b)
       real*8  a(3),b(3),cdot
       cdot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
       return
       end


c----------------------------------------------------------
c    RELEVANT ROUTINES START HERE
c-----------------------------------------------------------
c        THIS IS THE MAIN SUBROUTINE.          
         SUBROUTINE chainevolve ! THIS ROUTINE CHANGED TO A STAND ALONE ONE!
     &  (NN,XX,VX,MX,TIME,DELTAT,TOL,NEWREG,KSMX,soft,cmet,cl,Ixc,NBH,
     &   spini,TDRadius)
         INCLUDE 'ARCCOM2e2.CH'
c         INCLUDE 'mergearc.inc'
         common/collision/icollision,ione,itwo,iwarning
         REAL*8 XX(*),VX(*),MX(*),cmet(3),spini(3)!,CMXX(3),CMVX(3)
         logical newreg
         save

         tnext=time+deltat    
         tstep=DELTAT           ! use a local tstep rather than DELTAT
         nmerger = 0        ! no mergers yet

 10      CONTINUE
         CALL  ARC
     &  (NN,XX,VX,MX,TIME,tstep,TOL,NEWREG,KSMX,soft,cmet,cl,Ixc,NBH,
     &  spini)
C     &   CMXX,CMVX) ! likewise not needed here

            step_now = step     ! save step 

         IF (icollision.NE.0) THEN ! handle a collison

            nmerger             = nmerger + 1

c	    write(6,*)' merger happened!'

            CALL  Merge_i1_i2   ! merge the two particles

            newreg=.TRUE.       ! chain has changed
            NN=N                ! copy new chain
            DO i=1,NN
               MX(i)=M(i)
               DO k=1,3
                  xx(3*i-3+k)=x(3*i-3+k)
                  vx(3*i-3+k)=v(3*i-3+k)
               ENDDO
            ENDDO               ! done copying new chain

            tstep=tnext-time    ! set time step to continue
                                ! chain integration
c            IF ((abs(tstep).GT.1.e-6*deltat).AND.(NN.GT.1)) GOTO 10

         ENDIF

         step=step_now          ! restore the earlier step 

         RETURN
         END

       Subroutine MERGE_I1_I2
        include 'ARCCOM2e2.CH'
        REAL*8 SM(NMX),XR(NMX3),XDR(NMX3)
         COMMON/collision/icollision,Ione,Itwo,iwarning
         common/outputindex/index4output(200),N_ini
         save
           L=0
           DO I=1,ione-1
           SM(I)=M(I)
           DO  K=1,3
           XR(3*I-3+K)=X(3*I-3+K)
           XDR(3*I-3+K)=V(3*I-3+K)
           end do 
           end do
                       Myks=M(ione)
                       Mkax=M(itwo)
          SM(Ione)=M(Ione)+M(Itwo)
          DO 6 K=1,3
          XR(3*Ione-3+K)=(M(Ione)*X((Ione-1)*3+K)
     &     +M(Itwo)*X((Itwo-1)*3+K))/SM(Ione)
          XDR(3*Ione-3+K)=(M(Ione)*V((Ione-1)*3+K)
     &     +M(Itwo)*V((Itwo-1)*3+K))/SM(Ione)
6         CONTINUE

          DO I=Ione+1,Itwo-1
          sm(i)=m(i)
          do k=1,3
          XR(3*I-3+K)=X(3*I-3+k)
          XDR(3*I-3+K)=V(3*I-3+k)
          end do
          end do
          
          do i=Itwo,N-1
c          index4output(i)=i+1+N_ini-N
          index4output(i)=index4output(i+1)
          end  do
          
          DO I=Itwo+1,N
          sm(i-1)=m(i)
          do k=1,3
          XR(3*I-6+K)=X(3*I-3+k)
          XDR(3*I-6+K)=V(3*I-3+k)
          end do
          end do
          
C         MOVE THE REDUCED SYSTEM TO M,X,V
          L=0
c         New value of the number of bodies.
          N=N-1
          if(Itwo.le.NofBH)NofBH=NofBH-1 ! # of BH's reduced!
          if(N.eq.1)then
          write(6,*)' Only one body left!'
          STOP
          end if


ccc          write(6,*)' M ',(m(k),k=1,N+1)
          DO 8 I=1,N
          M(I)=SM(I)
          DO 7 K=1,3
          X(3*i-3+k)=XR(3*i-3+k)
          V(3*i-3+k)=XDR(3*i-3+k)
7         CONTINUE
8         CONTINUE
          icollision=0
ccc          write(6,*),' m ',(M(k),k=1,N)
c             write(6,*)' Merge:',ione,itwo,'  N, NBH=',N,NofBH
c     &     ,' masses ',Myks,Mkax
                write(67,*)' merge ',ione,itwo
                    ione=0
                    itwo=0
          RETURN
          END

         SUBROUTINE ARC
     &  (NN,XX,VX,MX,TIME,DELTAT,TOL,NEWREG,KSMX,soft,cmet,cl,Ixc,NBH,
     &  spini)
C     &   ,CMXX,CMVX)     ! do not use these here
c        BETTER TO USE CM-coords & vels for XX & VX and CMXX CMVX
c        FOR CM-position (needed in the Perturbations routine).
c-----------------------------------------------------------------
c        NOTE: some variables (eg. Energy and EnerGR are only in the
c        common. The internal NB-energy = ENERGY+EnerGR  (should be)
c        Energy= integrated E-value (excluding grav.radiation)
c        EnerGr= Energy radiated away (grav.radiation if Clight.ne.0.0)
C        CHAIN INTEGRATION. Perturbations & CM-motion included (in principle).
c        NN=# of bodies; XX=(cm)coords, VX=(cm)vels, MX=masses,
cc        CMXX=coords of CM, CMVX=vels of CM ! removed
c        TIME=time, dettaT=time interval
c        STEP=stepsize (set=0 initially)
c        NEWREG=.true. iff chain membership has changed
c        KSMX=max # of steps without return (use some large # )
c        soft =optional softening( U=1/sqrt(r**2+soft**2) )
c        cmet= 3-d vector that determines the method:
c         (1,0,0) =logH, (0,1,0)=TTL,(0,0,1)=DIFSY2 without t-tranformation
c        
c        cl=speed of light 
c        NOTE: cl=0 => no relativistic terms !!!
c        Ixc = 1 => exact time, =0 no exact time but return after CHTIME>DELTAT

         INCLUDE 'ARCCOM2e2.CH'
         COMMON/DerOfTime/GTIME
         COMMON/DIAGNOSTICS/GAMMA,H,IWR
         common/omegacoefficients/OMEC(NMX,NMX)
         common/collision/icollision,ione,itwo,iwarning
         common/itemaxcommon/aitemax,itemax,itemax_used

         REAL*8 G0(3),XX(*),VX(*),MX(*),cmet(3),spini(3)!,CMXX(3),CMVX(3)
         REAL*8 Y(1500),SY(1500)
         LOGICAL MUSTSWITCH,NEWREG
         data ntrue,nfalse,nwritten/3*0/
         save
c          Initial constants of motion
           if(newreg)then
           ntrue=ntrue+1
           else
           nfalse=nfalse+1
           end if
           if(ntrue.gt.nfalse+10.and.nwritten.eq.0)then
           nwritten=1
           write(6,*)char(7),char(7)
           write(6,*)' NEWREG should be set .TRUE. only' 
           write(6,*)' in the very beginning of a new simulation'
           write(6,*)' NOT at every step!! (May reduce accuracy!!)'
           write(6,*)' even if it may look like the contrary.'
           end if
           if(NN.gt.NMX)then
           write(6,*)' THIS CODE CAN HANDLE ONLY ',NMX,' BODIES '
           write(6,*)' Yuo are trying to use N=',NN
           write(6,*)' Try increasing NMX in ARCCOM2e2.CH '
           write(6,*)' and increase some (large) dimensions elsewhere'
           write(6,*)' in the same proportion.  STOPPING'
           STOP
           end if
c           if(cmet(1).eq.0.0 .and. cmet(2).ne.0.0)then
c           write(6,*)' In this version cmethod(1) should not  be zero'
c           write(6,*)' if cmethod(2).ne.0.0 '
c           write(6,*)cmet,' = cmethod(k),k=1,3 '
c           write(6,*)' STOPPING '
c           STOP
c           end if
           if(deltat.eq.0.0 .and. Ixc .eq.1)then
           write(6,*)' You cannot use DELTA=0 and Ixc=1 '
           write(6,*)' since then every output will be at time=0 '
           write(6,*)' STOPPING '
           STOP
           end if
           if(cmet(1)+cmet(2)+cmet(3).eq.0.0)then
           write(6,*)' You have not defined the time-transformation'
           write(6,*)cmet,' = cmethod(k),k=1,3 '
           write(6,*)' STOPPING '
           STOP
           end if

           CHTIME=0.0
           icollision=0
           Taika=TIME ! to common
           NofBH=NBH  ! - " -
           IF(NEWREG)THEN
           iwarning=0
            itemax=12
            itemax_used=0
           ee=soft**2  ! to common
           do k=1,3
           spin(k)=spini(k) ! SPIN
           cmethod(k)=cmet(k) ! -"-
           end do
           clight=cl    ! -"-
           N=NN
           mass=0
           DO I=1,N
           M(I)=MX(I)
           mass=mass+m(i)
           END DO

           MMIJ=0.0
           DO I=1,N-1
           DO J=I+1,N
           MMIJ=MMIJ+M(I)*M(J)
           END DO
           END DO
           MMIJ=MMIJ/(N*(N-1)/2.d0)
           DO I=1,3*N
           X(I)=XX(I)
           V(I)=VX(I)
           END DO
           IF(MMIJ.eq.0.0)THEN
           write(6,*)'You have at most one non-zero mass => t''=1/0 and'
           write(6,*)'this does not work'
           STOP
           END IF
           call FIND CHAIN INDICES
           IF(IWR.GT.0)WRITE(6,1232)time,(INAME(KW),KW=1,N)
           call INITIALIZE XC and WC
           call CONSTANTS OF MOTION(ENERGY,G0,ALAG)
           EnerGr=0 ! energy radiated away
           gtime=1/ALAG 
           do K=1,3
c           CMX(K)=0
c           CMV(K)=0
           end do
           call omegacoef           
           STIME=0.0
           NEWREG=.FALSE.
            WTTL=Wfunction()
           mmss=0
           do i=1,n-1
           do j=i+1,n
           mmss=mmss+m(i)*m(j)
           end do 
           end do
            call Take Y from XC WC (Y,Nvar)
            do i=1,Nvar
            SY(i)=0
            end do
            call Initial Stepsize(X,V,M,N,ee,step) ! New initial step determination
            EPS=TOL
            Ncall=0
           END IF ! NEWREG
           KSTEPS=0
           nzero=0
777        KSTEPS=KSTEPS+1
           call Take Y from XC WC (Y,Nvar)
           call Obtain Order of Y(SY)
           if(Ixc.eq.1) step=min(step,1.1d0*deltaT/GTIME)
c           EPS=tolerance
                                         stime=0
                                    f1=chtime-deltaT ! for exact time iteration
                                         d1=gtime
111        hold=step 
c           write(6,*)' step=',step,Nvar,stepi
           stime=0
           Ncall=Ncall+1
           if(Ncall.eq.Ncall/100*100)then
           old=step
           call Initial Stepsize(X,V,M,N,ee,step_new) ! New initial step determination
           if(iwr.gt.0)write(6,*)' step=',step,' old=',
     &      old,' new/old=',step_new/old
            step=max(old,step_new)
           end if
          call  DIFSYAB(Nvar,EPS,SY,step,stime,Y)
c-------------------------------------try to change the step (if zero)
            if(step.eq.0.0)then
                 nzero=nzero+1
                 write(6,*)' STEP=0!',char(7)
              if(nzero.gt.10)then
                 write(6,*)' STOP, step permanently =0 '
                 STOP
                 else
                  if(nzero.ne.nzero/2*2)then
                  step=hold*2.d0  ! try increase
                  else
                  step=hold*.3d0   ! try decrease
                  end if
                write(6,*)' Trying again with h=',step
                goto 111
              end if
                 else
                 nzero=0
            end if
c--------------------------------------------------------

           call Put Y to XC WC  (Y,Nvar) 
            call CHECK SWITCHING CONDITIONS(MUST SWITCH)
            IF(MUST SWITCH)THEN
            call Chain Transformation !
            WTTL=Wfunction() ! this may not be necessary, but probably OK.
            call Take Y from XC WC(Y,Nvar)
              IF(IWR.GT.0) WRITE(6,1232)time+chtime,(INAME(KW),KW=1,N)

1232         FORMAT(1X,g12.4,' I-CHAIN',20I3)
             END IF ! MUST SWITCH
                                          f2=chtime-deltaT ! for exact time iteration 
                                          d2=gtime 
                                          x1=-stime
                                          x2=0.0
              DLT=DELTAT*(1-1.d-6)
      IF(CHTIME.LT.DLT.and.(KSTEPS.lt.KSMX)
     &  .and.(icollision.eq.0))goto 777
         if(KSTEPS.lt.KSMX .and.Ixc.eq.1)then
        ! ITERATE TO EXACT OUTPUTTIME
        call Iterate2ExactTime(Y,Nvar,deltaT,f1,d1,f2,d2,x1,x2)
         end if
                call update x and v 
                 DO I=1,3*N
                 XX(I)=X(I)
                 VX(I)=V(I)
                 END DO
                 DO I=1,3
                 spini(I)=spin(I)
c                 CMXX(I)=CMX(I)
c                 CMVX(I)=CMV(I)
                 END DO
                 TIME=TIME+CHTIME
                 RETURN
        END
         subroutine Iterate2ExactTime(Y0,Nvar,deltaT,f1,d1,f2,d2,x1,x2)
         INCLUDE 'ARCCOM2e2.CH'
         COMMON/DerOfTime/GTIME
         REAL*8 Y(1500),SY(1500),Y0(*)
         data tiny/1.d-6/
         save
           it=0 
           hs=abs(x1-x2)
c                   write(6,*)' |x1-x2| =',hs,x1,deltat,x2
 1111     CONTINUE
c          it=it+1
         do i=1,nvar
         y(i)=y0(i)
         end do
         stime=0
         dx1=-f1/d1
         dx2=-f2/d2
         if(abs(dx1).lt.abs(dx2))then
         xnew=x1+dx1
         else
         xnew=x2+dx2
         end if
c          
         test=(x1-xnew)*(xnew-x2)
         if(test.lt.(-tiny*hs).or.(it+1).eq.(it+1)/5*5)then
         xnew=(x1+x2)/2 ! bisect if out of interval
c           write(6,*)' bisection ',it, x1-xnew,xnew-x2,test
c     &,' it x1-t,t-x2 tst'
          end if

          sfinal=xnew

         call Put Y to XC WC  (Y,Nvar)
c--------------------------------------------------------------------------
                  call Obtain Order of Y(SY)
c                   eps=tolerance EPS in common
                   steppi=0
                   do k=1,5
                   step=sfinal-stime
                   if(abs(step).gt.1.e-6*abs(hs).or.k.eq.1)then !!!!
                   steppi=step
                   call  DIFSYAB(Nvar,EPS,SY,step,stime,Y)
                   it=it+1
                   else
                   goto 222
                   end if

                   end do
222               continue
                   call Put Y to XC WC  (Y,Nvar)
                 call UPDATE X AND V
                 fnew=chtime-deltaT
                 dfnew=gtime
c         keep it bracketed
           if(f1*fnew.le.0.0)then
            f2=fnew
            d2=dfnew
            x2=xnew
            else
            f1=fnew
            d1=dfnew
            x1=xnew
           end if
        if((abs(deltaT-chtime).gt.1.e-5*deltat).and.(it.lt.100))
     &     goto 1111
c ONE FINAL STEP SHOULD BE HERE (if above not-so-accurate test)          
c--------------------------------------------------------------------
           do i=1,Nvar
           y0(i)=y(i)
           end do
           call Put Y to XC WC  (Y,Nvar)
           call UPDATE X AND V
c           write(6,*)it,chtime-deltaT, ' ite T_error ',steppi
           return
           end

        subroutine LEAPFROG(STEP,Leaps,stime)
        implicit real*8 (a-h,M,o-z)
        save
        CALL PUT V 2 W
        hs=step
        h2=hs/2
        call XCmotion(h2)
        stime=stime+h2
        do k=1,Leaps-1
        call WCmotion(hs)
        call XCmotion(hs)
        stime=stime+hs
        end do
        call WCmotion(hs)
        call XCmotion(h2)
        stime=stime+h2
        return
        end

        function Wfunction()
        INCLUDE 'ARCCOM2e2.CH'
         common/omegacoefficients/OMEC(NMX,NMX)
         save
                 OMEGA=0.0d0
         DO I=1,N-1
         DO J=I+1,N
         if(omec(i,j).ne.0.0)then
         RIJ=SQRT(SQUARE(X(3*I-2),X(3*J-2)))
         OMEGA=OMEGA+omec(i,j)/RIJ 
         end if
         END DO
         END DO
         Wfunction=OMEGA
         RETURN
         END
         subroutine omegacoef
        INCLUDE 'ARCCOM2e2.CH'
        common/omegacoefficients/OMEC(NMX,NMX)
        SAVE
        icount=0
        do i=1,N-1
        do j=i+1,N
c        if(1.e-3*mmij.gt.m(i)*m(j).and.cmethod(2).ne.0.0)then 
         if(m(i)+m(j).gt.0.0 .and. cmethod(2).ne.0.0)then
          OMEC(I,J)=mmij
          OMEC(J,I)=mmij
          icount=icount+1
        else
          OMEC(I,J)=0
          OMEC(J,I)=0
        end if
        end do
        end do
        if(icount.eq.0.0)cmethod(2)=0 ! all terms zero anyway
        return
        end

        SUBROUTINE XCMOTION(hs)
        INCLUDE 'ARCCOM2e2.CH'

         COMMON/IncrementCommon/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        COMMON/DerOfTime/G
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
         save
        Te=-ENERGY-EnerGR
         if(cmethod(1).ne.0.0d0)then
        call EVALUATE V(V,WC)
        do I=1,N
        I0=3*I-3
        Te=Te+M(I)*(V(I0+1)**2+V(I0+2)**2+V(I0+3)**2)/2
        end do
         end if ! cmethod(1).ne.0.0d0
        G=1/(Te*cmethod(1)+WTTL*cmethod(2)+cmethod(3)) ! = t'
               if(G.lt.0.0.and.iwr.gt.0)then
               write(6,*)1/G,' tdot <0 ! '
        return ! seriously wrong, but may work (this step gets rejected)
               end if
        dT= hs*G
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        XCinc(L+K)=XCinc(L+k)+WC(L+K)*dT
        XC(L+K)=XC(L+K)+WC(L+K)*dT
        END DO
        END DO
        CHTIMEinc=CHTIMEinc+dT
        CHTIME=CHTIME+dT
        do k=1,3
        CMXinc(k)=CMXinc(k)+dt*cmv(k)
        cmx(k)=cmx(k)+dt*cmv(k)
        end do
        RETURN
        END

        subroutine PUT V 2 W
       include 'ARCCOM2e2.CH'
       common/vwcommon/Ww(nmx3),WTTLw,cmvw(3),spinw(3)
       save
       do i=1,3*(N-1)
       Ww(i)=WC(I) 
       end do
       WTTLw=WTTL
       do k=1,3
       spinw(k)=spin(k)
       cmvw(k)=cmv(k)
       end do
       return
       end

cccccccccccccccc
c     Modified 6/23/2009 by D. Merritt (name too long for compiler)
ccccccccccccccc
c       subroutine COORDINATE DEPENDENT PERTURBATIONS(ACC) ! USER DEFINED
       subroutine COORDINATE DEPENDENT PERTURBATI(ACC) ! USER DEFINED
        INCLUDE 'ARCCOM2e2.CH'

c    IN STAND-ALONE-MODE THIS ROUTINE SIMPLY SETS: ACC(I)=0

C---  added common blocks in include file
c        common/perturbers/MP(100),XP(300),NPertur ! this MUST be changed by user
c        INCLUDE 'chaininterface.inc'

        real*8 ACC(*)!,dx(3)!,xppred(3,NPMX)   ! NMAX comes from inc-file
        REAL*8 S,S1,S2
        save
c       HERE ONE MUST EVALUATE THE ACCELERATIONS DUE TO THE PERTURBERS.
C       Physical positions and velocities (in the inertial coordinate)
C       system are in vectors X and V
C       (X(1)=X_1,X(2)=Y_1,X(3)=Z_1, X(4)=X_2, X(5)=Y_2,...)
C       After a call to this routine the EXAccerations
C       are assumed to be in the vector ACC.


          TrueTIME=Taika+CHTIME ! this is (?) the time to be used to obtain perturber coords
c                               ! note: chtime is measured from the beginning of this CHAIN call
c           these times are in common
          
C---  init acc
        DO  I=1,3*N
         ki=i-(i-1)/3*3
          ACC(I)=0!-1.d-3*(x(i)+cmx(ki)) ! test only REMOVE this and REPLACE with whatever you want to
        END DO
         if(1.eq.2) then !THIS MUST BE REMOVED FOR REAL USE
c        This is an example only:
        do ic=1,N               ! ic= chainparticle index

           ic0=3*ic-3

           do i=1,Npertur       ! i=perturber index

              S  = TrueTIME - t0p(i)  ! this is the step size
              S1 = 0.5d0*S            ! and some multiples of it
              S2 = S/3.d0             ! as used in the predictor

              rr=0
              do k=1,3          ! coordinate index

C---  here comes the predictor
c                 xppred(k,i) =  ((fpdot(k,i)*S2 + fp(k,i))*S1
c     $                         +  vp(k,i))*S   + xp(k,i)

C---  and predicted distance to chain particle 
c                 dx(k)=X(ic0+k)-xppred(k,i)
c                 rr=rr+dx(k)**2
              end do            ! k
              rrr=sqrt((rr+ee)**3.d0)

              do k=1,3
c                 ACC(ic0+k)=ACC(ic0+k)-MP(I)/rrr*dx(k) ! hope this is correct (not confirmed)
              end do            ! k
               
           end do               ! i
        end do                  ! ic         

        end if !1.eq.2
        return
        end

        subroutine Velocity Dependent Perturbations
     &   (dT,Va,spina,acc,dcmv,df,dfGR,dspin)
        INCLUDE 'ARCCOM2e2.CH'
        real*8 df(*),Va(*),dcmv(3),dfGR(*),dfR(nmx3),acc(nmx3)
        real*8 dspin(3),spina(3)
        save
c       h\E4iri\F6t        
        
                 if(Clight.ne.0.0)then ! include only if Clight set >0
        call Relativistic ACCELERATIONS(dfr,dfGR,Va,spina,dspin)
              else
               do i=1,3*n
                dfr(i)=0
                 dfgr(i)=0
                 end do
                 do k=1,3
                 dspin(k)=0
                 end do
                 end if
        do i=1,3*n
         df(i)=acc(i)+dfr(i)
        end do
        call reduce 2 cm(df,m,n,dcmv) 
        return
        end 
        SUBROUTINE CHECK SWITCHING CONDITIONS(MUSTSWITCH)
        INCLUDE 'ARCCOM2e2.CH'
        LOGICAL MUSTSWITCH
        DATA Ncall,NSWITCH/0,200000000/
        save
        MUST SWITCH=.FALSE.
        Ncall=Ncall+1
C       Switch anyway after every NSWITCHth step.
        IF(Ncall.GE.NSWITCH)THEN
        Ncall=0
        MUST SWITCH=.TRUE.
        RETURN
        END IF
C       Inspect the structure of the chain.
C       NOTE: Inverse values 1/r are used instead of r itself.
        ADISTI=0.5*(N-1)/RSUM
        LRI=N-1
        DO I=1,N-2
        DO J=I+2,N
        LRI=LRI+1
C       Do not inspect if 1/r is small.
        IF(RINV(LRI).GT.ADISTI)THEN
         IF(J-I.GT.2)THEN
C        Check for a dangerous long loop.
C          RINVMX=MAX(RINV(I-1),RINV(I),RINV(J-1),RINV(J))
           IF(I.GT.1)THEN
           RINVMX=MAX(RINV(I-1),RINV(I))
           ELSE
           RINVMX=RINV(1)
           END IF
           RINVMX=MAX(RINVMX,RINV(J-1))
           IF(J.LT.N)RINVMX=MAX(RINVMX,RINV(J))
           IF(RINV(LRI).GT.RINVMX)THEN ! 0.7*RINVMX may be more careful
           MUST SWITCH=.TRUE.
           Ncall=0
           RETURN
           END IF
         ELSE
C        Is this a triangle with smallest size not regularised?
           IF( RINV(LRI).GT.MAX(RINV(I),RINV(I+1)))THEN
           MUST SWITCH=.TRUE.
           Ncall=0
           RETURN
           END IF
         END IF
        END IF
        END DO
        END DO
        RETURN
        END
         SUBROUTINE FIND CHAIN INDICES
         INCLUDE 'ARCCOM2e2.CH'
        REAL*8 RIJ2(NMXM)
        INTEGER IC(NMX2),IJ(NMXM,2),IND(NMXM)
        LOGICAL USED(NMXM),SUC,LOOP
        save
        L=0
        DO I=1,N-1
        DO J=I+1,N
        L=L+1
        RIJ2(L)=SQUARE(X(3*I-2),X(3*J-2))
        IJ(L,1)=I
        IJ(L,2)=J
        USED(L)=.FALSE.
        END DO
        END DO
        call ARRANGE(L,RIJ2,IND)
        LMIN=1+NMX
        LMAX=2+NMX
        IC(LMIN)=IJ(IND(1),1)
        IC(LMAX)=IJ(IND(1),2)
        USED(IND(1))=.TRUE.
1        DO I=2,L
        LI=IND(I)
        IF( .NOT.USED(LI))THEN
        call CHECK CONNECTION(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
        IF(SUC)THEN
        USED(LI)=.TRUE.
        GOTO 2
        ELSE
        USED(LI)=LOOP
        END IF
        END IF
        END DO
2        IF(LMAX-LMIN+1.LT.N)GO TO 1
        L=0
        DO I=LMIN,LMAX
        L=L+1
        INAME(L)=IC(I)
        END DO
        RETURN
        END
        SUBROUTINE CHECK CONNECTION(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
         INCLUDE 'ARCCOM2e2.CH'
        INTEGER IC(*),ICC(2),IJ(NMXM,2)
        LOGICAL SUC,LOOP
        save
        SUC=.FALSE.
        LOOP=.FALSE.
        ICC(1)=IC(LMIN)
        ICC(2)=IC(LMAX)
        DO I=1,2
        DO J=1,2
        IF(ICC(I).EQ.IJ(LI,J))THEN
        JC=3-J
        LOOP=.TRUE.
        DO L=LMIN,LMAX
        IF(IC(L).EQ.IJ(LI,JC))RETURN
        END DO
        SUC=.TRUE.
        LOOP=.FALSE.
        IF(I.EQ.1)THEN
        LMIN=LMIN-1
        IC(LMIN)=IJ(LI,JC)
        RETURN
        ELSE
        LMAX=LMAX+1
        IC(LMAX)=IJ(LI,JC)
        RETURN
        END IF
        END IF
        END DO
        END DO
        RETURN
        END
        SUBROUTINE ARRANGE(N,Array,Indx)
        implicit real*8 (a-h,o-z)
        dimension Array(*),Indx(*)
        save
        do 11 j=1,N
        Indx(j)=J
11      continue
        if(N.lt.2)RETURN
        l=N/2+1
        ir=N
10      CONTINUE
        IF(l.gt.1)THEN
        l=l-1
        Indxt=Indx(l)
        q=Array(Indxt)
        ELSE
        Indxt=Indx(ir)
        q=Array(Indxt)
        Indx(ir)=Indx(1)
        ir=ir-1
        IF(ir.eq.1)THEN
        Indx(1)=Indxt
        RETURN
        END IF
        END IF
        i=l
        j=l+l
20      IF(j.le.ir)THEN
            IF(j.lt.ir)THEN
               IF(Array(Indx(j)).lt.Array(Indx(j+1)))j=j+1
            END IF
            IF(q.lt.Array(Indx(j)))THEN
               Indx(i)=Indx(j)
               i=j
               j=j+j
            ELSE
               j=ir+1
            END IF
         GOTO 20
         END IF
         Indx(i)=Indxt
         GO TO 10
         END
        SUBROUTINE INITIALIZE XC AND WC
        INCLUDE 'ARCCOM2e2.CH'
        save
C        Center of mass
        DO K=1,3
        CMX(K)=0.0
        CMV(K)=0.0
        END DO
        MASS=0.0
        DO I=1,N
        L=3*(I-1)
        MC(I)=M(INAME(I)) ! masses along the chain
        MASS=MASS+MC(I)
         DO K=1,3
         CMX(K)=CMX(K)+M(I)*X(L+K)
         CMV(K)=CMV(K)+M(I)*V(L+K)
         END DO
        END DO
        DO K=1,3
        CMX(K)=CMX(K)/MASS
        CMV(K)=CMV(K)/MASS
        END DO
c       Rearange according to chain indices.
        DO I=1,N
        L=3*(I-1)
        LF=3*INAME(I)-3
         DO K=1,3
         XI(L+K)=X(LF+K)
         VI(L+K)=V(LF+K)
         END DO
        END DO

C       Chain coordinates & vels ! AND INITIAL `WTTL'
        WTTL=0            !  initialize W 
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        XC(L+K)=XI(L+K+3)-XI(L+K)
        WC(L+K)=VI(L+K+3)-VI(L+K)
        END DO
        END DO
        RETURN
        END
        SUBROUTINE UPDATE X AND V
        INCLUDE 'ARCCOM2e2.CH'
        REAL*8 X0(3),V0(3)
        save
C        Obtain physical variables from chain quantities.

        DO K=1,3
        XI(K)=0.0
        VI(k)=0.0
        X0(K)=0.0
        V0(k)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        VI(L+3+K)=VI(L+K)+WC(L+K)
        XI(L+3+K)=XI(L+K)+XC(L+K)
        END DO
        END DO
        DO I=1,N
        L=3*(I-1)
        DO K=1,3
        V0(K)=V0(K)+VI(L+K)*MC(I)/MASS
        X0(K)=X0(K)+XI(L+K)*MC(I)/MASS
        END DO
        END DO
C        Rearrange according to INAME(i) and add CM.
        DO I=1,N
        L=3*(I-1)
        LF=3*(INAME(I)-1)
        DO K=1,3
        X(LF+K)=XI(L+K)-X0(K)!+CMX(K) ! CM-coords
        V(LF+K)=VI(L+K)-V0(K)!+CMV(K) ! CM-vels
        END DO
        END DO
        RETURN
        END
        SUBROUTINE CHAIN TRANSFORMATION
        INCLUDE 'ARCCOM2e2.CH'
        REAL*8 XCNEW(NMX3),WCNEW(NMX3)
        INTEGER IOLD(NMX)
        save
        L2=3*(INAME(1)-1)
        DO K=1,3
        X(L2+K)=0.0
        END DO
C       Xs are needed when determining new chain indices.
        DO I=1,N-1
        L=3*(I-1)
        L1=L2
        L2=3*(INAME(I+1)-1)
        DO K=1,3
        X(L2+K)=X(L1+K)+XC(L+K)
        END DO
        END DO
C        Store the old chain indices.
        DO I=1,N
        IOLD(I)=INAME(I)
        END DO

C       Find new ones.
        call FIND CHAIN INDICES

C       Construct new chain coordinates. Transformation matrix
C       (from old to new) has only coefficients -1, 0 or +1.
        DO I=1,3*(N-1)
        XCNEW(I)=0.0
        WCNEW(I)=0.0
        END DO
        DO ICNEW=1,N-1
C       Obtain K0 &  K1 such that iold(k0)=iname(icnew)
c                                 iold(k1)=iname(icnew+1)
        LNEW=3*(ICNEW-1)
        DO I=1,N
        IF(IOLD(I).EQ.INAME(ICNEW))K0=I
        IF(IOLD(I).EQ.INAME(ICNEW+1))K1=I
        END DO
        DO ICOLD=1,N-1
        LOLD=3*(ICOLD-1)
        IF( (K1.GT.ICOLD).AND.(K0.LE.ICOLD))THEN
C       ADD
        DO K=1,3
        XCNEW(LNEW+K)=XCNEW(LNEW+K)+XC(LOLD+K)
        WCNEW(LNEW+K)=WCNEW(LNEW+K)+WC(LOLD+K)
        END DO
        ELSEIF( (K1.LE.ICOLD).AND.(K0.GT.ICOLD) )THEN
C        SUBTRACT
        DO K=1,3
        XCNEW(LNEW+K)=XCNEW(LNEW+K)-XC(LOLD+K)
        WCNEW(LNEW+K)=WCNEW(LNEW+K)-WC(LOLD+K)
        END DO
        END IF
        END DO
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,3*(N-1)   !!!!!!!!!!!!!!!!!!
        xc(i)=xcnew(i)   !!!!!!!!!!!!!!!!!!!
        wc(i)=wcnew(i)
        END DO           !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       Auxiliary quantities.
        MASS=0.0
        DO I=1,N
        MC(I)=M(INAME(I))
        MASS=MASS+MC(I)
        END DO

        RETURN
        END

        FUNCTION SQUARE(X,Y)
        implicit real*8 (a-h,m,o-z)
        REAL*8 X(3),Y(3),SQUARE
        common/softening/ee,cmethod(3),clight,NofBH ! only ee needed here
        save
        SQUARE=(X(1)-Y(1))**2+(X(2)-Y(2))**2+(X(3)-Y(3))**2+ee
        RETURN
        END
        
        SUBROUTINE DIFSYAB(N,EPS,S,h,t,Y)!,Jmax)
        implicit real*8 (a-h,o-z)
c       N=coordin. m\E4\E4r\81E(=3*NB)
c       F=funktion nimi (FORCE)
        parameter (NMX=1500,NMX2=2*NMX,nmx27=nmx2*7) ! NMX=MAX(N),N=3*NB
        REAL*8 Y(N),YR(NMX2),YS(NMX2),y0(NMX)
     +  ,DT(NMX2,7),D(7),S(N),EP(4)
        LOGICAL KONV,BO,KL,GR
        DATA EP/.4D-1,.16D-2,.64D-4,.256D-5/
        data dt/nmx27*0.0d0/
        save
        Jmax=10 ! JMAX set here
        IF(EPS.LT.1.D-14)EPS=1.D-14
        IF(N.gt.NMX)write(6,*) ' too many variables!', char(7)
        if(jmax.lt.4)write(6,*)' too small Jmax (=',jmax,')'
        JTI=0
        FY=1
        redu=0.8d0
        Odot7=0.7
        do i=1,N
        y0(i)=y(i)
        s(i)=max(abs(y0(i)),s(i))
        end do
10      tN=t+H
        BO=.FALSE.
C
        M=1
        JR=2
        JS=3
        DO  J=1,Jmax! 10

        do i=1,N
        ys(i)=y(i)
        s(i)=max(abs(ys(i)),s(i)) 
        end do
C

        IF(BO)then
        D(2)=1.777777777777778D0
        D(4)=7.111111111111111D0
        D(6)=2.844444444444444D1
        else
        D(2)=2.25D0
        D(4)=9.D0
        D(6)=36.0D0
        end if

        IF(J.gt.7)then
        L=7
        D(7)=6.4D1
        else
        L=J
        D(L)=M*M
        end if

        KONV=L.GT.3
           subH=H/M
           call SubSteps(Y0,YS,subH,M) ! M substeps of size H/M.
        KL=L.LT.2
        GR=L.GT.5
        FS=0.



        DO  I=1,N 
        V=DT(I,1)
        C=YS(I)
        DT(I,1)=C
        TA=C

        IF(.NOT.KL)THEN
        DO  K=2,L
        B1=D(K)*V
        B=B1-C
        W=C-V
        U=V
        if(B.ne.0.0)then
        B=W/B
        U=C*B
        C=B1*B
        end if
        V=DT(I,K)
        DT(I,K)=U
        TA=U+TA
        END DO ! K=2,L
        SI=max(S(I),abs(TA))
        IF(DABS(YR(I)-TA).GT.SI*EPS) KONV=.FALSE.
        IF(.NOT.(GR.OR.SI.EQ.0.D0))THEN
        FV=DABS(W)/SI
        IF(FS.LT.FV)FS=FV
        END IF
        END IF ! .NOT.KL.
        YR(I)=TA
        END DO ! I=1,N

c       end of I-loop
        IF(FS.NE.0.D0)THEN
        FA=FY
        K=L-1
        FY=(EP(K)/FS)**(1.d0/FLOAT(L+K))
ccccccccccccccccccc
c     Modified 6/23/2009 by D. Merritt
ccccccccccccccccccc
c        FY=min(FY,1.4) !1.4 ~ 1/0.7 ; where 0.7 = initial reduction factor
        FY=min(FY,1.4d0) !1.4 ~ 1/0.7 ; where 0.7 = initial reduction factor
        IF(.NOT.((L.NE.2.AND.FY.LT.Odot7*FA).OR.FY.GT.Odot7))THEN
        H=H*FY
               JTI=JTI+1
               IF(JTI.GT.25)THEN
               H=0.0
               RETURN
               END IF
        GO TO 10 ! Try again with a smaller step.
        END IF
        END IF

        IF(KONV)THEN
        t=tN
        H=H*FY
        DO  I=1,N
        Y(I)=YR(I)+y0(i) !!!!!!!
        END DO
        RETURN
        END IF

        D(3)=4.D0
        D(5)=1.6D1
        BO=.NOT.BO
        M=JR
        JR=JS
        JS=M+M
        END DO ! J=1,Jmax
        redu=redu*redu+.001d0 ! square the reduction factor (but minimum near 0.001)
        H=H*redu 
        GO TO 10 ! Try again with smaller step.
        END

        subroutine SubSteps(Y0,Y,H,Leaps)
        implicit real*8 (a-h,m,o-z)
        real*8 Y(*),Y0(*)!,ytest(1000)
        common/softening/ee,cmethod(3),Clight,NofBh
        COMMON/collision/icollision,Ione,Itwo,iwarning
        save
        icollision=0
        call Put Y to XC WC  (Y0,Nvar) ! Y -> XC, WTTL, WC
        call Initialize increments 2 zero
        call  LEAPFROG(H,Leaps,stime) ! advance 
c        call Take Y from XC WC (Ytest,Nvars) ! take new results -> Y
        call take increments 2 Y(y)
c        do i=1,Nvar
c        write(6,*)i,Y0(i)+Y(i)-Ytest(I),' i yNew-Yold '
c        end do
        return
        end
        subroutine Initialize increments 2 zero
        INCLUDE 'ARCCOM2e2.CH'
        COMMON/IncrementCommon/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        do i=1,3*(N-1)
        XCinc(i)=0
        WCinc(i)=0
        end do
        do k=1,3
        CMXinc(k)=0
        CMVinc(k)=0
        spin inc(k)=0
        end do
        WTTLinc=0
        ENERGYinc=0
        EnerGRinc=0
        CHTIMEinc=0
        return
        end

        subroutine Take Increments 2 Y(Y)
        INCLUDE 'ARCCOM2e2.CH'
        COMMON/IncrementCommon/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)

        real*8 Y(*)
        save
        L=1
        Y(L)=CHTIMEinc
        do i=1,3*(N-1)
        L=L+1   
        Y(L)=XCinc(I)
        end do
        L=L+1
        Y(L)=WTTLinc
        do i=1,3*(N-1)
        L=L+1
        Y(L)=WCinc(I)
        end do
        do i=1,3
        L=L+1
        Y(L)=CMXinc(I)
        end do
        do i=1,3
        L=L+1
        Y(L)=CMVinc(I)
        end do
        L=L+1
        Y(L)=ENERGYinc
        L=L+1
        Y(L)=EnerGRinc
        do k=1,3
        L=L+1
        Y(L)=spin inc(k)
        end do
c        Nvar=L  
        RETURN
        END

        subroutine Put Y to XC WC (Y,Lmx)
         INCLUDE 'ARCCOM2e2.CH'
        real*8 Y(*)
        save
        L=1
        CHTIME=Y(L)
        do i=1,3*(N-1)
        L=L+1
        XC(I)=Y(L)
        end do
        L=L+1
        WTTL=Y(L)
        do i=1,3*(N-1)
        L=L+1
        WC(I)=Y(L)
        end do
        do i=1,3
        L=L+1
        CMX(I)=Y(L)
        end do
        do i=1,3
        L=L+1
        CMV(I)=Y(L)
        end do
        L=L+1
        ENERGY=Y(L)
        L=L+1
        EnerGR=Y(L)
        do k=1,3
        L=L+1
        spin(k)=Y(L)
        end do
        Lmx=L
        RETURN
        END
        subroutine Take Y from XC WC (Y,Nvar)
         INCLUDE 'ARCCOM2e2.CH'
        real*8 Y(*)
        save
        L=1
        Y(L)=CHTIME
        do i=1,3*(N-1)
        L=L+1   
        Y(L)=XC(I)
        end do
        L=L+1
        Y(L)=WTTL
        do i=1,3*(N-1)
        L=L+1
        Y(L)=WC(I)
        end do
        do i=1,3
        L=L+1
        Y(L)=CMX(I)
        end do
        do i=1,3
        L=L+1
        Y(L)=CMV(I)
        end do
        L=L+1
        Y(L)=ENERGY
        L=L+1
        Y(L)=EnerGR
        do k=1,3
        L=L+1
        Y(L)=spin(k)
        end do
        Nvar=L  
        RETURN
        END
        subroutine Obtain Order of Y(SY)
        INCLUDE 'ARCCOM2e2.CH'
        real*8 SY(*)
        save
        w_old=0.1d0
        w_new=1-w_old
        L=1
        SY(L)=ABS(CHTIME)*w_new+sy(L)*w_old
        SR=0
        XCmin=1.d99
        UPO=0
        do i=1,N-1
        i0=3*i-3
        XCA=abs(XC(I0+1))+abs(XC(I0+2))+abs(XC(I0+3))
        SR=SR+XCA
        UPO=UPO+MMIJ/XCA
        XCmin=min(XCA,XCmin)
         do k=1,3 
         L=L+1
         SY(L)=XCA*w_new+sy(L)*w_old
         end do ! k
        end do  ! I
        L=L+1
        SY(L)=(abs(WTTL*1.e2)+mass**2/XCmin)*w_new+sy(L)*w_old
        SW0=sqrt(abs(Energy/mass))
        SW=0
        do i=1,N-1
        i0=3*i-3
        WCA=abs(WC(I0+1))+abs(WC(I0+2))+abs(WC(I0+3))
        SW=SW+WCA
        do k=1,3
        L=L+1

        if(WCA.ne.0.0)then
        SY(L)=WCA*w_new+sy(L)*w_old
        else
        SY(L)=SW0*w_new+sy(L)*w_old
        end if
        end do ! k
        end do ! i

        L=1
        do i=1,N-1
        i0=3*i-3
         do k=1,3 
         L=L+1
         if(SY(L).eq.0.0)SY(L)=SR/N*w_new+sy(L)*w_old
         end do ! k
        end do  ! I
        L=L+1 ! WTTL
        do i=1,N-1
        i0=3*i-3
        do k=1,3
        L=L+1
        if(SY(L).eq.0.0)SY(L)=(SW/N+sqrt(UPO/mass))*w_new+sy(L)*w_old
c        if(SY(L).eq.0.0)SY(L)=1
        end do ! k
        end do ! i


        CMXA=abs(cmx(1))+abs(cmx(2))+abs(cmx(3))+1.d-2*SR/N
        CMVA=abs(cmv(1))+abs(cmv(2))+abs(cmv(3))+1.d-2*SW/N

        do i=1,3
        L=L+1
        SY(L)=CMXA*w_new+sy(L)*w_old ! cmx
        end do

        do i=1,3
        L=L+1
        SY(L)=CMVA*w_new+sy(L)*w_old ! cmv
        end do

        L=L+1
        SY(L)=(ABS(ENERGY)+0.1*UPO)*w_new+sy(L)*w_old ! E
        L=L+1
        SY(L)=SY(L-1)*w_new+sy(L)*w_old
        if(SY(1).eq.0.0)SY(1)=(sqrt(sr/mass)*sr*1.d-2)*w_new+sy(1)*w_old ! time
        do k=1,3
        L=L+1
        SY(L)=1. ! spin comps. 
        end do

        RETURN
        END

        SUBROUTINE EVALUATE X
        INCLUDE 'ARCCOM2e2.CH'
        REAL*8 X0(3)
        save
C        Obtain physical variables from chain quantities.

        DO K=1,3
        XI(K)=0.0
        X0(K)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        XI(L+3+K)=XI(L+K)+XC(L+K)
        END DO
        END DO
        DO I=1,N
        L=3*(I-1)
        DO K=1,3
        X0(K)=X0(K)+XI(L+K)*MC(I)/MASS
        END DO
        END DO
C        Rearrange according to INAME(i) and add CM.
        DO I=1,N
        L=3*(I-1)
        LF=3*(INAME(I)-1)
        DO K=1,3
        X(LF+K)=XI(L+K)-X0(K)!+CMX(K) ! CM-coords
        END DO
        END DO
        RETURN
        END
        SUBROUTINE EVALUATE V(VN,WI)
        INCLUDE 'ARCCOM2e2.CH'
        REAL*8 V0(3),VN(*),WI(*)
        save
C        Obtain physical V's from chain quantities.

        DO K=1,3
        V0(k)=0.0
        VI(k)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        VI(L+3+K)=VI(L+K)+WI(L+K)!WC(L+K)
        END DO
        END DO
        DO I=1,N
        L=3*(I-1)
        DO K=1,3
        V0(K)=V0(K)+VI(L+K)*MC(I)/MASS
        END DO
        END DO
C        Rearrange according to INAME(i) and add CM.
        DO I=1,N
        L=3*(I-1)
        LF=3*(INAME(I)-1)
        DO K=1,3
        VN(LF+K)=VI(L+K)-V0(K)!+CMV(K)
        V(LF+K)=VN(LF+K) ! 
        END DO
        END DO
        RETURN
        END

       SUBROUTINE Relativistic ACCELERATIONS(ACC,ACCGR,Va,spina,dspin)
        INCLUDE 'ARCCOM2e2.CH'
        REAL*8 ACC(*),dX(3),dW(3),dF(3),Va(*),ACCGR(*),dfGR(3),dsp(3)
     &  ,spina(3),dspin(3)
          common/collision/icollision,ione,itwo,iwarning
        common/notneeded/rijnotneeded
                 common/deeveet/dv2(3),dv4(3),dv5(3)
        save
        Cl=Clight! SPEED OF LIGHT 
C       INITIALIZE THE relativistic acceration(s) here. 
        DO  I=1,3*N
        ACC(I)=0.0
        ACCGR(I)=0.0
        end do
        do k=1,3
        dspin(k)=0
        end do
        DO IK=1,N
        I=INAME(IK)
        I3=3*I
        I2=I3-1
        I1=I3-2
        DO  JK=IK+1,N
        J=INAME(JK)
        IF(min(i,j).le.NofBH)THEN  ! only BH - BH, max->min => BH*
c	IF(1.eq.1)THEN  ! only BH - BH, max->min => BH*
        J3=J+J+J
        J2=J3-1
        J1=J3-2
        if(JK.NE.IK+1)THEN
        dx(1)=X(J1)-X(I1)
        dx(2)=X(J2)-X(I2)
        dx(3)=X(J3)-X(I3)
        dw(1)=Va(J1)-Va(I1)
        dw(2)=Va(J2)-Va(I2)
        dw(3)=Va(J3)-Va(I3)
        ELSE
        K1=3*IK-2
        K2=K1+1
        K3=K2+1
        dx(1)=XC(K1)
        dx(2)=XC(K2)
        dx(3)=XC(K3)
        dw(1)=Va(J1)-Va(I1) 
        dw(2)=Va(J2)-Va(I2)
        dw(3)=Va(J3)-Va(I3)
        END IF
        vij2=dw(1)**2+dw(2)**2+dw(3)**2

c       This (cheating) avoids vij>cl and produces only O(1/c^6) 'errors'.
        do k=1,3
        dw(k)=dw(k)/(1+(vij2/cl**2)**4)**.125d0 !  avoid V_ij > c !!
c        dw(k)=dw(k)/(1+(vij2/cl**2)**2)**.25d0 ! not so good
        end do
        
        vij2=dw(1)**2+dw(2)**2+dw(3)**2
        RS=2.d0*(m(i)+m(j))/CL**2

        RIJ2=dx(1)**2+dx(2)**2+dx(3)**2
        rij=sqrt(rij2)
        rdotv=dx(1)*dw(1)+dx(2)*dw(2)+dx(3)*dw(3)
        Ii=min(i,j)
        Jx=max(i,j)
         call Relativistic
     & Terms(Ii,dX,dW,rij,rdotv,vij2,m(Ii),m(Jx),cl,DF,dfGR,spina,dsp)
            RS=2.d0*(m(i)+m(j))/CL**2

c      if(rij.lt.4.*RS.and.iwarning.lt.2)
!       IF(RIJ.LT.TDRADIUS.AND.IWARNING.LT.2)
!     &  WRITE(6,*)' NEAR COLLISION: R/RS',RIJ/RS,I,J
!     & ,RIJ,RS ! diagno
c            if(rij.lt.4*RS)then!
c       if(rij.lt.TDRadius)then!
	if((rij.lt.0.1).and.((i.eq.1).or.(j.eq.1)))then
            iwarning=iwarning+1
            icollision=1   ! collision indicator
            ione=min(i,j)
            itwo=max(i,j)
            return
            end if
 	IF((RIJ.LT.0.046).AND.((I.EQ.2).OR.(J.EQ.2)))THEN
 			IWARNING=IWARNING+1
 			ICOLLISION=1   ! COLLISION INDICATOR
 			IONE=MIN(I,J)
 			ITWO=MAX(I,J)
 			RETURN
 			END IF
         do k=1,3
         dspin(k)=dspin(k)+dsp(k)
         end do
        ACC(I1)=ACC(I1)+m(j)*dF(1) ! here I assume action = reaction
        ACC(I2)=ACC(I2)+m(j)*dF(2) ! which is not really true for 
        ACC(I3)=ACC(I3)+m(j)*dF(3) ! relativistic terms (but who cares)
        ACC(J1)=ACC(J1)-m(i)*dF(1)
        ACC(J2)=ACC(J2)-m(i)*dF(2)
        ACC(J3)=ACC(J3)-m(i)*dF(3)
c        Grav.Rad.-terms
        ACCgr(I1)=ACCgr(I1)+m(j)*dFgr(1) ! here I assume action = reaction
        ACCgr(I2)=ACCgr(I2)+m(j)*dFgr(2) ! which is not really true for 
        ACCgr(I3)=ACCgr(I3)+m(j)*dFgr(3) ! relativistic terms (but who cares)
        ACCgr(J1)=ACCgr(J1)-m(i)*dFgr(1)
        ACCgr(J2)=ACCgr(J2)-m(i)*dFgr(2)
        ACCgr(J3)=ACCgr(J3)-m(i)*dFgr(3)

                      END IF
        END DO ! J
        END DO ! I
        
        RETURN
        END

         subroutine 
     & Relativistic terms(I1,X,V,r,rdotv,v2,m1,m2,c,DV,DVgr,spina,dspin)
         implicit real*8 (a-h,m,n,o-z)
         real*8 X(3),V(3),DV(3),n(3),ny,nv,m1,m2,m
         real*8 dv2(3),dv3(3),dv4(3),dv5(3),dvgr(3),spina(3),dspin(3)
         data beta,gamma/1.d0,1.d0/
         save
          m=m1+m2

          my=m1*m2/m

          ny=my/m
          n(1)=x(1)/r
          n(2)=x(2)/r
          n(3)=x(3)/r
 
          nv=rdotv/r
          v4=v2*v2
           r2=r*r
           
                         IF(1.eq.1)THEN
           do i=1,3
          dv2(i)=m/c**2*n(i)/r2*(m/r*(2*(beta+gamma)+2*ny) ! 1/c**2 terms
     &    -v2*(gamma+3*ny)+3*ny/2*nv**2)
     &    +m*v(i)*nv/c**2/r**2*(2*gamma+2-2*ny)
c           dv2(i)=0.d0
           end do
         
          do i=1,3
          dv4(i)=1/c**4*(                               ! 1/c**4 terms
     & +ny*m*n(i)/r2*(-2*v4+1.5d0*v2*nv**2*(3-4*ny)-15*nv**4/8*(1-3*ny))
     & +m**2*n(i)/r**3*(v2/2*ny*(11+4*ny)+2*nv**2*(1+ny*(12+3*ny)))
     & +ny*m*v(i)/r**2*(8*v2*nv-3*nv**3/2*(3+2*ny))
     & -m**2/2/r**3*v(i)*nv*(4+43*ny)-m**3*n(i)/r**4*(9+87*ny/4) )
c            dv4(i) = 0.d0
           end do
                  ELSE
                  do k=1,3
                  dv2(k)=0
                  dv4(k)=0
                  end do
                  END IF
           do i=1,3                         
           dv5(i)=ny/c**5*(   ! gravitational radiation terms
     &    -8*m**2/r**3/5*(v(i)*(v2+3*m/r)-n(i)*nv*(3*v2+17*m/3/r)))
              dv5(i) = 0.d0
           end do
         IF(I1.eq.1)then
         call gopu_SpinTerms(X,V,r,M1,m2,c,spina,dv3,dspin) ! spinterms ->dv3
         else
         do k=1,3
         dv3(k)=0
         dspin(k)=0
         end do
         end if
           do i=1,3
           dv(i)=-1/m*(dv2(i)+dv3(i)+dv4(i)+dv5(i))
           dvgr(i)=-1/m*dv5(i)
           end do
          return
          end

        subroutine Reduce2cm(x,m,nb,cm)
        implicit real*8 (a-h,m,o-z)
        real*8 x(*),m(*),cm(3)
        save
        cm(1)=0
        cm(2)=0
        cm(3)=0
        sm=0
        do i=1,nb
        sm=sm+m(i)
        do k=1,3
        cm(k)=cm(k)+m(i)*x(k+3*(i-1))
        end do
        end do
        do k=1,3
        cm(k)=cm(k)/sm
        end do
        do i=1,nb
        do k=1,3
        x(k+3*(i-1))=x(k+3*(i-1))-cm(k)
        end do
        end do
        return
        end

       subroutine cross(a,b,c)
       real*8 a(3),b(3),c(3)
       save
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end


      subroutine gopu_SpinTerms(X,V,r,M1,m2,c,alpha,dv3,dalpha)
      implicit real*8 (a-h,m,n,o-z)
      real*8 x(3),v(3),dv3(3),n(3)
      real*8 dalpha(3),w(3),alpha(3)
      real*8 nxa(3),vxa(3),J(3)
      save
                   ! This routine assumes: The BH mass M1>>m2. Spin of
                   ! m2 is neglected.
       do k=1,3
       n(k)=x(k)/r 
       end do
       m=m1+m2
       eta=m1*m2/m**2
       SQ=sqrt(1-4*eta)
       Aq=-12/(1+sq)
       Bq= -6/(1+sq)-3
       Cq=1+6/(1+sq)
       rdot=cdot(n,v)
       call cross(n,v,w)
       anxv=cdot(alpha,w)
       call cross(n,alpha,nxa)
       call cross(v,alpha,vxa)
       do k=1,3
       dv3(k)=-m1**2/(c*r)**3*
     & (Aq*anxv*n(k)+rdot*Bq*nxa(k)+Cq*vxa(k))
       end do
       coeff=eta*m/(c*r)**2*(3/(1+sq)+.5d0)
       call cross(w,alpha,dalpha)
        do k=1,3
        dalpha(k)=coeff*dalpha(k)
        end do
c  C.Will Q2-terms
        sjj=0
        do k=1,3
        j(k)=M1**2/c*alpha(k)
        sjj=sjj+j(k)**2
        end do
        sj=sqrt(sjj)
        if(sj.ne.0.0)then  ! if sj=0, then J(k)=0 and Q-term =0 anyway
        do k=1,3
        j(k)=j(k)/sj
        end do
        end if
        Q2=-sjj/M1/c**2 !  X=X_j-X_i in this code
        do k=1,3
       dv3(k)=dv3(k)  ! add Quadrupole terms
     &  +1.5*Q2/r**4*(n(k)*(5*cdot(n,j)**2-1)-2*cdot(n,j)*j(k))
        end do
       return
       end

        SUBROUTINE Initial Stepsize(X,V,M,NB,ee,step)
        IMPLICIT real*8 (A-H,m,O-Z)
        DIMENSION X(*),V(*),M(*)
        save
        T=0.0
        U=0.0
        RMIN=1.D30
        mass=M(NB)
        time_step2=1.e30
        DO I=1,NB-1
        mass=mass+M(I)
        DO J=I+1,Nb
        MIJ=M(I)*M(J)
        KI=(I-1)*3
        KJ=(J-1)*3
        xx=X(KI+1)-X(KJ+1)
        yy=X(KI+2)-X(KJ+2)
        zz=X(KI+3)-X(KJ+3)
        R2=xx*xx+yy*yy+zz*zz+ee
        vx=V(KI+1)-V(KJ+1)
        vy=V(KI+2)-V(KJ+2)
        vz=V(KI+3)-V(KJ+3)
        vv=vx*vx+vy*vy+vz*vz
        R1=Sqrt(R2)
        time_step2=min(time_step2,R2/(.5*vv+(M(I)+M(J))/R1)) ! ~2B radius of convergence^2
        U=U+MIJ/R1
        T=T+MIJ*(vx*vx+vy*vy+vz*vz)
        END DO
        END DO
        T=T/(2*mass)
        ENERGY=T-U
        Alag=T+U
        STEP=0.5*U*sqrt(time_step2)        
        RETURN
        END

        subroutine elmnts
     & (x,v,m,a,e,mo,inc,Om,oo,alfa,q,tq)
c       NOTE: wrong results can be produced in exeptional situations
c       where some angles are undefined in terms of the expressions used.
c       This may happen in exactly planar, rectilinear .. orbits
c       Troubles can often be avoided by a very small 'perturbation' of x and/or v.
        implicit real*8 (a-h,m,o-z)
        parameter(rad=180.d0/3.141592653589793d0 )
        real*8 x(3),w(3),v(3),inc,jx,jy,jz
        save
        mu=sqrt(m)
        do k=1,3
        w(k)=v(k)/mu
        end do
        r=sqrt(x(1)**2+x(2)**2+x(3)**2)
        w2=w(1)**2+w(2)**2+w(3)**2
        eta=x(1)*w(1)+x(2)*w(2)+x(3)*w(3)
        alfa=2/r-w2
        zeta=1-alfa*r

c       areal velocity vector (jx,jy,jz)
        jx=x(2)*w(3)-x(3)*w(2)
        jy=x(3)*w(1)-x(1)*w(3)
        jz=x(1)*w(2)-x(2)*w(1)
        d=sqrt(jx*jx+jy*jy+jz*jz)

c       eccentricity vector (ex,ey,ez)
        ex=w(2)*jz-w(3)*jy-x(1)/r
        ey=w(3)*jx-w(1)*jz-x(2)/r
        ez=w(1)*jy-w(2)*jx-x(3)/r

        e=sqrt(ex*ex+ey*ey+ez*ez)
        b=sqrt(jx*jx+jy*jy)
        inc=atn2(b,jz)*rad
        Om=atn2(jx,-jy)*rad
        oo=atn2(ez*d,ey*jx-ex*jy)*rad
        a=1/alfa
        sqaf=sqrt(abs(alfa))
        q=d*d/(1+e) 
        too=oot(alfa,eta,zeta,q,e,sqaf)
        tq=too/mu
        mo=too*sqaf**3*rad
        return
        end
        function atn2(s,c)
        implicit real*8 (a-h,o-z)
        parameter(twopi=2*3.141592653589793d0)
        save
        atn2=atan2(s,c)
        if(atn2.lt.0.0)atn2=atn2+twopi
        return
        end
        function oot(alfa,eta,zeta,q,e,sqaf) ! oot=pericentre time
c       alfa=1/a; eta=sqrt(a) e sin(E); zeta=e Cos(E),
c       q=a(1-e), e=ecc, sqaf=sqrt(|a|)
        implicit real*8 (a-h,o-z)
        parameter(tiny=1.d-18)
        save
        if(zeta.gt.0.0)then 
c        ellipse (near peri), parabola or hyperbola.
         ecc=max(e,tiny)
         X=eta/ecc
         Z=alfa*X*X
         oot=X*(q+X*X*g3(Z))
        else
c       upper half of an elliptic orbit.
        oot=(atan2(eta*sqaf,zeta)/sqaf-eta)/alfa
        end if
         return
         end
        function g3(z)
        implicit real*8 (a-h,o-z)
        common/mita/zero
        save
        if(z.gt.0.025d0)then ! elliptic
        x=sqrt(z)
        g3 = (asin(x)-x)/x**3
        elseif(z.lt.-0.025d0)then ! hyperbolic
        x = sqrt(-z)
        g3 = (log(x+sqrt(1+x*x))-x )/x/z
        else ! Pade approximant for small  |z|
c       g3 = (1/6.d0-19177*z/170280 + 939109*z*z/214552800)/
c     &  (1-7987*z/7095 + 54145*z*z/204336)
       g3 = (1+6*(-19177*z/170280 + 939109*z*z/214552800))/
     &  (6*(1-7987*z/7095 + 54145*z*z/204336))
        zero=0
        end if
        return
        end

        SUBROUTINE CONSTANTS OF MOTION(ENE_NB,G,Alag)
c        IMPLICIT real*8 (A-H,m,O-Z)
c        DIMENSION G(3)
        include 'ARCCOM2e2.CH'
         real*8 g(3)
        common/justforfun/Tkin,Upot,dSkin,dSpot
        save
c       Contants of motion in the centre-of-mass system.        
        T=0.0
        U=0.0
        G(1)=0.
        G(2)=0.
        G(3)=0.
        RMIN=1.D30
c        mass=M(N)
        DO Ik=1,N-1
        I=INAME(IK)      ! along the chain
c        mass=mass+M(I)
        DO Jk=Ik+1,N
        J=INAME(JK)      !  -"-
        MIJ=M(I)*M(J)
        KI=(I-1)*3
        KJ=(J-1)*3
        IF(JK.NE.IK+1)THEN
        xx=X(KI+1)-X(KJ+1)
        yy=X(KI+2)-X(KJ+2)
        zz=X(KI+3)-X(KJ+3)
        vx=V(KI+1)-V(KJ+1)
        vy=V(KI+2)-V(KJ+2)
        vz=V(KI+3)-V(KJ+3)
        ELSE
        K1=3*IK-2
        K2=K1+1
        K3=K2+1
        XX=XC(K1)   ! use chain vectors when possible
        YY=XC(K2)   ! (this often reduces roundoff)
        ZZ=XC(K3)
        VX=WC(K1)
        VY=WC(K2)
        VZ=WC(K3)
        END IF

        R2=xx*xx+yy*yy+zz*zz+ee

        U=U+MIJ/SQRT(R2)
        T=T+MIJ*(vx*vx+vy*vy+vz*vz)
        G(1)=G(1)+MIJ*(yy*vz-zz*vy)
        G(2)=G(2)+MIJ*(zz*vx-xx*vz)
        G(3)=G(3)+MIJ*(xx*vy-yy*vx)
        END DO
        END DO
        T=T/(2*mass)
        G(1)=G(1)/mass
        G(2)=G(2)/mass
        G(3)=G(3)/mass
        ENE_NB=T-U
        Alag=T+U
        Tkin=T ! to justforfun
        Upot=U ! to justforfun
        OmegaB=Wfunction()
        dSkin=cmethod(1)*(T-ENERGY-ENERGR)+cmethod(2)*WTTL+cmethod(3)
        dSpot=cmethod(1)*U+cmethod(2)*OmegaB+cmethod(3)
        RETURN
        END

      SUBROUTINE FIND BINARIES(time)  ! this is a toy analysis routine
      INCLUDE 'ARCCOM2e2.CH'
      REAL*8 XX(3),W(3)
C     SEARCH FOR BINARIES [diagnostics only]
      save 
      DO I=1,N-1
         DO J=I+1,N
            LI=3*(I-1)
            LJ=3*(J-1)
            OM=1./SQRT(M(I)+M(J))
            DO K=1,3
               XX(K)=X(LI+K)-X(LJ+K)
               W(K) =(V(LI+K)-V(LJ+K))*OM
            END DO
            R2=XX(1)**2+XX(2)**2+XX(3)**2
            ETA=XX(1)*W(1)+XX(2)*W(2)+XX(3)*W(3)
            W2=W(1)**2+W(2)**2+W(3)**2
            R=SQRT(R2)
            OA=2./R-W2
            ZETA=1.-OA*R
            ECC2=ZETA**2+OA*ETA**2
            ECC=SQRT(ECC2)
            OA0=2.*(N-1)/(RSUM+1.E-20)
            IF(OA.GT.OA0)THEN
               WRITE(88,123)time,I,J,1./OA,ECC
               call flush(88)
            END IF
         END DO
      END DO
123   FORMAT
     &(1x,F12.1,' BINARY:(',I3,',',I3,')'
     & ,' A=',1P,G12.2,' e=',0P,f10.4)
      RETURN
      END

        SUBROUTINE  WCMOTION(hs)
        INCLUDE 'ARCCOM2e2.CH'
         COMMON/IncrementCommon/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
       common/vwcommon/Ww(nmx3),WTTLw,cmvw(3),spinw(3)
         common/omegacoefficients/OMEC(NMX,NMX)
         common/apuindex/ikir
        COMMON/DerOfTime/G
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
         real*8 FC(NMX3),XAUX(3),acc(nmx3)
         real*8 F(NMX3),!df(nmx3),dfGR(nmx3),
     &   GOM(nmx3)!,dcmv(3),Va(nmx3),afc(nmx3),dfE(3),dspin(3)
         save
         call EVALUATE X 
         RSUM=0.0
         OMEGA=0.0d0 
         U=0
         do i=1,3*N
         f(i)=0
         GOM(i)=0
         end do
         DO I=1,N-1
         L=3*(I-1)
         RIJL2=xc(L+1)**2+xc(L+2)**2+xc(L+3)**2+ee
         RIJL=SQRT(RIJL2)
C        Evaluate RSUM for decisionmaking.
         RSUM=RSUM+RIJL
         RINV(I)=1.d0/RIJL
         U=U+MC(I)*MC(I+1)*RINV(I)
         A=RINV(I)**3
         i0=3*i-3
         j=i+1
         j0=3*j-3
          omeker=omec(iname(i),iname(j))
         do K=1,3
          AF=A*XC(I0+K)
         f(I0+k)=f(i0+k)+MC(J)*AF
         f(j0+k)=f(j0+k)-MC(I)*AF
      if(cmethod(2).ne.0.0d0.and.omeker.ne.0.0)then
         GOM(I0+k)=GOM(I0+k)+AF*omeker
         GOM(J0+k)=GOM(J0+k)-AF*omeker
          end if
         end do
         if(cmethod(2).ne.0.0.and.omeker.ne.0.0)then
         OMEGA=OMEGA+omeker*RINV(I) 
         end if
         END DO

         LRI=N-1
C       Physical coordinates
        DO K=1,3
        XI(K)=0.0
        END DO
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        XI(L+3+K)=XI(L+K)+XC(L+K) 
        END DO
        END DO
C        Non-chained contribution
        DO I=1,N-2
        LI=3*(I-1)
        DO J=I+2,N  
        LJ=3*(J-1)
        RIJ2=0.0+ee
          IF(J.GT.I+2)THEN
           DO K=1,3
           XAUX(K)=XI(LJ+K)-XI(LI+K)
           RIJ2=RIJ2+XAUX(K)**2
           END DO
           ELSE
           DO K=1,3
           XAUX(K)=XC(LI+K)+XC(LI+K+3)
           RIJ2=RIJ2+XAUX(K)**2
           END DO
          END IF
        RIJ2INV=1/RIJ2
        LRI=LRI+1
        RINV(LRI)=SQRT(RIJ2INV)
          U=U+MC(I)*MC(J)*RINV(LRI)
          omeker=omec(iname(i),iname(j))
          if(omeker.ne.0.0.and.cmethod(2).ne.0.0)then
          OMEGA=OMEGA+omeker*RINV(LRI)
          end if
          DO K=1,3
          A=RINV(LRI)**3*XAUX(K)
          f(LI+K)=f(LI+K)+MC(J)*A 
          f(LJ+K)=f(LJ+K)-MC(I)*A
      if(cmethod(2).ne.0.0d0.and.omeker.ne.0.0)then
          GOM(LI+K)=GOM(LI+K)+A*omeker
          GOM(LJ+K)=GOM(LJ+K)-A*omeker
            end if
          END DO
         END DO ! J=I+2,N
        END DO  ! I=1,N-2
         dT=hs/(U*cmethod(1)+OMEGA*cmethod(2)+cmethod(3)) ! time interval
cccccccccccccc
c     Modified 6/23/2009 by D. Merritt
ccccccccccccc
c           call Coordinate Dependent Perturbations (acc) 
           call Coordinate Dependent Perturbati (acc) 
                 do i=1,n-1
                 do k=1,3
                 L=3*(i-1)
                 FC(L+k)=f(3*i+k)-f(3*i+k-3)
                 end do
                 end do
         IF(clight.gt.0.0)then       ! V-dependent ACC 
         call  V_jump(Ww,spinw,cmvw,WTTLw,WC,spin,FC,acc,dt/2
     &,gom,energyj,energrj,1) ! Auxiliary W (=Ww) etc
         call V_jump(WC,spin,cmv,WTTL,Ww,spinw,FC,acc,dt
     s,gom,energy,energr,2)   ! 'true' W  etc
         call  V_jump(Ww,spinw,cmvw,WTTLw,WC,spin,FC,acc,dt/2
     &,gom,energyj,energrj,3) ! Auxiliary W (=Ww) ets
         ELSE ! c>0
        call V_jACConly(WC,cmv,WTTL,FC,acc,dt,
     &  gom,energy,energrj)  ! here ACC depends ONLY on COORDINATES 
         END IF
        RETURN
        END


        subroutine V_jump(WCj,spinj,cmvj,wttlj,WCi,spini,FCj,acc,dt,
     &  gom,energyj,energrj,ind)
        include 'ARCCOM2e2.CH'
         COMMON/IncrementCommon/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        REAL*8 wcj(*),fcj(*),df(nmx3),dcmv(3),afc(nmx3),gom(*)
     &  ,dfe(nmx3),dfgr(nmx3),dspin(3),spinj(3),cmvj(3),wci(nmx3)
     &  ,spini(3),acc(*)
        save
        call EVALUATE V(V,WCi) 
c adding V-dependent perts.
              if(clight.gt.0.0)then
              call Velocity Dependent Perturbations
     &         (dT,V,spini,acc,dcmv,df,dfGR,dspin)
              else
              do i=1,3*n
              df(i)=acc(i)
              end do
              end if
                 do i=1,n-1
                 L=3*I-3
                 I1=3*INAME(I)-3
                 I2=3*INAME(I+1)-3
                 do k=1,3
                 afc(L+k)=df(I2+k)-df(I1+k) 
                 end do
                 end do
          IF(IND.EQ.2)THEN 
        dotE=0
        dotEGR=0
        do I=1,N
        I0=3*I-3
        do k=1,3
        dfE(k)=df(i0+k)-dfGR(i0+k)!
        end do
        dotE=dotE+! NB-Energy change (without Grav.Rad.)
     &    M(I)*(V(I0+1)*dfE(1)+V(I0+2)*dfE(2)+V(I0+3)*dfE(3)) ! %
        do k=1,3
        dfE(k)=dfGR(I0+k)
        end do
        dotEGR=dotEGR+ ! radiated energy 
     &    M(I)*(V(I0+1)*dfE(1)+V(I0+2)*dfE(2)+V(I0+3)*dfE(3))
        end do
                  ENERGYj=ENERGYj+dotE*dT
                  EnerGrj=EnerGRj+dotEGR*dT
                  if(ind.eq.2)then
                  ENERGYinc=ENERGYinc+dotE*dt
                  EnerGRinc=EnerGRinc+dotEGR*dT
                  end if !ind.eq.2
          END IF ! IND=2        
               if(cmethod(2).ne.0.0d0)then
        dotW=0
        do I=1,N
        k0=3*I-3
        i0=3*iname(i)-3
        dotW=dotW+
     &  (V(I0+1)*GOM(k0+1)+V(I0+2)*GOM(K0+2)+V(I0+3)*GOM(K0+3))
        end do

                  WTTLj=WTTLj+dotW*dT 
                if(ind.eq.2) WTTLinc=WTTLinc+dotW*dT
                end if ! cmethod(2).ne.0.0
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        if(ind.eq.2)WCinc(L+K)=WCinc(L+K)+(FCj(L+K)+afc(L+K))*dT 
        WCj(L+K)=WCj(L+K)+(FCj(L+K)+afc(L+K))*dT 
        END DO
        END DO

        do k=1,3
        spinj(k)=spinj(k)+dT*dspin(k)
        cmvj(k)=cmvj(k)+dT*dcmv(k)
        end do
        if(ind.eq.2)then
                do k=1,3
        spin inc(k)=spin inc(k)+dT*dspin(k)
        cmv inc(k)=cmv inc(k)+dT*dcmv(k)
        end do
        end if ! ind.eq.2
        RETURN
        END


        subroutine V_jACConly(WCj,CMVj,WTTLj,FC,acc,dt,
     &  gom,energyj,energrj)
        include 'ARCCOM2e2.CH'
         COMMON/IncrementCommon/WTTLinc,XCinc(NMX3),WCinc(NMX3),
     & CMXinc(3),CMVinc(3),ENERGYinc,Energrinc,CHTIMEinc,spin inc(3)
        REAL*8 wcj(*),fc(*),dcmv(3),afc(nmx3),gom(*)
     &  ,dfe(nmx3),cmvj(3),acc(*),WCi(NMX3)
        save
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        WCi(L+K)=WC(L+K)+FC(L+K)*dT/2  !( no inc here!)
        END DO
        END DO
        call EVALUATE V(V,WCi) 
        call reduce 2 cm(acc,m,n,dcmv) 
        do I=1,3*N
        V(I)=V(I)+acc(I)*dT/2 ! average Velocity
        end do
c adding V-dependent perts.
                 do i=1,n-1
                 L=3*I-3
                 I1=3*INAME(I)-3
                 I2=3*INAME(I+1)-3
                 do k=1,3
                 afc(L+k)=acc(I2+k)-acc(I1+k) ! CHAIN vector accelerations 
                 end do
                 end do
        dotE=0
        dotEGR=0
        do I=1,N
        I1=3*I-2
        do k=1,3
        dfE(k)=acc(i0+k)
        end do
        dotE=dotE+M(I)*cdot(V(I1),acc(i1))   
        end do
                  ENERGYj=ENERGYj+dotE*dT
                  EnerGrj=EnerGRj+dotEGR*dT

                 ENERGYinc=ENERGYinc+dotE*dT
                 EnerGRinc=EnerGRinc+dotEGR*dT

               if(cmethod(2).ne.0.0d0)then
        dotW=0
        do I=1,N
        k1=3*I-2
        i1=3*iname(i)-2
        dotW=dotW+cdot(V(I1),GOM(K1))
        end do
                  WTTLinc=WTTLinc+dotW*dT
                  WTTLj=WTTLj+dotW*dT 
                end if ! cmethod(2).ne.0.0
        DO I=1,N-1
        L=3*(I-1)
        DO K=1,3
        WCinc(L+K)=WCinc(L+K)+(FC(L+K)+afc(L+K))*dT
        WCj(L+K)=WCj(L+K)+(FC(L+K)+afc(L+K))*dT 
        END DO
        END DO

        do k=1,3
        cmv inc(k)=cmv inc(k)+dT*dcmv(k)
        cmvj(k)=cmvj(k)+dT*dcmv(k)
        end do
        RETURN
        END
