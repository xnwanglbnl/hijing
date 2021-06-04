        PROGRAM d_AU
        DIMENSION SCIP(300,300),RNIP(300,300),SJIP(300,300),JTP(3),
     &            IPCOL(90000),ITCOL(90000)
        COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
C
        COMMON/HIJCRDN/YP(3,300),YT(3,300)
        COMMON/HIJGLBR/NELT,NINT,NELP,NINP
        COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
        COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
        COMMON/HISTRNG/NFP(300,15),PP(300,25),NFT(300,15),PT(300,25)
        COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
        COMMON/HIJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
     &          K2SG(900,100),PXSG(900,100),PYSG(900,100),
     &          PZSG(900,100),PESG(900,100),PMSG(900,100)
        COMMON/HIJJET4/NDY,IADY(900,2),KFDY(900),PDY(900,5)
C*********************************************************************************
        CHARACTER FRAME*8,PROJ*8,TARG*8,FLNAME*10,HSTFILE*80
        double precision dNdy(-100:1000,3), dNdeta(-100:100,3)
     &    ,centrality(4), eta, y, absP,dy,  bb(4)
        data centrality/0.0, 0.145, 0.255, 1.00/
        data bb/0.0, 3.89, 5.21, 14.0/
C********* b value for 0.0%, 20%, 40%, 100%
C********* correspondingly in the definition of STAR
C********* http://iopscience.iop.org/1742-6596/5/1/003
        dy=0.1
        dNdy=0.0
        dNdeta=0.0

        WRITE(*,*) 'Enter the file name of the result.'
        READ(*,*) HSTFILE
        WRITE(*,*) 'Enter the number of events.'
        READ(*,*) N_EVENT
c        WRITE(*,*) 'Enter the energy.'
c        READ(*,*) EFRM
        EFRM=200.0

        open(10,file=HSTFILE)
C        EFRM=7000.0
        IHPR2(4)=0
        IHPR2(10)=0
        IHPR2(6)=1
        FRAME='CMS'
        PROJ='A'
        TARG='A'
        IAP=2
        IZP=1
        IAT=197
        IZT=79
c        IAT=203
c        IZT=82
        CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
        Do Icent=1,3
c         BMIN=SQRT(centrality(Icent))*(HIPR1(34)+HIPR1(35))*0.5
c         BMAX=SQRT(centrality(Icent+1))*(HIPR1(34)+HIPR1(35))*0.5

         Bmin=bb(Icent)
         Bmax=bb(Icent+1)
         percent=0.0
         partn=0.0
         avgb=0.0
         DO 1000 ,K=1,N_EVENT
                IF(REAL(K)/REAL(N_EVENT)
     &             .GT.PERCENT/100.0)THEN
                    WRITE(*,*) PERCENT,'%', Icent
                    PERCENT=PERCENT+1.0
                ENDIF
100     CALL HIJING (FRAME,BMIN,BMAX)
          nparti=np+nt
          partn=partn+real(nparti)/real(n_event)
          avgb=avgb+HINT1(19)/real(n_event)
          do 500, I=1,NATT
            IF(LUCHGE(KATT(I,1)) .EQ. 0) GOTO 500
            ABSP=SQRT(PATT(i,1)**2+PATT(i,2)**2+PATT(i,3)**2)
            ETA=0.5*LOG((ABSP+PATT(i,3))/(ABSP-PATT(i,3)))
            Y=0.5*LOG((PATT(i,4)+PATT(i,3))/(PATT(i,4)-PATT(i,3)))

            if(abs(Y)<=10.0) then
            IbinY=Int(Y/dy+1000.5)-1000
             dNdy(IbinY,Icent)=dNdy(IbinY,Icent)+1.0/dble(N_event)/dy
            endif
            if(abs(eta)<=10.0) then
            Ibineta=Int(eta/dy+1000.5)-1000
             dNdeta(Ibineta,Icent)=dNdeta(Ibineta,Icent)
     &                             +1.0/dble(N_event)/dy
            endif
500       continue
1000     CONTINUE
         write(*,*)partn, avgb
        Enddo

        do I=-100,100
          write(10,15) dble(I)*dy, dNdy(I,:), dNdeta(I,:)
        enddo
15    format(12(d14.7,1x))
        END


