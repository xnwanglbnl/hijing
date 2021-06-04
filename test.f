      program test
      CHARACTER FRAME*8,PROJ*8,TARG*8,dum*32
      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/
      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/
      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
      SAVE  /HIMAIN2/
      COMMON/HIMAIN3/VATT(130000,4)
      SAVE  /HIMAIN3/
      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
C
C       
C....initialize HIJING for Au+Au collisions at c.m. energy of 200 GeV:
C      WRITE(*,*)  'random number seed'
      READ(*,*) dum,NSEED
C      WRITE(*,*) 'frame(LAB,CMS), enegy-frame'
      READ(*,*)  dum,FRAME,EFRM
C      WRITE(*,*) 'Proj, Targ(A,P,PBAR)'
      READ(*,*)  dum,PROJ,TARG
C      WRITE(*,*)  'A,Z of proj; A,Z of targ'
      READ(*,*)  dum,IAP,IZP,IAT,IZT
C      WRITE(*,*) 'number of events'
      READ(*,*) dum,N_EVENT
      READ(*,*) dum,IHPR2(10)
C
C
C*** initialize HIJING
        CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
        write(*,*),'Sjet=', HINT1(11),'mb','Stot=',HINT1(13),'mb'
C
C....set BMIN=0 and BMAX=0.0 for central interaction
      BMIN=0.0
      BMAX=0.0
C....generating N_EVENT events of central AA interaction:
      DO 200 IE=1,N_EVENT
         CALL HIJING(FRAME,BMIN,BMAX)
         WRITE(*,*) IE,NATT,EATT
         do 100 i=1,NATT
            write(*,*) PATT(i,4),VATT(i,4),
     &           sqrt(VATT(i,1)**2+VATT(i,2)**2),VATT(i,3)
 100        continue
 200  continue
      STOP
      END 
