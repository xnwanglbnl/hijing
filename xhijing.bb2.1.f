c**********************************************************************
	PROGRAM XHIJINGbb20
cmg gyulassy +topor, includes new version 2 of JJbar via diquark pt kick
C gyulassy +topor, include SCF effects 
C 
c       
c       This program runs the HIJING/BB subroutine which simulates A+B 
c       nuclear collisions and places its output into various 
c       histograms (dN/dy,dN/dpt**2, etc). 
c       
c       Modified from a previous version by M. Gyulassy
c       12/98 - S.E.Vance.
c
c       To be put back into the scheme, the phi correlations of the last
c       event and the eta_cut. 
c       
c**********************************************************************       

	double precision dd1,dd2
	integer ippad,itpad
	character frame*8,targ*8,proj*8,hstfile*80,numfile*80
 	character descr*48, chtit*80
	dimension pthbk(219)
	dimension id(100)
	dimension xhb(400),yhb(400),exhb(400),eyhb(400)
c       
c       HIJING common blocks
c       
	common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
	save  /hiparnt/
     
       
        common/hijdat/hidat0(10,10),hidat(10)
        save  /hijdat/
        common/hijcrdn/yp(3,300),yt(3,300)
        save  /hijcrdn/
        COMMON/HIJGLBR/NELT,NINT,NELP,NINP
        SAVE  /HIJGLBR/      
 
	common/himain1/natt,eatt,jatt,nt,np,n0,n01,n10,n11
	save  /himain1/
	common/himain2/katt(130000,4),patt(130000,4)
	save  /himain2/
 	common/hijjet1/npj(300),kfpj(300,500),pjpx(300,500),pjpy(300,500)
     1    ,pjpz(300,500),pjpe(300,500),pjpm(300,500)
     2    ,ntj(300),kftj(300,500),pjtx(300,500),pjty(300,500)
     3    ,pjtz(300,500),pjte(300,500),pjtm(300,500)
	save  /hijjet1/
 	common/hijjet2/nsg,njsg(900),iasg(900,3),k1sg(900,100),k2sg(900,100)
     1    ,pxsg(900,100),pysg(900,100),pzsg(900,100),pesg(900,100)
     2    ,pmsg(900,100)
	save  /hijjet2/
c add topor 
        COMMON/HIJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
	SAVE  /HIJJET4/
	COMMON/RANSEED/NSEED
	SAVE  /RANSEED/


c       


	common/histrng/nfp(300,15),pp(300,15),nft(300,15),pt(300,15)
	save  /histrng/
c       
	common/lujets/n_lu,k_lu(9000,5),p_lu(9000,5),v_lu(9000,5)   
	save  /lujets/
	common/ludat1/mstu(200),paru(200),mstj(200),parj(200)
	save  /ludat1/
c       
        COMMON/HIJB/nb,ib(15),rb(15),sb(5),kb(10,900),kfb(10,15),
     &  pb(10,15),ncb(15),nbjet(10),kbjet(10,100),pbjet(10,100,5)   
        SAVE  /HIJB/

C old-topor     &       pb(10,15),ncb(15)
c       
cmg
      common/JJbarFrag/ Ifudge, Fudgefac, Probab
cmg

        COMMON/NCOMPARE/nexp
        SAVE /NCOMPARE/
c
c	PAW common blocks
c       
 	common/pawc/hmemor(300000)    
c       
c       Read input from data file
c       
c       Energy per nucleon pair (GeV) and frame (CMS or LAB)    
	read (5,*) descr,efrm,frame  
	write(6,*) descr,efrm,' ',frame
c       Type of projectile particle (P,PBARor A), its charge number and 
c       its atomic number. 
	read (5,*) descr,proj,izp,iap
	write(6,*) descr,proj,izp,iap
c       Type of target particle (P,PBARor A), its charge number and 
c       its atomic number.     
	read (5,*) descr,targ,izt,iat
	write(6,*) descr,targ,izt,iat
c       Impact parameters for the collision
	read (5,*) descr,bmin,bmax
	write(6,*) descr,bmin,bmax
c       Number of Events
	read (5,*) descr,n_event
	write(6,*) descr,n_event
c       Seed for the Random Number Generator 
	read (5,*) descr,nseed
	write(6,*) descr,nseed
c       Names of the histogram and numbers file. 
	read (5,*) descr,hstfile,numfile
	write(6,*) descr,hstfile,numfile
	write(6,*) ' '
cmg
c       Flag and cross sec for B junction
	read (5,*) descr,IflagJ, sigJ 
	write(6,*) descr,IflagJ, sigJ
	write(6,*) ' '
c       Flag and cross sec for BBar junction simulation
	read (5,*) descr, Ifudge, sigJJ 
	write(6,*) descr, Ifudge, sigJJ
	write(6,*) ' '
c       pt -> pt * Fudgefac if BBar occurs with Probab determined in
c       hijing subroutine (see cmg comments there)
	read (5,*) descr, Fudgefac
	write(6,*) descr, Fudgefac
	write(6,*) ' '
cmg

cmg

	open(23,file=numfile,status='unknown')

c       set options here
c       These are the parameters set by M. Gyulassy and V. Topor Pop
	mg = 1
	if (mg.eq.1) then
c       IHPR2(8):(D=10)The maximum number of jet production per 
c               nucleon-nucleon interaction. When IHPR2(8)=0, jet
c               production will be turned off. When IHPR2(8)<0, the
c               number of jets production will be fixed for each
c               NN collision at a value of abs(IHPR2(8)).
	   ihpr2(8) = 10

c       IHPR2(3):(D=0)switch for one hard scattering with specified PT
c		PT=HIPR1(10) per event.
c		=1:ordinary hard processes;
c		=2:only direct photon production.
c               =3:only heavy quark production 
C ************ihpr2(3) = 0
C test triger on charm 
C               ihpr2(3)=3              
                       ihpr2(3) = 0
                           

c       HIPR1(10):(D=-2.25 GeV) To specify the value of PT for one hard
c		scattering generated per event(we call them sample 
c		jets), which could be used to study jet in AA 
c		collision. If HIPR1(10) is negative, its absolute value
c		gives the low limit of the PT of the sample jets.
c               p_t triger jet
C *****	   hipr1(10) = 30.0
C *****topor test charm 
C ******for inclusive charm production we have to reset hipr1(10) to 0.0
C ***** One can also trig larger pt QQbar production by giving HIPR1(10)
C ***** a non vanishing value  
C *****************            hipr1(10) = 0.0

                       hipr1(10) = 30.0

c	IHPR2(6):(D=1)switch for the nuclear effect on the parton 
c		distribution function; 
c            shadow flag  yes=1 no =0
 	   ihpr2(6) = 1

c       IHPR2(4):(D=1)switch for jet quenching in the the exicited 
c               nuclear matter,i.e. jets only interact with those 
c               nuclear matter that has suffered collisions.
c         quench flag yes=1 no=0
	   ihpr2(4) = 1

c	IHPR2(5):(D=1)switch for the pt kick due to soft initial state 
c		interaction both for valence quarks or diquarks and 
c		gluons. 
c          p_t kick flag
	   ihpr2(5) = 1

c	IHPR2(2):(D=3)switch for initial and final state radiation in 
c		the hard scattering.
c		=0: both initial and final radiation off;
c		=1: initial radiation on and final radiation off;
c		=2: initial radiation off and final radiation on;
c		=3: both initial and final radiation on.
c        radiation flag
	   ihpr2(2) = 3

c	HIPR1(13):(D=1.0 fm)The mean free path of a jet to interact 
c		when it goes through an excited nuclear matter.
c           quench lambda 
	   hipr1(13) = 1.

c	HIPR1(14):(D=2.0 GeV/fm) The string tension of a gluon jet. 
c		It is the energy loss dE/dx=HIPR1(14) of the gluon 
cc		jet when the string is streched. For quark jet 
c		dE/dx=0.5*HIPR1(14). 
c                quench dE/dx

C**************       	      hipr1(14) = 2.0

                         hipr1(14) = 2.0

c	HIPR1(17):(D=2.0) K factor for the differential jet cross sections
C topor increase to 3.5 for RHIC(Eskola et al)
C *************            HIPR1(17) = 3.5
C topor increase to 6.0 for SPS
                HIPR1(17) = 2.0 

c	HIPR1(23):(D=2.0) 
C       A separate K factor for the differential jet cross sections
C       for heavy quark and direct photon         
C *************            HIPR1(23) = 4.0 test

                HIPR1(23) = 2.0

c	HIPR1(7):(D=1.5) 
C       Quark mass for the differential jet cross sections
C       for heavy quark    
C ************** HIPR1(7) = test mc=1.25 /sqrt(2)      
C *************  HIPR1(7) = test mc=1.25 /sqrt(3) 
C ************** HIPR1(7) = test mc=1.25 /sqrt(5)
C ***************HIPR1(7) = test mc=1.25 /sqrt(1.5)

C *************  HIPR1(7) = test mc=1.25 /sqrt(3)

                 HIPR1(7) = 1.25


c	HIPR1(12):(D=1.0 fm)The transverse distance between a 
c		traversing jet and an excited nucleon(or string system)
c 		below which they will interact and the jet will lose
c		energy and momentum to the nucleon.
c              quench sigma radius
	   hipr1(12) = 1.0

C***********************************************************
C topor Modified Cronin parameters for SPS energies
C topor hipr1(19) (D= 0.1);hipr1(20) (D= 1.4);hipr1(21) (D= 1.6)
C           hipr1(19) = 0.1
C           hipr1(20) = 1.4
C           hipr1(21) = 1.6
C topor *****************************************************


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c               ' no single diffrac flag     '     0
c            ihpr2(13)(D=1) Option to turn on single
c           diffractive reactions'

C              IHPR2(13)=0
C test nsdoff
C               IHPR2(13) = 1
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	HIPR1(29):(D=0.4 fm) the minimum distance between two nucleons
c 		inside a nucleus when the coordinates of all the
c		nucleons in a nucleus are initialized.
c                 hard core radius
	   hipr1(29) = 0.4

c	IHPR2(20):(D=1)Option to turn off(IHPR2(20)=0) the 
c		fragmentation.
c                fragmentation on off 
	   ihpr2(20) = 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c             ' pi0 decay, ioflag ' 0 0    
c                 ihpr2(12)=0

C            IHPR2(12)=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c	HIPR1(11):(D=2.0 GeV) Minimum pt of a jet which will interact 
c		with excited nuclear matter, or when the pt of a jet 
c		is smaller than HIPR1(11) it will stop interacting 
c		further more.
c        ' min glue energy for quench '    2.0
	   hipr1(11) = 2.0

c	HIPR1(8):(D=2.0 GeV) Minimum pt of hard or semihard scatterings.
c       ' mini jet p0                '    2.0 130/ and 2.18 200
c topor HIPR1(8)= 1.7 at sqrt{s_NN}=17 GeV
c topor HIPR1(8)= 2.0 at sqrt{s_NN}=130 GeV
c topor HIPR1(8)= 2.18 at sqrt{s_NN}=200 GeV
c topor HIPR1(8)= 3.5 at  sqrt{s_NN}=5500 GeV  
c topor for other energies use Xin Nian Wang parametrization 
C topor  Phys. Lett. B 527, 85(2002). 

           hipr1(8) = 2.20
	

c	HIPR1(30):(D=2*HIPR1(31)=57.0 mb) the inclusive cross 
c		section sigma_s for soft interaction. The default
c		value sigma_s=sigma_0 is used to insure the geometrical
c		scaling of pp interaction cross section at low energies.
	   HIPR1(30) = 57.0

c	HIPR1(31):(D=28.5 mb) the cross section sigam_0 which 
c		characterizes the geometrical size of a nucleon
c		(HIPR1(31)=pi*b0**2). The default value is only
c		for high-energy limite(sqrt{s}>200 GeV).
	   hipr1(31) = 28.5

c       HIPR1(32):(D=3.90) The parameter in the scaled eikonal function 
c		of nucleon Xi(R)=A**2*(A*R)**3*K3(A*R),A=HIPR1(32) 
c		which is a Fourier transform of a dipole formtor 
	   hipr1(32) = 3.90
	end if 
	nbab = 1
c       set parameters for the baryon junctions. 
	if (nbab.eq.1) then
c
c       string fragmentation parameters, see manual for Jetset v.7.2. 
c       parj(1):(D=0.10),parj(2):(D=0.30),parj(3):(D=0.40),parj(4):(D=0.05),
c       parj(11):(D=0.5)


c

C topor parameters p1=0.07 PYTHIA 
C          parj(1)=0.07
C          parj(2)=0.30
C          parj(3)=0.40
C          parj(4)=0.05
C          parj(11)=0.5
C          parj(21)=parj(21)


C initial parameters , modified by irope1=1

         parj(1)=0.02
         parj(2)=0.30
         parj(3)=0.40
         parj(4)=0.05
         parj(11)=0.5


c       irope1, rkappa1:
c	The variable rope indicates a larger string tension (SCF).  
C        The variable
c	rkappa1 indicates the value of the new tension divided by the value
c	of the old tension.
c
	       
              irope1 = 1

C topor  rkappa1  and parameter f3 from *dat file
C  topor have different values for pp, dau, AuAu
C in *.dat file f3=1  without JJbar loops
C in *.dat file f3=3 with JJbar loops
C*************************************************
C topor             rkappa1=1.0 (without SCF effects)
C topor             rkappa1=1.5 ; f3=1 for pp@200 GEV
C******************************************************
C topor             rkappa1=2.0 ;f3=3 

C topor       for dAu , AuAu (non-strange+strange particles)

c topor             rkappa1=3.0; f3=3 for Au+Au multi-strange particles
C**********************************************************************

               rkappa1= 3.0

	   if (irope1.eq.1) then
	      parj(1) = parj(1)**(1.0/rkappa1)
            
              parj(2) = parj(2)**(1.0/rkappa1)

              parj(3)=parj(3)**(1.0/rkappa1)

	      parj(4) = 1.0/3.0*((3.0*parj(4))**(1.0/rkappa1))

              parj(21) = parj(21)*sqrt(rkappa1)
	   end if
 
           write(6,*) 'parj(1,2,3,4,21)',parj(1),parj(2),parj(3),
     &     parj(4),parj(21)

           write(6,*) 'hipr1(17)=', hipr1(17)
           write(6,*) 'hipr1(23)=', hipr1(23)
           write(6,*) 'hipr1(7)=', hipr1(7)
c
c       ib(1):  0,1 for switching off/on the baryon stopping interactions
c
	   ib(1) = IflagJ
c       
c       ib(2):  0,1 for switching off/on the 1 baryon stopping interaction
c       rb(2):  cross section (mb) of the 1 baryon stopping component
c       
	   ib(2) = IflagJ
	   rb(2) = sigJ
c       
c       rb(4):   minimum mass of the Y string configuration 
c       
	   rb(4) = 5.0
c       
c       ib(4):   pT kick for the junction baryon  
c       
	   ib(4) = 1
c
c       ib(5):  switch for a full print out of the baryon junction effect
c
	   ib(5) = 0
c
c	ib(6): 0,1 for switching off/on the baryon - anti-baryon production
c       rb(6): cross section (mb) of the baryon - anti-baryon production
c       
	   ib(6) = Ifudge
	   rb(6) = sigJJ
c       
c       rb(7): minimum mass of the baryon - anti-baryon configuration. 
c       
	   rb(7) = 6.0
c
c       sb(4): percentage of the junctions which remain junctions in a 
c       subsequent baryon junction exchange interaction.
c
	   sb(4) = 0.9

c
c	ib(11):  baryon stopping with diquark breaking
c	rb(11):  cross section (mb) for the diquark breaking
c	
	ib(11) = 0
	rb(11) = 0.0
c	
	end if
        
c

c	nbje = 0
c       set parameters for the baryon junctions. 
c	if (nbje.eq.1) then
c
c       string fragmentation parameters, see manual for Jetset v.7.2. 
c       parj(1):(D=0.10),parj(2):(D=0.30),parj(3):(D=0.40),parj(4):(D=0.05),
c       parj(11):(D=0.5)
c
c           parj(1)=0.07
c           parj(2)=0.30
c           parj(3)=0.40
c           parj(4)=0.05
c           parj(11)=0.5
c
c
c       ib(1):  0,1 for switching off/on the baryon stopping interactions
c
c	   ib(1) = 1
c       
c       ib(2):  0,1 for switching off/on the 1 baryon stopping interaction
c       rb(2):  cross section (mb) of the 1 baryon stopping component
c       
c	   ib(2) = 1
c	   rb(2) = 18.0
c       
c       rb(4):   minimum mass of the Y string configuration 
c       
c	   rb(4) = 5.0
c       
c       ib(4):   pT kick for the junction baryon  
c       
c	   ib(4) = 1
c
c       ib(5):  switch for a full print out of the baryon junction effect
c
c	   ib(5) = 0
c
c	ib(6): 0,1 for switching off/on the baryon - anti-baryon production
c       rb(6): cross section (mb) of the baryon - anti-baryon production
c       
c	   ib(6) = 0
c	   rb(6) = 0.0
c       
c       rb(7): minimum mass of the baryon - anti-baryon production. 
c       
c	   rb(7) = 0.0
c
c	ib(11):  baryon stopping with diquark breaking
c	rb(11):  cross section (mb) for the diquark breaking 
c	
c	ib(11) = 0
c	rb(11) = 0.0
c	
c	end if


 



 
c       The cut on eta for the pt distributions. 
	eta_cut = 0.5

c
c
c       scenario with which to compare (relates to the decay of the
c       baryons).
	nexp = 2

c
	call hlimit(300000)
	
c       initialize HIJING
	
	call hijset(efrm,frame,proj,targ,iap,izp,iat,izt)
	
	write(6,*),'HIJING was initialized with following cros secs (mb):'
	write(6,*),'HINT1(10)=',HINT1(10),' cross sec for jet production'
        write(6,*),'HINT1(11)=',HINT1(11),' inclusive jet cross section'
	write(6,*),'HINT1(12)=',HINT1(12),' inelastic cross section'
	write(6,*),'HINT1(13)=',HINT1(13),' total cross section'
	write(6,*),' '

c	
c       create histograms for the event geometry.  These numbers of the 
c       histograms range from 20-99.
c topor add nb histo
c           call hbook1(31,'NCB(i)',20.0,0.0,20,0.0)
C         call hbook1(32,'NCB(2)',1000.0,0.0,10,0.0)
C         call hbook1(33,'NCB(3)',1000.0,0.0,10,0.0)
C         call hbook1(34,'NCB(4)',1000.0,0.0,10,0.0)
C         call hbook1(31,'NCB(i)',20.0,0.0,20,0.0)
C         call hbook1(25,'BJ(4) or JBARJ(6)',20.0,0.0,10,0.0)  
	ippad = ( 0.1*ihnt2(1))
	itpad = ( 0.1*ihnt2(3)) 
	if(ihnt2(1).gt.1) then  
       
	   call hbook1(50,'projectile participants',ihnt2(1)+ippad,0.0,
     1	  real(ihnt2(1)+ippad),0.0)
	   call hidopt(50,'logy')
	   ipt=ihnt2(1)+ihnt2(3)+ippad+itpad
	   call hbook1(52,'total participants',ipt,0.0,real(ipt),0.0)
	   call hidopt(52,'logy')
	   call hbook1(53,'proj eveto/ebeam',ihnt2(1)+ippad,0.0,
     1	   real(ihnt2(1)+ippad),0.0)
	   call hidopt(53,'logy')
	   xxmax = 2.4*ihnt2(1)**0.3333
	   delrr = 0.2
	   nnmax = xxmax/delrr
	   xxmax = nnmax*delrr
	   call hbook1(55,'rho(r), proj',nnmax,0.0,xxmax,0.0)
	   call hbook1(58,'dn(r), proj',nnmax,0.0,xxmax,0.0)
	endif
	if(ihnt2(3).gt.1) then
	   call hbook1(51,'participant targ',ihnt2(3)+itpad,0.0,
     1	  real(ihnt2(3)+itpad),0.0)
	   call hidopt(51,'logy')
	   call hbook1(54,'targ eveto/etarg',ihnt2(3)+itpad,0.0,
     1	   real(ihnt2(3)+itpad),0.0)
	   call hidopt(54,'logy')
	   xxmax = 2.4*ihnt2(3)**0.3333
	   delrr = 0.2
	   nnmax = xxmax/delrr
	   xxmax = nnmax*delrr
	   call hbook1(56,'rho(r), targ',nnmax,0.0,xxmax,0.0)
	   call hbook1(59,'dn(r),  targ',nnmax,0.0,xxmax,0.0)
	   call hbook1(57,'number of binary',100,0.0,
     1	  1.5*ihnt2(3)**1.3333,0.0)
	endif

c
c       Create histograms for pt distributions.
c       These histograms range from 100-199. 
c
	do 10 i=1,21
 	   pthbk(i)=(i-1)*0.05
 10	continue
	do 20 i=1,6
	   pthbk(21+i)=1.0+0.1*i
 20	continue
	do 30 i=1,192
	   pthbk(27+i)=1.6+i*0.2
 30	continue
	n_ptmx=218

C \beg topor
c       pt distributions for the photons.
 	call hbookb(115,'dN/dpt**2, photon(res decay not pi0)',
     1	      n_ptmx,pthbk,0.)
	call hidopt(115,'logy')
	call hbookb(116,'dN/dpt**2, direct photon',n_ptmx,pthbk,0.)
	call hidopt(116,'logy')
	call hbookb(117,'dN/dpt**2, pi0',n_ptmx,pthbk,0.)
	call hidopt(117,'logy')
	call hbookb(118,'dN/dpt**2, greekro0',n_ptmx,pthbk,0.)
	call hidopt(118,'logy')
	call hbookb(1118,'dN/dpt**2, greekomega',n_ptmx,pthbk,0.)
	call hidopt(1118,'logy')         
	call hbookb(119,'dN/dpt**2, greeketa',n_ptmx,pthbk,0.)
	call hidopt(119,'logy')
	call hbookb(1119,'dN/dpt**2, greeketaprim',n_ptmx,pthbk,0.)
	call hidopt(1119,'logy')
	call hbookb(1133,'dN/dpt**2, greekphi',n_ptmx,pthbk,0.)
	call hidopt(1133,'logy')

c       pt distributions for the mesons.
	call hbookb(101,'dN/dpt**2, charged',n_ptmx,pthbk,0.)
	call hidopt(101,'logy')
	call hbookb(120,'dN/dpt**2, pi+',n_ptmx,pthbk,0.)
	call hidopt(120,'logy')
	call hbookb(121,'dN/dpt**2, pi-',n_ptmx,pthbk,0.)
	call hidopt(121,'logy')
	call hbookb(122,'dN/dpt**2, K+',n_ptmx,pthbk,0.)
	call hidopt(122,'logy')
	call hbookb(123,'dN/dpt**2, K-',n_ptmx,pthbk,0.)
	call hidopt(123,'logy')

        call hbookb(126,'dN/dpt**2, Ks^0',n_ptmx,pthbk,0.)
	call hidopt(126,'logy')
        call hbookb(127,'dN/dpt**2, Ks^0~',n_ptmx,pthbk,0.)
	call hidopt(127,'logy')
        
        call hbookb(140,'dN/dpt**2, D+(+-)',n_ptmx,pthbk,0.)
	call hidopt(140,'logy')
        call hbookb(1140,'dN/dpt**2, resD+*(+-)',n_ptmx,pthbk,0.)
	call hidopt(1140,'logy')
        


        call hbookb(141,'dN/dpt**2, D0(+-)',n_ptmx,pthbk,0.)
	call hidopt(141,'logy')
        call hbookb(1141,'dN/dpt**2, resD0*(+-)',n_ptmx,pthbk,0.)
	call hidopt(1141,'logy')

        call hbookb(142,'dN/dpt**2, D_s+(+-)',n_ptmx,pthbk,0.)
	call hidopt(142,'logy')        
        call hbookb(1142,'dN/dpt**2, resD_s+*(+-)',n_ptmx,pthbk,0.)
	call hidopt(1142,'logy') 

        call hbookb(143,'dN/dpt**2, B+(+-)',n_ptmx,pthbk,0.)
	call hidopt(143,'logy')
        call hbookb(1143,'dN/dpt**2, resB+*(+-)',n_ptmx,pthbk,0.)
	call hidopt(1143,'logy')
         

        call hbookb(144,'dN/dpt**2, B0(+-)',n_ptmx,pthbk,0.)
	call hidopt(144,'logy')
        call hbookb(1144,'dN/dpt**2, resB0*(+-)',n_ptmx,pthbk,0.)
	call hidopt(1144,'logy') 

        call hbookb(145,'dN/dpt**2, B_s0(+-)',n_ptmx,pthbk,0.)
	call hidopt(145,'logy')
        call hbookb(1145,'dN/dpt**2, resB_s0*(+-)',n_ptmx,pthbk,0.)
	call hidopt(1145,'logy')
     
        call hbookb(146,'dN/dpt**2, L_c+',n_ptmx,pthbk,0.)
	call hidopt(146,'logy')
         
        call hbookb(147,'dN/dpt**2, B_c+',n_ptmx,pthbk,0.)
	call hidopt(147,'logy')
        
        call hbookb(1147,'dN/dpt**2, resB_c+*',n_ptmx,pthbk,0.)
	call hidopt(1147,'logy')

        call hbookb(148,'dN/dpt**2, J/PSI',n_ptmx,pthbk,0.)
	call hidopt(148,'logy') 

        call hbookb(149,'dN/dpt**2, Upsilon',n_ptmx,pthbk,0.)
	call hidopt(149,'logy')

c       pt distributions for the baryons. 
	call hbookb(150,'dN/dpt**2, p',n_ptmx,pthbk,0.)
	call hidopt(150,'logy')
	call hbookb(151,'dN/dpt**2, p~',n_ptmx,pthbk,0.)
	call hidopt(151,'logy')


        call hbookb(154,'dN/dpt**2, Lambda',n_ptmx,pthbk,0.)
	call hidopt(154,'logy')
        call hbookb(155,'dN/dpt**2, Lambda~',n_ptmx,pthbk,0.)
	call hidopt(155,'logy')

        call hbookb(156,'dN/dpt**2, Sigma-',n_ptmx,pthbk,0.)
	call hidopt(156,'logy')
        call hbookb(157,'dN/dpt**2, Sigma~+',n_ptmx,pthbk,0.)
	call hidopt(157,'logy')
        call hbookb(158,'dN/dpt**2, Sigma0',n_ptmx,pthbk,0.)
	call hidopt(158,'logy')
        call hbookb(159,'dN/dpt**2, Sigma~0',n_ptmx,pthbk,0.)
	call hidopt(159,'logy')
        call hbookb(160,'dN/dpt**2, Sigma+',n_ptmx,pthbk,0.)
	call hidopt(160,'logy')
        call hbookb(161,'dN/dpt**2, Sigma~-',n_ptmx,pthbk,0.)
	call hidopt(161,'logy')       

        call hbookb(162,'dN/dpt**2, Xi-',n_ptmx,pthbk,0.)
	call hidopt(162,'logy')
        call hbookb(163,'dN/dpt**2, Xi~+',n_ptmx,pthbk,0.)
        call hidopt(163,'logy')
        call hbookb(164,'dN/dpt**2, Xi0',n_ptmx,pthbk,0.)
        call hidopt(164,'logy')
        call hbookb(165,'dN/dpt**2, Xi~0',n_ptmx,pthbk,0.)
        call hidopt(165,'logy')      

        call hbookb(166,'dN/dpt**2, Omega-',n_ptmx,pthbk,0.)
        call hidopt(166,'logy')
        call hbookb(167,'dN/dpt, Omega~+',n_ptmx,pthbk,0.)
        call hidopt(167,'logy')   
    
        call hbookb(103,'dN/dpt**2, ch+',n_ptmx,pthbk,0.)
        call hidopt(103,'logy')
        call hbookb(104,'dN/dpt**2, ch-',n_ptmx,pthbk,0.)
        call hidopt(104,'logy')
c /end topor   




c       Create histograms for the rapidity distributions. 
c       These histograms range from 200-299.
c
	dd1=hint1(1)
	dd2=dsqrt(1.0-0.078/dd1**2)
	ymax=int(0.5*dlog((1.0+dd2)/(1.0-dd2)))
	if(frame.eq.'LAB') then
	   ymin=-2.0
	   ymax=2.0*ymax
	else
	   ymin=-ymax	   
        endif
	nymx=(ymax-ymin)/0.2
	w_y=(ymax-ymin)/real(nymx)

C \beg topor
c       Create the rapidity distributions of the q's, g's and photons. 
	call hbook1(201,'dN/dy, charged',nymx,ymin,ymax,0.0)
        call hbook1(202,'dN/dy, neutral',nymx,ymin,ymax,0.0)
        call hbook1(203,'dN/dy, ch+',nymx,ymin,ymax,0.0)
        call hbook1(204,'dN/dy, ch-',nymx,ymin,ymax,0.0)
	CALL HBOOK1(205,'dN/dy,  gluons',NYMX,YMIN,YMAX,0.0)
	CALL HBOOK1(206,'dN/dy,  q+qbar',NYMX,YMIN,YMAX,0.0)
	CALL HBOOK1(207,'dN/dy,  qq+qqbar',NYMX,YMIN,YMAX,0.0)
	CALL HBOOK1(210,'dN/dy,  c+cbar',NYMX,YMIN,YMAX,0.0)

	call hbook1(215,'dN/dy, photons(res. decay not pi0'
     1   ,nymx,ymin,ymax,0.0)
        call hbook1(216,'dN/dy, direct photons',nymx,ymin,ymax,0.0)
	call hbook1(217,'dN/dy, pi0',nymx,ymin,ymax,0.0)	
        call hbook1(218,'dN/dy, greekro0',nymx,ymin,ymax,0.0)
        call hbook1(1218,'dN/dy, greekomega',nymx,ymin,ymax,0.0)
        call hbook1(219,'dN/dy, greeketa',nymx,ymin,ymax,0.0)
        call hbook1(1219,'dN/dy, greeketa-prim',nymx,ymin,ymax,0.0)
         call hbook1(2233,'dN/dy, greekphi',nymx,ymin,ymax,0.0)        
        	

c       Create the rapidity distributions of the mesons.  

	call hbook1(220,'dN/dy, pi+',nymx,ymin,ymax,0.0)
	call hbook1(221,'dN/dy, pi-',nymx,ymin,ymax,0.0)
	call hbook1(222,'dN/dy, K+',nymx,ymin,ymax,0.0)
	call hbook1(223,'dN/dy, K-',nymx,ymin,ymax,0.0)
	call hbook1(224,'dN/dy, K0',nymx,ymin,ymax,0.0)
	call hbook1(225,'dN/dy, K~0',nymx,ymin,ymax,0.0)
	call hbook1(226,'dN/dy, K_S0',nymx,ymin,ymax,0.0)
	call hbook1(227,'dN/dy, K_s0~',nymx,ymin,ymax,0.0)
	call hbook1(228,'dN/dy, K*0',nymx,ymin,ymax,0.0)

	call hbook1(240,'dN/dy, D(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(2240,'dN/dy, resD(+,-)*',nymx,ymin,ymax,0.0)

        call hbook1(241,'dN/dy, D0(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(2241,'dN/dy, resD0*(+,-)',nymx,ymin,ymax,0.0)
        
         call hbook1(242,'dN/dy, D_s+(+,-)',nymx,ymin,ymax,0.0)
         call hbook1(2242,'dN/dy, resD_s+*(+,-)',nymx,ymin,ymax,0.0)


        call hbook1(243,'dN/dy, B(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(2243,'dN/dy, resB*(+,-)',nymx,ymin,ymax,0.0)
 
        call hbook1(244,'dN/dy, B0(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(2244,'dN/dy, resB0*(+,-)',nymx,ymin,ymax,0.0) 
        
        call hbook1(245,'dN/dy, B_s0(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(2245,'dN/dy, resB_s0*(+,-)',nymx,ymin,ymax,0.0)

         call hbook1(246,'dN/dy, L_c+(+,-)',nymx,ymin,ymax,0.0)
      
         call hbook1(247,'dN/dy, B_c+(+,-)',nymx,ymin,ymax,0.0)
         call hbook1(2247,'dN/dy, resB_c+*(+,-)',nymx,ymin,ymax,0.0) 
                  
         call hbook1(248,'dN/dy, J/PSI',nymx,ymin,ymax,0.0)
         
         call hbook1(249,'dN/dy, Upsilon',nymx,ymin,ymax,0.0)

c       Create the rapidity distributions of the baryons.  
	call hbook1(250,'dN/dy, p',nymx,ymin,ymax,0.0)
	call hbook1(251,'dN/dy, p~',nymx,ymin,ymax,0.0)
	call hbook1(252,'dN/dy, n',nymx,ymin,ymax,0.0)
	call hbook1(253,'dN/dy, n~',nymx,ymin,ymax,0.0)
	call hbook1(254,'dN/dy, Lambda',nymx,ymin,ymax,0.0)
	call hbook1(255,'dN/dy, Lambda~',nymx,ymin,ymax,0.0)	
	call hbook1(256,'dN/dy, Sigma-',nymx,ymin,ymax,0.0)
	call hbook1(257,'dN/dy, Sigma~+',nymx,ymin,ymax,0.0)
	call hbook1(258,'dN/dy, Sigma0',nymx,ymin,ymax,0.0)
	call hbook1(259,'dN/dy, Sigma~0',nymx,ymin,ymax,0.0)
	call hbook1(260,'dN/dy, Sigma+',nymx,ymin,ymax,0.0)
	call hbook1(261,'dN/dy, Sigma~-',nymx,ymin,ymax,0.0)
	call hbook1(262,'dN/dy, Xi-',nymx,ymin,ymax,0.0)
	call hbook1(263,'dN/dy, Xi~+',nymx,ymin,ymax,0.0)
	call hbook1(264,'dN/dy, Xi0 ',nymx,ymin,ymax,0.0)
	call hbook1(265,'dN/dy, Xi~0 ',nymx,ymin,ymax,0.0)
	call hbook1(266,'dN/dy, Omega-',nymx,ymin,ymax,0.0)
	call hbook1(267,'dN/dy, Omega~+',nymx,ymin,ymax,0.0)




c
c       Create histograms for the pseudo-rapidity distributions. 
c       These histograms range from 300-399.
c
	call hbook1(301,'dN/deta,  charged',nymx,ymin,ymax,0.0)
        call hbook1(302,'dN/deta,  neutral',nymx,ymin,ymax,0.0)
        call hbook1(303,'dN/deta,  ch+',nymx,ymin,ymax,0.0)
        call hbook1(304,'dN/deta,  ch-',nymx,ymin,ymax,0.0)
	call hbook1(315,'dN/deta,  photons',nymx,ymin,ymax,0.0)
        call hbook1(317,'dN/deta,  pi0',nymx,ymin,ymax,0.0)
c       
c       Create the pseudo-rapidity distributions of the mesons. 
	call hbook1(320,'dN/deta, pi+',nymx,ymin,ymax,0.0)
	call hbook1(321,'dN/deta, pi-',nymx,ymin,ymax,0.0)
	call hbook1(322,'dN/deta, K+',nymx,ymin,ymax,0.0)
	call hbook1(323,'dN/deta, K-',nymx,ymin,ymax,0.0)
	call hbook1(324,'dN/deta, K0',nymx,ymin,ymax,0.0)
	call hbook1(325,'dN/deta, K~0',nymx,ymin,ymax,0.0)
	call hbook1(326,'dN/deta, K_S0',nymx,ymin,ymax,0.0)
	call hbook1(327,'dN/deta, K_s0~',nymx,ymin,ymax,0.0)
	call hbook1(328,'dN/deta, K*0',nymx,ymin,ymax,0.0)
	
        call hbook1(340,'dN/deta, D+(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(3340,'dN/deta, resD+*(+,-)',nymx,ymin,ymax,0.0)
 
        call hbook1(341,'dN/deta, D0(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(3341,'dN/deta, resD0*(+,-)',nymx,ymin,ymax,0.0)

        call hbook1(342,'dN/deta, D_s+(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(3342,'dN/deta, resD_s+*(+,-)',nymx,ymin,ymax,0.0)
 
        call hbook1(343,'dN/deta, B+(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(3343,'dN/deta, resB+*(+,-)',nymx,ymin,ymax,0.0)

        call hbook1(344,'dN/deta, B0(+,-)',nymx,ymin,ymax,0.0)
         call hbook1(3344,'dN/deta, resB0*(+,-)',nymx,ymin,ymax,0.0)

        call hbook1(345,'dN/deta, B_s0(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(3345,'dN/deta, resB_s0*(+,-)',nymx,ymin,ymax,0.0)
        
        call hbook1(346,'dN/deta, Lc+(+,-)',nymx,ymin,ymax,0.0) 
        call hbook1(347,'dN/deta, Bc+(+,-)',nymx,ymin,ymax,0.0)
        call hbook1(3347,'dN/deta, resBc+*(+,-)',nymx,ymin,ymax,0.0)

        call hbook1(348,'dN/deta, J/PSI',nymx,ymin,ymax,0.0)
         
        call hbook1(349,'dN/deta, Upsilon',nymx,ymin,ymax,0.0)

c       Create the pseudo-rapidity distributions of the baryons.  
	call hbook1(350,'dN/deta, p',nymx,ymin,ymax,0.0)
	call hbook1(351,'dN/deta, p~',nymx,ymin,ymax,0.0)
	call hbook1(352,'dN/deta, n',nymx,ymin,ymax,0.0)
	call hbook1(353,'dN/deta, n~',nymx,ymin,ymax,0.0)
	call hbook1(354,'dN/deta, Lambda',nymx,ymin,ymax,0.0)
	call hbook1(355,'dN/deta, Lambda~',nymx,ymin,ymax,0.0)
        call hbook1(356,'dN/deta, Sigma-',nymx,ymin,ymax,0.0)
        call hbook1(357,'dN/deta, Sigma~-',nymx,ymin,ymax,0.0)
        call hbook1(358,'dN/deta, Sigma0',nymx,ymin,ymax,0.0)
        call hbook1(359,'dN/deta, Sigma~0',nymx,ymin,ymax,0.0)
        call hbook1(360,'dN/deta, Sigma+',nymx,ymin,ymax,0.0)
        call hbook1(361,'dN/deta, Sigma~-',nymx,ymin,ymax,0.0)
        call hbook1(362,'dN/deta, Xi-',nymx,ymin,ymax,0.0)
        call hbook1(363,'dN/deta, Xi~+',nymx,ymin,ymax,0.0)
        call hbook1(364,'dN/deta, Xi0',nymx,ymin,ymax,0.0)
        call hbook1(365,'dN/deta, Xi0~',nymx,ymin,ymax,0.0)
        call hbook1(366,'dN/deta, Omega-',nymx,ymin,ymax,0.0)
        call hbook1(367,'dN/deta, Omega~+',nymx,ymin,ymax,0.0)
c
c       Create the det/dy distributions
	call hbook1(400,'det/dy',nymx,ymin,ymax,0.0)
	call hbook1(401,'det/dy char',nymx,ymin,ymax,0.0)
	call hbook1(402,'det/dy neut',nymx,ymin,ymax,0.0)
        call hbook1(403,'det/dy ch+',nymx,ymin,ymax,0.0)
        call hbook1(404,'det/dy ch-',nymx,ymin,ymax,0.0)



c        Create the det/deta distributions
          call hbook1(500,'det/deta',nymx,ymin,ymax,0.0)
          call hbook1(501,'det/deta char',nymx,ymin,ymax,0.0)
          call hbook1(502,'det/deta neut',nymx,ymin,ymax,0.0)
          call hbook1(503,'det/deta ch+',nymx,ymin,ymax,0.0)
          call hbook1(504,'det/deta ch-',nymx,ymin,ymax,0.0)
c\end topor
c         
	w_event=1.0/real(n_event)
	
	r_rms=0.0
	ptot=0.0
	detot=0.0
	antot=0.0
	n_ev0=0
	i_out=0
c      
c       set the initial number of the various particles species to 0. 
c
	ttot = 0.0
       
	tchar = 0.0
        tcharn=0.0
        tcharp=0.0 
	tneut = 0.0
	tpi0 = 0.0
	tpip = 0.0
	tpim = 0.0
	tKp  = 0.0
	tKm  = 0.0
	tK0  = 0.0
	taK0 = 0.0
	tK0S = 0.0
	tKstar = 0.0
	tD   = 0.0
        tDres = 0.0
        tD0  = 0.0
        tD0res = 0.0
        tDs  = 0.0
        tDsres = 0.0
        tlc = 0.0
        tB   = 0.0
        tBres = 0.0
        tB0  = 0.0
        tB0res =0.0 
        tBs  = 0.0 
        tBsres = 0.0
        tBc = 0.0
        tBcres = 0.0 
        tJpsi = 0.0
        tUpsi =0.0
	tp   = 0.0
	tap  = 0.0
	tn   = 0.0
	tan  = 0.0
	tlam = 0.0
	talam = 0.0
	tsig0  = 0.0
	tasig0 = 0.0
	tsigm = 0.0
	tasigp = 0.0	
	tsigp = 0.0
	tasigm = 0.0	
	txi0 = 0.0
	taxi0 = 0.0
	txim = 0.0
	taxip = 0.0
	tomem = 0.0
	taomep = 0.0

c       set up the weights for the ET distributions.  
	call hijing(frame,bmin,bmax)
	      
	avn=(1.97+0.42*alog(hint1(1))+0.592*(alog(hint1(1)))**2)*4.0
	if(ihnt2(1).gt.1.or.ihnt2(3).gt.1 ) avn=natt-ihnt2(1)-ihnt2(3)
	etmax=avn*(0.60-0.030*alog(hint1(1))+
     1    0.0053*(alog(hint1(1)))**2)*1.5
	etmax=max(40.0,etmax)
	avn=max(40.0,avn)
	w_et=2.0
	w_an=2.0
	netmx=etmax/w_et
	navmx=avn/w_an
	w_et=etmax/real(netmx)
	w_an=avn/real(navmx)
	
c       record various properties of the event. 
	call hbook1(70,'P(n)',navmx,0.0,avn,0.)
	call hbook1(71,'P(n),jet=0',navmx,0.0,avn,0.)
	call hbook1(72,'P(n),jet=1',navmx,0.0,avn,0.)
	call hbook1(73,'P(n),jet>1',navmx,0.0,avn,0.)
c
	call hbook1(81,'dsigma/dEt',netmx,0.0,etmax,0.0)
	call hidopt(81,'logy')
	call hbook1(82,'pt(eta_cut) vs nch',navmx,0.0,avn,0.)
	call hbook1(83,'nch(eta_cut) vs nch',navmx,0.0,avn,0.)



c
	do 2000 i_event=1,n_event
c
	   if(iap.lt.10) icount = 1000
	   if(iap.ge.10) icount = 100
	   if(mod(i_event,icount).eq.0) write(6,*) ' i_out = ',i_event
c       
	   nch_tot=0
	   nch_cor=0
	   pt_cor=0.0
	   ettot=0.0
	   vetop=0.0
	   vetot=0.0
	   
 100	   continue
	   
	   call hijing(frame,bmin,bmax) 
	   ncol = n0	   
	   if (ncol.gt.0) then
	      nhit = nhit + 1
	   else if (ncol.le.0) then
	      nmiss = nmiss + 1
	      goto 100
	   end if 
	   
	   delrr = 0.2	   
c       record the initial distributions of the nucleons. 
	   do 120 i=1, iap
	      rp=sqrt(yp(1,i)**2+yp(2,i)**2+yp(3,i)**2)
	      nbin=rp/delrr 
	      rbin=0.5*delrr + nbin*delrr
	      call hfill(55,rp,0.0,w_event/(12.566*delrr*rbin**2))
	      call hfill(58,rp,0.0,1.)
	      r_rms=r_rms+rp**2
 120	   continue
	   do 121 i=1, IAT
	      RT=SQRT(YT(1,I)**2+YT(2,I)**2+YT(3,I)**2)
	      nbin=rt/delrr 
	      rbin=0.5*delrr + nbin*delrr
	      call hfill(56,rt,0.0,w_event/(12.566*delrr*rbin**2))
	      call hfill(59,rt,0.0,1.)
	      r_rms=r_rms+rt**2
 121	   CONTINUE
c
	   call hfill(57,real(n0+n01+n10+n11),0.0,w_event)
C       topor histo for nb -- baryon junction and JBARJ 

C	   call hfill(25,real(nb),0.0,w_event)
C           call hfill(31,real(ncb(i)),0.0,w_event)
c           call hfill(32,real(ncb(2)),0.0,w_event)
c           call hfill(33,real(ncb(3)),0.0,w_event)
c           call hfill(34,real(ncb(4)),0.0,w_event)
c       compute the number of charge particles per event. 
	   do 150 i=1,natt
	      if(luchge(katt(i,1)).ne.0) nch_tot=nch_tot+1
 150	   continue
	   nch_prod=nch_tot-izp-izt 
	   
c       loop through the particles produced from one event. 
	   do 1000 i=1,natt
	      
c       determine the kinematic properties of each particle. 
	      xulm=ulmass(katt(i,1))
	      ptr2=patt(i,1)**2+patt(i,2)**2
	      xulmt=sqrt(ptr2 + xulm**2)
         
	      
	      rap=100.
	      if (xulmt.gt.0.) then
		 if (patt(i,3).gt.0.) then
		    rap=alog((patt(i,4)+patt(i,3))/xulmt)
		 else
		    rap=alog(xulmt/(patt(i,4)-patt(i,3)) )
		 end if
	      end if
	      
c       determine the veto energies for the target and the projectile.
	      if(rap.gt.hint1(4)-0.02) then
		 vetop=vetop+patt(i,4)
	      else if(rap.lt.hint1(5)+0.02) then
		 vetot=vetot+patt(i,4)
	      endif
c       
c       skip over those particles which have not interacted. 
	      if(katt(i,2).eq.0.or.katt(i,2).eq.10) go to 1000
	      
	      ptr=sqrt(ptr2)
	      pab=sqrt(patt(i,3)**2+ptr2)
	      et=patt(i,4)*ptr/max(pab,0.0000001)
	      phi=ulangl(patt(i,1),patt(i,2))
	      theta=acos(patt(i,3)/pab)
	      
c       determine eta
	      if (patt(i,3).gt.0.) then
		 etann=pab+patt(i,3)
		 etadd=ptr
	      else
		 etann=ptr
		 etadd=pab-patt(i,3)
	      end if
	      
	      eta=100.
	      if (etadd.gt.0.) eta=alog(etann/etadd)
	      
c       Note that the following scheme breaks down for large eta
c	etann=max(abs(pab+patt(i,3)),0.000001)
c       etadd=max(abs(pab-patt(i,3)),0.000001)
c       eta=0.5*alog(etann/etadd)
	      
	   
c       determine xf and xe
	      xf=2.0*sqrt(abs(patt(i,4)**2-patt(i,3)**2))*
     1	   sinh(rap-hint1(3))/hint1(1)
	      xe=2.0*sqrt(abs(patt(i,4)**2-patt(i,3)**2))*
     1	   cosh(rap-hint1(3))/hint1(1)
	      
c       
c       record histograms for the charged, neutral and total
c       
	      call hfill(400,rap,0.0,et*w_event/w_y)
c topor  dE_t/deta 
           
              call hfill(500,eta,0.0,et*w_event/w_y)
        	      ttot = ttot + w_event
                     ptot=ptot+ptr     
	      if (luchge(katt(i,1)).eq.0) then
                 call hfill(202,rap,0.0,w_event/w_y)
		 call hfill(302,eta,0.0,w_event/w_y)
                 call hfill(402,rap,0.0,et*w_event/w_y)
                 call hfill(502,eta,0.0,et*w_event/w_y)
		 tneut = tneut + w_event
	      else
		 if(abs(eta).le.eta_cut) then
		    call hfill(101,ptr,0.0,w_event/2.0/max(ptr,0.00001))
		    nch_cor=nch_cor+1
		    pt_cor=pt_cor+ptr
		 endif		 
		 call hfill(201,rap,0.0,w_event/w_y)
		 call hfill(301,eta,0.0,w_event/w_y)
		 call hfill(401,rap,0.0,et*w_event/w_y)
                 call hfill(501,eta,0.0,et*w_event/w_y)
                 tchar = tchar + w_event
              
	      end if	

	      

c topor      ch-      

               if (luchge(katt(i,1)).lt.0) then
           
                 call hfill(204,rap,0.0,w_event/w_y)
		 call hfill(304,eta,0.0,w_event/w_y)
                 call hfill(404,rap,0.0,et*w_event/w_y)
                 call hfill(504,eta,0.0,et*w_event/w_y)
                     tcharn = tcharn + w_event       
             
                 
              
            if(abs(eta).le.eta_cut.and.luchge(katt(i,1)).lt.0) then 
            call hfill(104,ptr,0.0,w_event/2.0/max(ptr,0.0000001))
                endif
                endif 
                
                                    

c topor      ch+      

                 if (luchge(katt(i,1)).gt.0) then
                 call hfill(203,rap,0.0,w_event/w_y)
		 call hfill(303,eta,0.0,w_event/w_y)
                 call hfill(403,rap,0.0,et*w_event/w_y)
                 call hfill(503,eta,0.0,et*w_event/w_y)
                tcharp = tcharp + w_event
               
              

                 
           if(abs(eta).le.eta_cut.and.luchge(katt(i,1)).gt.0) then 
    		call hfill(103,ptr,0.0,w_event/2.0/max(ptr,0.0000001))
                endif
                 endif 
                   
c       
c       
c       record the gammas:
c 
	      if (katt(i,1).eq.22) then
c       gammas from all resonances except the pi0 decay !!
		 call hfill(215,rap,0.0,w_event/w_y)
		 call hfill(315,eta,0.0,w_event/w_y)
		 if (katt(i,2).ne.40) then       
		    if (abs(rap).le.eta_cut) call hfill(115,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))
		 else
c       gammas from direct QCD proc
		    if(abs(rap).le.eta_cut) call hfill(116,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))
		 endif
 
            
 
                                              



c       
c       pi0
	      else if (abs(katt(i,1)).lt.1000) then
		 if (katt(i,1).eq.111) then
                    call hfill(217,rap,0.0,w_event/w_y)
		    call hfill(317,eta,0.0,w_event/w_y)
                     tpi0 = tpi0 + w_event
 		    if (abs(rap).le.eta_cut) call hfill(117,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))



		    
c       pi+
		 else if (katt(i,1).eq.211) then
		    call hfill(320,eta,0.0,w_event/w_y)
		    call hfill(220,rap,0.0,w_event/w_y)
                   tpip = tpip + w_event


		    if(abs(rap).le.eta_cut) call hfill(120,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

		    
c       pi-
		 else if(katt(i,1).eq.-211) then
		    call hfill(321,eta,0.0,w_event/w_y)
		    call hfill(221,rap,0.0,w_event/w_y)
                  tpim = tpim + w_event

   		    if(abs(rap).le.eta_cut) call hfill(121,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))



C       



c       greekro0
		 else if(katt(i,1).eq.113) then

		    call hfill(218,rap,0.0,w_event/w_y)

   		    if(abs(rap).le.eta_cut) call hfill(118,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))


c       greekomega
		 else if(katt(i,1).eq.223) then

		    call hfill(1218,rap,0.0,w_event/w_y)

   		    if(abs(rap).le.eta_cut) call hfill(1118,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))


c       greeketa
		 else if(katt(i,1).eq.221) then

		    call hfill(219,rap,0.0,w_event/w_y)

   		    if(abs(rap).le.eta_cut) call hfill(119,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))


c       greeketaprim
		 else if(katt(i,1).eq.331) then

		    call hfill(1219,rap,0.0,w_event/w_y)

   		    if(abs(rap).le.eta_cut) call hfill(1119,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))


c       greekphi
		 else if(katt(i,1).eq.333) then

		    call hfill(2233,rap,0.0,w_event/w_y)

   		    if(abs(rap).le.eta_cut) call hfill(1133,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))



		    
		
c       Ks^0
		 else if(katt(i,1).eq.310) then
		    call hfill(226,rap,0.0,w_event/w_y)
                    call hfill(326,eta,0.0,w_event/w_y)
		    tK0 = tK0 + w_event


		    if(abs(rap).le.eta_cut) call hfill(126,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.00001))


c       Ks^0~0
		 else if(katt(i,1).eq.-310) then
		    call hfill(227,rap,0.0,w_event/w_y)
                     call hfill(327,eta,0.0,w_event/w_y)
		    taK0 = taK0 + w_event

     		    if(abs(rap).le.eta_cut) call hfill(127,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.00001))              



c       K+
		 else if(katt(i,1).eq.321) then
		    call hfill(322,eta,0.0,w_event/w_y)
		    call hfill(222,rap,0.0,w_event/w_y)
                      tKp = tKp + w_event

		    if (abs(rap).le.eta_cut) call hfill(122,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.00001))

		   
c       K-
		 else if(katt(i,1).eq.-321) then
		    call hfill(323,eta,0.0,w_event/w_y)
		    call hfill(223,rap,0.0,w_event/w_y)
                  tKm = tKm + w_event

		    if (abs(rap).le.eta_cut) call hfill(123,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.00001))

		   


c       K0
		 else if(katt(i,1).eq.311) then
                    call hfill(324,eta,0.0,w_event/w_y)
		    call hfill(224,rap,0.0,w_event/w_y)
		    tK0 = tK0 + w_event
		     
		    
c       K*0
		 else if(katt(i,1).eq.313) then
                    call hfill(328,rap,0.0,w_event/w_y)
		    call hfill(228,rap,0.0,w_event/w_y)
		    tKstar = tKstar + w_event

c       D+(+,-)
		 else if(abs(katt(i,1)).eq.411) then
		    call hfill(340,eta,0.0,w_event/w_y)
		    call hfill(240,rap,0.0,w_event/w_y)
                    tD = tD + w_event
		    if(abs(rap).le.eta_cut) call hfill(140,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       resD+*(+,-)
		 else if(abs(katt(i,1)).eq.413) then
		    call hfill(3340,eta,0.0,w_event/w_y)
		    call hfill(2240,rap,0.0,w_event/w_y)
                    tDres = tDres + w_event
		    if(abs(rap).le.eta_cut) call hfill(1140,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))
       

c       D0(+,-)
		 else if(abs(katt(i,1)).eq.421) then
		    call hfill(341,eta,0.0,w_event/w_y)
		    call hfill(241,rap,0.0,w_event/w_y)
                    tD0 = tD0 + w_event
		    if(abs(rap).le.eta_cut) call hfill(141,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       resD0*(+,-)
		 else if(abs(katt(i,1)).eq.423) then
		    call hfill(3341,eta,0.0,w_event/w_y)
		    call hfill(2241,rap,0.0,w_event/w_y)
                    tD0res = tD0res + w_event
		    if(abs(rap).le.eta_cut) call hfill(1141,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))


c       D_s+(+,-)
		 else if(abs(katt(i,1)).eq.431) then
		    call hfill(342,eta,0.0,w_event/w_y)
		    call hfill(242,rap,0.0,w_event/w_y)
                    tDs = tDs + w_event
		    if(abs(rap).le.eta_cut) call hfill(142,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       resD_s+*(+,-)
		 else if(abs(katt(i,1)).eq.433) then
		    call hfill(3342,eta,0.0,w_event/w_y)
		    call hfill(2242,rap,0.0,w_event/w_y)
                    tDsres = tDsres + w_event
		    if(abs(rap).le.eta_cut) call hfill(1142,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))


c       B+(+,-)
		 else if(abs(katt(i,1)).eq.521) then
		    call hfill(343,eta,0.0,w_event/w_y)
		    call hfill(243,rap,0.0,w_event/w_y)
                    tB = tB + w_event
		    if(abs(rap).le.eta_cut) call hfill(143,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))     

c       resB+*(+,-)
		 else if(abs(katt(i,1)).eq.523) then
		    call hfill(3343,eta,0.0,w_event/w_y)
		    call hfill(2243,rap,0.0,w_event/w_y)
                    tBres = tBres + w_event
		    if(abs(rap).le.eta_cut) call hfill(1143,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001)) 


c       B0(+,-)
		 else if(abs(katt(i,1)).eq.511) then
		    call hfill(344,eta,0.0,w_event/w_y)
		    call hfill(244,rap,0.0,w_event/w_y)
                    tB0 = tB0 + w_event
		    if(abs(rap).le.eta_cut) call hfill(144,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))  

c       resB0*(+,-)
		 else if(abs(katt(i,1)).eq.513) then
		    call hfill(3344,eta,0.0,w_event/w_y)
		    call hfill(2244,rap,0.0,w_event/w_y)
                    tB0res = tB0res + w_event
		    if(abs(rap).le.eta_cut) call hfill(1144,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))  
   

c       B_s0(+,-)
		 else if(abs(katt(i,1)).eq.531) then
		    call hfill(345,eta,0.0,w_event/w_y)
		    call hfill(245,rap,0.0,w_event/w_y)
                    tBs = tBs + w_event
		    if(abs(rap).le.eta_cut) call hfill(145,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       resB_s0*(+,-)
		 else if(abs(katt(i,1)).eq.533) then
		    call hfill(3345,eta,0.0,w_event/w_y)
		    call hfill(2245,rap,0.0,w_event/w_y)
                    tBsres = tBsres + w_event
		    if(abs(rap).le.eta_cut) call hfill(1145,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       B_c+(+,-)
		 else if(abs(katt(i,1)).eq.541) then
		    call hfill(347,eta,0.0,w_event/w_y)
		    call hfill(247,rap,0.0,w_event/w_y)
                    tBc = tBc + w_event
		    if(abs(rap).le.eta_cut) call hfill(147,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       resB_c+*(+,-)
		 else if(abs(katt(i,1)).eq.543) then
		    call hfill(3347,eta,0.0,w_event/w_y)
		    call hfill(2247,rap,0.0,w_event/w_y)
                    tBcres = tBcres + w_event
		    if(abs(rap).le.eta_cut) call hfill(1147,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       J/PSI(+,-)
		 else if(abs(katt(i,1)).eq.443) then
		    call hfill(348,eta,0.0,w_event/w_y)
		    call hfill(248,rap,0.0,w_event/w_y)
                    tJpsi = tJpsi + w_event
		    if(abs(rap).le.eta_cut) call hfill(148,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))

c       Upsilon(+,-)
		 else if(abs(katt(i,1)).eq.553) then
		    call hfill(349,eta,0.0,w_event/w_y)
		    call hfill(249,rap,0.0,w_event/w_y)
                    tUpsi = tUpsi + w_event
		    if(abs(rap).le.eta_cut) call hfill(148,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000000001))


		 end if 

c       
c       
c       record the baryons. 
c      
	      else if (abs(katt(i,1)).gt.1000) then
c       p
		 if(katt(i,1).eq.2212) then
		    call hfill(350,eta,0.0,w_event/w_y)
		    call hfill(250,rap,0.0,w_event/w_y)
                 tp = tp + w_event

		    if(abs(rap).le.eta_cut) call hfill(150,ptr,0.0,
     1		      w_event/2.0/max(ptr,0.0000001))

		   
c       p~
		 else if(katt(i,1).eq.-2212) then
		    call hfill(351,eta,0.0,w_event/w_y)
		    call hfill(251,rap,0.0,w_event/w_y)
              tap = tap + w_event

		    if(abs(rap).le.eta_cut) call hfill(151,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

		    
c       n
		 else if (katt(i,1).eq.2112) then
		    call hfill(252,rap,0.0,w_event/w_y)
                    call hfill(352,eta,0.0,w_event/w_y)
		    tn = tn + w_event
c       n~
		 else if (katt(i,1).eq.-2112) then
		    call hfill(253,rap,0.0,w_event/w_y)
                    call hfill(353,eta,0.0,w_event/w_y)
		    tan = tan + w_event
c       Lambda
		 else if (katt(i,1).eq.3122) then
		    call hfill(254,rap,0.0,w_event/w_y)
                    call hfill(354,eta,0.0,w_event/w_y)
                    tlam = tlam + w_event

		    if(abs(rap).le.eta_cut) call hfill(154,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

		    
c       Lambda~
		 else if (katt(i,1).eq.-3122) then
		    call hfill(255,rap,0.0,w_event/w_y)
                    call hfill(355,rap,0.0,w_event/w_y)
		 talam = talam + w_event 

		    if(abs(rap).le.eta_cut) call hfill(155,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

c       Lambda_C(+,-)
		 else if (abs(katt(i,1)).eq.4122) then
		    call hfill(246,rap,0.0,w_event/w_y)
                    call hfill(346,rap,0.0,w_event/w_y)
		 tlc = tlc + w_event 

		    if(abs(rap).le.eta_cut) call hfill(146,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))


c       Sigma-
		 else if (katt(i,1).eq.3112) then
		    call hfill(256,rap,0.0,w_event/w_y)
                    call hfill(356,eta,0.0,w_event/w_y)
             		    tsigm = tsigm + w_event

		    if(abs(rap).le.eta_cut) call hfill(156,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

               

c       Sigma~+
		 else if (katt(i,1).eq.-3112) then
		    call hfill(257,rap,0.0,w_event/w_y)
                     call hfill(357,eta,0.0,w_event/w_y)
		    tasigp = tasigp + w_event


		    if(abs(rap).le.eta_cut) call hfill(157,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

                       

c       Sigma0
		 else if (katt(i,1).eq.3212) then
		    call hfill(258,rap,0.0,w_event/w_y)
                     call hfill(358,eta,0.0,w_event/w_y)
		    tsig0 = tsig0 + w_event

		    if(abs(rap).le.eta_cut) call hfill(158,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

                       
                 
c       Sigma~0
		 else if (katt(i,1).eq.-3212) then
		    call hfill(259,rap,0.0,w_event/w_y)
                      call hfill(359,eta,0.0,w_event/w_y) 
		    tasig0 = tasig0 + w_event

		    if(abs(rap).le.eta_cut) call hfill(159,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

                       


c       Sigma+
		 else if (katt(i,1).eq.3222) then
		    call hfill(260,rap,0.0,w_event/w_y)
                    call hfill(360,eta,0.0,w_event/w_y)
		    tsigp = tsigp + w_event

		    if(abs(rap).le.eta_cut) call hfill(160,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

c       Sigma~-
		 else if (katt(i,1).eq.-3222) then
		    call hfill(261,rap,0.0,w_event/w_y)
                   call hfill(361,eta,0.0,w_event/w_y)
		    tasigm = tasigm + w_event

		    if(abs(rap).le.eta_cut) call hfill(161,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))


		    
		    
c       Xi-
		 else if (katt(i,1).eq.3312) then
		    call hfill(262,rap,0.0,w_event/w_y)
                    call hfill(362,eta,0.0,w_event/w_y)
                    txim = txim + w_event

                  if(abs(rap).le.eta_cut) call hfill(162,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

		    
c       Xi~+
		 else if (katt(i,1).eq.-3312) then
		    call hfill(263,rap,0.0,w_event/w_y)
                     call hfill(363,eta,0.0,w_event/w_y)
                     taxip = taxip + w_event

		    if(abs(rap).le.eta_cut) call hfill(163,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

		   
c       Xi0
		 else if (katt(i,1).eq.3322) then
		    call hfill(264,rap,0.0,w_event/w_y)
                    call hfill(364,eta,0.0,w_event/w_y)
		    txi0 = txi0 + w_event

		    if(abs(rap).le.eta_cut) call hfill(164,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))



c       Xi~0
		 else if (katt(i,1).eq.-3322) then
		    call hfill(265,rap,0.0,w_event/w_y)
                    call hfill(365,eta,0.0,w_event/w_y)
		    taxi0 = taxi0 + w_event

		    if(abs(rap).le.eta_cut) call hfill(165,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))



c       Omega-
		 else if (katt(i,1).eq.3334) then
		    call hfill(266,rap,0.0,w_event/w_y)
                    call hfill(366,eta, 0.0,w_event/w_y)
                       tomem = tomem + w_event

		    if(abs(rap).le.eta_cut) call hfill(166,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))

		   
c       Omega~+
		 else if (katt(i,1).eq.-3334) then
		    call hfill(267,rap,0.0,w_event/w_y)
                     call hfill(367,eta,0.0,w_event/w_y)
                    taomep = taomep + w_event	 

		    if(abs(rap).le.eta_cut) call hfill(167,ptr,0.0,
     1		 w_event/2.0/max(ptr,0.0000001))




		 end if

	      end if 


 1000	   continue

c       record the totals of the event
	   call hfill(81,ettot,0.0,w_event/w_et)
	   call hfill(82,real(nch_cor),0.0,pt_cor*w_event/w_an)
	   call hfill(83,real(nch_cor),0.0,nch_cor*w_event/w_an)
	   xn=nch_tot-0.5
	   call hfill(70,xn,0.0,w_event/w_an)
	   detot=detot+eatt-hint1(1)*(ihnt2(1)+ihnt2(3))/2.0
	   antot=antot+nch_tot*w_event
	   if(ihpr2(8).ne.0) then
	      if(jatt.eq.0) call hfill(71,xn,0.0,w_event/w_an)
	      if(jatt.eq.1) call hfill(72,xn,0.0,w_event/w_an)
	      if(jatt.gt.1) call hfill(73,xn,0.0,w_event/w_an)
	   end if 
	   
c	record the geometry of the event.
	   if(ihnt2(1).gt.1) then
	      call hfill(50,real(np),0.0,w_event)
	      call hfill(52,real(np+nt),0.0,w_event)
	      call hfill(53,vetop/hint1(6),0.0,w_event)
	   endif
	   if(ihnt2(3).gt.1) then
	      call hfill(51,real(nt),0.0,w_event)
	      call hfill(54,vetot/hint1(7),0.0,w_event)
	   endif
C       
	   an_cut=an_cut+nch_cor*w_event
	   pt_cut_tot=pt_cut_tot+pt_cor
	
c       record the parton statistics
	   
c       partons embeddedin the projectile
	   do 1200 i=1,iap
	      do 1210 j=1,npj(i)
c	gluon rap distrib
		 if (abs(kfpj(i,j)).eq.21) then
		    call hfill(205,rap,0.0,w_event/w_y)
		 end if
		 
c	q+qbar rap distrib
		 if (abs(kfpj(i,j)).le.6.and.abs(kfpj(i,j)).ge.1) then
		    call hfill(206,rap,0.0,w_event/w_y)
		 end if
		 
c	qq+qqbar rap distrib
		 if (abs(kfpj(i,j)).gt.1000) then
		    call hfill(207,rap,0.0,w_event/w_y)
		 end if
		 
c	c+cbar rap distrib
		 if (abs(kfpj(i,j)).eq.4) then
		    call hfill(210,rap,0.0,w_event/w_y)
		 end if
		 
 1210	      continue
 1200	   continue
	   
c       partons embedded in the target
	   do 1300 i=1,iat
	      do 1310, j=1,ntj(i)
c	gluon rap distrib
		 if (abs(kftj(i,j)).eq.21) then
		    call hfill(205,rap,0.0,w_event/w_y)
		 end if
		 
c	q+qbar rap distrib
		 if (abs(kftj(i,j)).le.6.and.abs(kftj(i,j)).ge.1) then
		    call hfill(206,rap,0.0,w_event/w_y)
		 end if
		 
c	qq+qqbar rap distrib
		 if (abs(kftj(i,j)).gt.1000) then
		    call hfill(207,rap,0.0,w_event/w_y)
		 end if
		 
c	c+cbar rap distrib
		 if (abs(kftj(i,j)).eq.4) then
		    call hfill(210,rap,0.0,w_event/w_y)
		 end if
	      
 1310	      continue
 1300	   continue
	   
c       partons from the isolated non nucleon strings
	   
	   do 1400 i=1,nsg
	      do 1410 j=1,njsg(i)

c	gluon rap distrib
		 if (abs(k2sg(i,j)).eq.21) then
		    call hfill(205,rap,0.0,w_event/w_y)
		 end if 
		 
c	q+qbar rap distrib
		 if (abs(k2sg(i,j)).le.6.and.abs(k2sg(i,j)).ge.1) then
		    call hfill(206,rap,0.0,w_event/w_y)
		 end if
c	qq+qqbar rap distrib
		 if (abs(k2sg(i,j)).gt.1000) then
		    call hfill(207,rap,0.0,w_event/w_y)
		 end if
c	c+cbar rap distrib
		 if (abs(k2sg(i,j)).eq.4) then
		    call hfill(210,rap,0.0,w_event/w_y)
		 end if
 1410	      continue
 1400	   continue
c       
	   
 2000	continue
	
c       
c \beg topor  include the bin width of the pT histograms

       

	id(1) = 101 
        id(2) = 103
        id(3) = 104
        id(4) = 115
        id(5) = 116
	id(6) = 117
	id(7) = 120
	id(8) = 121
	id(9) = 122
	id(10) =123
        id(11) =126
        id(12) =127
	id(13) =140
        id(14) =141
        id(15) =142 
        id(16) =143
        id(17) =144
        id(18) =145

	id(19) =150
	id(20) =151
        id(21) =154
        id(22) =155
        id(23) =156
        id(24) =157
        id(25) =158
        id(26) =159
        id(27) =160
        id(28) =161 
        id(29) =162
        id(30) =163
        id(31) =164
        id(32) =165
        id(33) =166
        id(34) =167
        id(35) =118
        id(36) =119
        id(37) =118
        id(38) =1118
        id(39) =1119
        id(40) =1133
        
        id(41) =1140
        id(42) =1141
        id(43) =1142
        id(44) =1143
        id(45) =1144
        id(46) =1145
        id(47) =146
        id(48) =1147
        id(49) =147
        id(50) =1147
        id(51) =148
        id(52) =149    
	jjmax = 52
	

        call hbookb(20,'DPT',n_ptmx,pthbk,0.0)
	call hidopt(20,'logy')
	do 2500 i=1,n_ptmx
	   ptbin=(pthbk(i+1)+pthbk(i))/2.0
	   dptbin=pthbk(i+1)-pthbk(i)
	   call hfill(20,ptbin,0.0,1.0/dptbin)
 2500	continue
	do 2700 i=1,jjmax
	   call hopera(id(i),'*',20,id(i),1.0,1.0)
 2700	continue

	do 3200 jj=1,jjmax
	   call hgive(id(jj),chtit,nx,xmi,xma,ny,ymi,yma,nwt,loc)
	   sum=hsum(id(jj))
	   call hrebin(id(jj),xhb,yhb,exhb,eyhb,nx,1,nx)	
	   do 3098 i=1,nx
	      nnx=nx-i+1
 3098	      if (yhb(nnx).ne.0.) go to 3099
 3099	      continue
	      do 3100 i=1,nnx
		 xxbin=xhb(i)
		 xxbin=(PTHBK(I+1)+PTHBK(I))/2.0
 3100	      continue
 3200	   continue
	   
c \end topor
	
C*******************************************
 4001	FORMAT(1X,I10,3X,'Ntot=',I10,3X,'Nch=',I10,
     &	3X,'Njet=',I10/1X,'***************************'/)
	
	PTOT=PTOT*W_EVENT/ANTOT
	PT_CUT=PT_CUT_TOT*W_EVENT/AN_CUT
	
	R_RMS=SQRT(R_RMS/FLOAT(IAP+IAT)*W_EVENT)
	write(6,*) ' geometry rms radius av =',r_rms
	write(6,*) ' DETOT =',DETOT*w_event,' ANTOT =',ANTOT,' <PT> =',PTOT
	write(6,*) 'AN_CUT=', AN_CUT,'  PT_CUT=',PT_CUT
	write(6,*) 'E_CM=',HINT1(1),' I_DIPOLE=',IHPR2(1),' I_RAD=',IHPR2(2)
	write(6,*) 'Q_CUT=',HIPR1(5),' PTRAD_CUT=',HIDAT(2),' P_RAD=',HIDAT(3)
	write(6,*) 'MCUT=',HIPR1(1),' PT_INIRAD',HIDAT(1),' CRS_SD=',HIDAT(4)
	write(6,*) 'PTRMS=',HIPR1(2),' a=',PARJ(41),' b=',PARJ(42)
	write(6,*) ' JMAX=',IHPR2(8),' HINT1(11)',HINT1(11),'Stot',HINT1(13)
	write(6,*) 'HINT1(10) =',HINT1(10),'Sinel =',HINT1(12)
	write(6,*) 'S_Trig =',HINT1(60),' S_trig_t =',HINT1(59)
        
        write(6,*) 'HINT1(58) =', HINT1(58) 
        write(6,*) 'HINT1(59) =', HINT1(59)
        write(6,*) 'HINT1(60) =', HINT1(60)
        write(6,*) 'HINT1(61) =', HINT1(61)
        write(6,*) 'HINT1(62) =', HINT1(62)
	write(6,*) 'HINT1(63) =', HINT1(63)
        write(6,*) 'HINT1(64) =', HINT1(64)
        write(6,*) 'HINT1(65) =', HINT1(65)



	write(6,*) ' '
	write(6,*) ' NHIT= ',NHIT,' OUT OF NEVENT =', nhit+nmiss
	write(6,*) ' '
	
c       reaction cross section in mb
	sigreac=3.1415*(bmax*bmax-bmin*bmin)*10.*float(nhit)/float(n_event)
	write(6,*) ' reaction cross section mb = ',sigreac
	write(6,*) ' '

	write(23,*) 'Total numbers of the various particles: '
	write(23,*) ' '
	write(23,*) 'total   =  ',ttot	
	write(23,*) 'charged =  ',tchar
	write(23,*) 'neutral =  ',tneut
	write(23,*) ' ' 
	write(23,*) 'h-      =  ',tcharn
        write(23,*) 'h+      = ',tcharp
	write(23,*) ' '
	write(23,*) ' **** Mesons **** '
	write(23,*) 'pi0      = ',tpi0
	write(23,*) 'pi+      = ',tpip
	write(23,*) 'pi-      = ',tpim
	write(23,*) 'K0       = ',tK0
	write(23,*) 'K~0      = ',taK0
	write(23,*) 'K+       = ',tKp
	write(23,*) 'K-       = ',tKm
	write(23,*) ' **** Baryons **** '
	write(23,*) 'p        = ',tp
	write(23,*) 'p~       = ',tap
	write(23,*) 'n        = ',tn
	write(23,*) 'n~       = ',tan
	write(23,*) 'Lambda   = ',tlam
	write(23,*) 'Lambda~  = ',talam
	write(23,*) 'Sigma-   = ',tsigm
	write(23,*) 'Sigma~+  = ',tasigp
	write(23,*) 'Sigma0   = ',tsig0
	write(23,*) 'Sigma~0  = ',tasig0
	write(23,*) 'Sigma+   = ',tsigp
	write(23,*) 'Sigma~-  = ',tasigm
	write(23,*) 'Xi-      = ',txim
	write(23,*) 'Xi~+     = ',taxip
	write(23,*) 'Xi0      = ',txi0
	write(23,*) 'Xi~0     = ',taxi0
	write(23,*) 'Omega-   = ',tomem
	write(23,*) 'Omega~+  = ',taomep
c - charm histo
        write(23,*) 'D+(+-)  = ',tD
        write(23,*) 'resD+*(+-)  = ',tDres
        write(23,*) 'D0(+-)  = ',tD0
        write(23,*) 'resD0*(+-)  = ',tD0res
        write(23,*) 'Ds+(+-)  = ',tDs
        write(23,*) 'resDs+*(+-)  = ',tDsres
        write(23,*) 'lambda_c(+-) = ',tlc
        write(23,*) 'B+(+-)  = ',tB
        write(23,*) 'resB+*(+-)  = ',tBres
        
        write(23,*) 'B0(+-)  = ',tB0
        
        write(23,*) 'resB0*(+-)  = ',tB0res
        write(23,*) 'B_c+(+-)  = ',tBc
        write(23,*) 'resB_c+*(+-)  = ',tBcres
        write(23,*) 'J/PSI(+-)  = ',tJpsi
        write(23,*) 'Upsilon(+-)  = ',tUpsi      


	
	call hrput(0,hstfile,'N')
	
	close(23)
	
	END
	
	






C topor ************random number generator ran2.f

          FUNCTION ran2(idum)
          INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,
     &     IQ2,IR1,IR2,NTAB,NDIV
          REAL ran2,AM,EPS,RNMX 
          PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,
     &     IMM1=IM1-1, 
     &     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, 
     &     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

C Long period (> 2 * 1018) random number generator of 
C L'Ecuyer with Bays-Durham shuffle and 
C added safeguards. Returns a uniform random deviate 
C between 0.0 and 1.0 (exclusive of the endpoint values). 
C Call with idum a negative integer to initialize; thereafter, 
C do not alter idum between successive 
C deviates in a sequence. RNMX should approximate the largest 
C floating value that is less than 1. 

           INTEGER idum2,j,k,iv(NTAB),iy 
           SAVE iv,iy,idum2 
           DATA idum2/123456789/, iv/NTAB*0/, iy/0/ 

           if (idum.le.0) then 
C Initialize.

           idum=max(-idum,1)

C Be sure to prevent idum = 0. 

            idum2=idum 
         do  j=NTAB+8,1,-1 

C Load the shuffle table (after 8 warm-ups).

           k=idum/IQ1


          idum=IA1*(idum-k*IQ1)-k*IR1 
          if (idum.lt.0) idum=idum+IM1 
          if (j.le.NTAB) iv(j)=idum 
          enddo 
C   11   continue
            iy=iv(1)
            endif
           k=idum/IQ1 

C Start here when not initializing. 

          idum=IA1*(idum-k*IQ1)-k*IR1

C Compute idum=mod(IA1*idum,IM1) 
C without overflows by Schrage's method.

            if (idum.lt.0) idum=idum+IM1 
             k=idum2/IQ2 
             idum2=IA2*(idum2-k*IQ2)-k*IR2

C Compute idum2=mod(IA2*idum2,IM2) likewise. 

             if (idum2.lt.0) idum2=idum2+IM2 
              j=1+iy/NDIV

C Will be in the range 1:NTAB. 

              iy=iv(j)-idum2 

C Here idum is shuffled, idum and idum2 are combined to generate output.

              iv(j)=idum 
               if(iy.lt.1)iy=iy+IMM1 
               ran2=min(AM*iy,RNMX) 

C Because users don't expect endpoint values. 

      return 
      END










                         




























