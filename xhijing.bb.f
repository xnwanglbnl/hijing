c**********************************************************************
	PROGRAM XHIJING
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
	character frame*8,targ*8,proj*8,hstfile*16,numfile*16
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
c       
	common/himain1/ natt,eatt,jatt,nt,np,n0,n01,n10,n11
	save  /himain1/
	common/himain2/katt(130000,4),patt(130000,4)
	save  /himain2/
 	common/hijjet1/npj(300),kfpj(300,500),pjpx(300,500),pjpy(300,500)
	1    ,pjpz(300,500),pjpe(300,500),pjpm(300,500)
	2    ,ntj(300),kftj(300,500),pjtx(300,500),pjty(300,500)
	3    ,pjtz(300,500),pjte(300,500),pjtm(300,500)
	save  /hijjet1/
 	common/hijjet2/nsg,njsg(900),iasg(900,3),k1sg(900,100)
	1    ,k2sg(900,100),pxsg(900,100),pysg(900,100),pzsg(900,100)
	2    ,pesg(900,100),pmsg(900,100)
	save  /hijjet2/
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
     &       pb(10,15),ncb(15)
        SAVE  /HIJB/
c       
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
	write(6,*) descr,nseedn
c       Names of the histogram and numbers file. 
	read (5,*) descr,hstfile,numfile
	write(6,*) descr,hstfile,numfile
	write(6,*) ' '

	open(23,file=numfile,status='unknown')

c       set options here
c       These are the parameters set by Miklos Gyulassy
	mg = 1
	if (mg.eq.1) then
c       IHPR2(8):(D=10)The maximum number of jet production per 
c               nucleon-nucleon interaction. When IHPR2(8)=0, jet
c               production will be turned off. When IHPR2(8)<0, the
c               number of jets production will be fixed for each
c               NN collision at a value of abs(IHPR2(8)).
	   ihpr2(8) = 20

c       IHPR2(3):(D=0)switch for one hard scattering with specified PT
c		PT=HIPR1(10) per event.
c		=1:ordinary hard processes;
c		=2:only direct photon production.
	   ihpr2(3) = 0

c       HIPR1(10):(D=-2.25 GeV) To specify the value of PT for one hard
c		scattering generated per event(we call them sample 
c		jets), which could be used to study jet in AA 
c		collision. If HIPR1(10) is negtive, its absolute value
c		gives the low limit of the PT of the sample jets.
	   hipr1(10) = 30.0

c	IHPR2(6):(D=1)switch for the nuclear effect on the parton 
c		distribution function.
	   ihpr2(6) = 0

c       IHPR2(4):(D=1)switch for jet quenching in the the exicited 
c               nuclear matter,i.e. jets only interact with those 
c               nuclear matter that has suffered collisions.
	   ihpr2(4) = 0

c	IHPR2(5):(D=1)switch for the pt kick due to soft initial state 
c		interaction both for valence quarks or diquarks and 
c		gluons.
	   ihpr2(5) = 1

c	IHPR2(2):(D=3)switch for initial and final state radiation in 
c		the hard scattering.
c		=0: both initial and final radiation off;
c		=1: initial radiation on and final radiation off;
c		=2: initial radiation off and final radiation on;
c		=3: both initial and final radiation on.
	   ihpr2(2) = 3

c	HIPR1(13):(D=1.0 fm)The mean free path of a jet to interact 
c		when it goes through an excited nuclear matter.
	   hipr1(13) = 1.

c	HIPR1(14):(D=2.0 GeV/fm) The string tension of a gluon jet. 
c		It is the energy loss dE/dx=HIPR1(14) of the gluon 
cc		jet when the string is streched. For quark jet 
c		dE/dx=0.5*HIPR1(14).
	   hipr1(14) = 1.

c	HIPR1(12):(D=1.0 fm)The transverse distance between a 
c		traversing jet and an excited nucleon(or string system)
c 		below which they will interact and the jet will lose 
c		energy and momentum to the nucleon.
	   hipr1(12) = 1.0

c	HIPR1(29):(D=0.4 fm) the minimum distance between two nucleons
c 		inside a nucleus when the coordinates of all the
c		nucleons in a nucleus are initialized.
	   hipr1(29) = 0.0

c	IHPR2(20):(D=1)Option to turn off(IHPR2(20)=0) the 
c		fragmentation.
	   ihpr2(20) = 1

c	HIPR1(11):(D=2.0 GeV) Minimum pt of a jet which will interact 
c		with excited nuclear matter, or when the pt of a jet 
c		is smaller than HIPR1(11) it will stop interacting 
c		further more.
	   hipr1(11) = 1.0

c	HIPR1(8):(D=2.0 GeV) Minimum pt of hard or semihard scatterings.
	   hipr1(8) = 2.0
	

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
	   parj(1) = 0.02
	   parj(2) = 0.23
	   parj(3) = 0.70
	   parj(4) = 0.05
	   parj(11) = 0.5

c
c       irope, rkappa:
c	The variable rope indicates a larger string tension.  The variable
c	kappa indicates the value of the new tension divided by the value
c	of the old tension.
c
	   irope = 0
	   rkappa = 1.4
	   if (irope.eq.1) then
	      parj(1) = parj(1)**(1.0/rkappa)
	      parj(2) = parj(2)**(1.0/rkappa)
	      parj(3) = parj(3)**(1.0/rkappa)
	      parj(4) = 1.0/3.0*((3.0*parj(4))**(1.0/rkappa))
	      parj(21) = parj(21)*sqrt(rkappa)
	   end if 

c
c       ib(1):  0,1 for switching off/on the baryon stopping interactions
c
	   ib(1) = 1
c       
c       ib(2):  0,1 for switching off/on the 1 baryon stopping interaction
c       rb(2):  cross section (mb) of the 1 baryon stopping component
c       
	   ib(2) = 1
	   rb(2) = 17.0
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
	   ib(6) = 1
	   rb(6) = 6.0
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
c
	nbje = 0
c       set parameters for the baryon junctions. 
	if (nbje.eq.1) then
c
c       string fragmentation parameters, see manual for Jetset v.7.2. 
c       parj(1):(D=0.10),parj(2):(D=0.30),parj(3):(D=0.40),parj(4):(D=0.05),
c       parj(11):(D=0.5)
c
	   parj(1) = 0.05
	   parj(2) = 0.23
	   parj(3) = 0.70
	   parj(4) = 0.05
	   parj(11) = 0.5

c
c       irope, rkappa:
c	The variable rope indicates a larger string tension.  The variable
c	kappa indicates the value of the new tension divided by the value
c	of the old tension.
c
	   irope = 0
	   rkappa = 1.45
	   if (irope.eq.1) then
	      parj(1) = parj(1)**(1.0/rkappa)
	      parj(2) = parj(2)**(1.0/rkappa)
	      parj(3) = parj(3)**(1.0/rkappa)
	      parj(4) = 1.0/3.0*((3.0*parj(4))**(1.0/rkappa))
	      parj(21) = parj(21)*sqrt(rkappa)
	   end if 

c
c       ib(1):  0,1 for switching off/on the baryon stopping interactions
c
	   ib(1) = 1
c       
c       ib(2):  0,1 for switching off/on the 1 baryon stopping interaction
c       rb(2):  cross section (mb) of the 1 baryon stopping component
c       
	   ib(2) = 1
	   rb(2) = 18.0
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
	   ib(6) = 0
	   rb(6) = 0.0
c       
c       rb(7): minimum mass of the baryon - anti-baryon production. 
c       
	   rb(7) = 0.0
c
c	ib(11):  baryon stopping with diquark breaking
c	rb(11):  cross section (mb) for the diquark breaking 
c	
	ib(11) = 0
	rb(11) = 0.0
c	
	end if


 
c       The cut on eta for the pt distributions. 
	eta_cut = 1.0

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
c
	ippad = ( 0.1*ihnt2(1))
	itpad = ( 0.1*ihnt2(3)) 
	if(ihnt2(1).gt.1) then 
	   call hbook1(50,'projectile participants',ihnt2(1)+ippad,0.0,
	1	real(ihnt2(1)+ippad),0.0)
	   call hidopt(50,'logy')
	   ipt=ihnt2(1)+ihnt2(3)+ippad+itpad
	   call hbook1(52,'total participants',ipt,0.0,real(ipt),0.0)
	   call hidopt(52,'logy')
	   call hbook1(53,'proj eveto/ebeam',ihnt2(1)+ippad,0.0,
	1	real(ihnt2(1)+ippad),0.0)
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
	1	real(ihnt2(3)+itpad),0.0)
	   call hidopt(51,'logy')
	   call hbook1(54,'targ eveto/etarg',ihnt2(3)+itpad,0.0,
	1	real(ihnt2(3)+itpad),0.0)
	   call hidopt(54,'logy')
	   xxmax = 2.4*ihnt2(3)**0.3333
	   delrr = 0.2
	   nnmax = xxmax/delrr
	   xxmax = nnmax*delrr
	   call hbook1(56,'rho(r), targ',nnmax,0.0,xxmax,0.0)
	   call hbook1(59,'dn(r),  targ',nnmax,0.0,xxmax,0.0)
	   call hbook1(57,'number of binary',100,0.0,
	1	1.5*ihnt2(3)**1.3333,0.0)
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

c       pt distributions for the photons.
 	call hbookb(115,'dN/dpt**2, photon(res decay not pi0)',
	1	      n_ptmx,pthbk,0.)
	call hidopt(115,'logy')
	call hbookb(116,'dN/dpt**2, direct photon',n_ptmx,pthbk,0.)
	call hidopt(116,'logy')
	call hbookb(117,'dN/dpt**2, pi0',n_ptmx,pthbk,0.)
	call hidopt(117,'logy')

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
        call hbookb(140,'dN/dpt**2, D(+-)',n_ptmx,pthbk,0.)
	call hidopt(140,'logy')

c       pt distributions for the baryons. 
	call hbookb(150,'dN/dpt**2, p',n_ptmx,pthbk,0.)
	call hidopt(150,'logy')
	call hbookb(151,'dN/dpt**2, p~',n_ptmx,pthbk,0.)
	call hidopt(151,'logy')
c
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

c       Create the rapidity distributions of the q's, g's and photons. 
	call hbook1(201,'dN/dy, charged',nymx,ymin,ymax,0.0)

	CALL HBOOK1(205,'dN/dy,  gluons',NYMX,YMIN,YMAX,0.0)
	CALL HBOOK1(206,'dN/dy,  q+qbar',NYMX,YMIN,YMAX,0.0)
	CALL HBOOK1(207,'dN/dy,  qq+qqbar',NYMX,YMIN,YMAX,0.0)
	CALL HBOOK1(210,'dN/dy,  c+cbar',NYMX,YMIN,YMAX,0.0)

	call hbook1(215,'dN/dy, photons',nymx,ymin,ymax,0.0)
	call hbook1(217,'dN/dy, pi0',nymx,ymin,ymax,0.0)	
	
c       Create the rapidity distributions of the mesons.  

	call hbook1(220,'dN/dy, pi+',nymx,ymin,ymax,0.0)
	call hbook1(221,'dN/dy, pi-',nymx,ymin,ymax,0.0)
	call hbook1(222,'dN/dy, K+',nymx,ymin,ymax,0.0)
	call hbook1(223,'dN/dy, K-',nymx,ymin,ymax,0.0)
	call hbook1(224,'dN/dy, K0',nymx,ymin,ymax,0.0)
	call hbook1(225,'dN/dy, K~0',nymx,ymin,ymax,0.0)
	call hbook1(226,'dN/dy, K_S0',nymx,ymin,ymax,0.0)
c	call hbook1(227,'dN/dy, K_L0',nymx,ymin,ymax,0.0)
	call hbook1(228,'dN/dy, K*0',nymx,ymin,ymax,0.0)
	call hbook1(240,'dN/dy, D(+,-)',nymx,ymin,ymax,0.0)

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
	call hbook1(315,'dN/deta,  photons',nymx,ymin,ymax,0.0)
c       
c       Create the pseudo-rapidity distributions of the mesons. 
	call hbook1(320,'dN/deta, pi+',nymx,ymin,ymax,0.0)
	call hbook1(321,'dN/deta, pi-',nymx,ymin,ymax,0.0)
	call hbook1(322,'dN/deta, K+',nymx,ymin,ymax,0.0)
	call hbook1(323,'dN/deta, K-',nymx,ymin,ymax,0.0)
	call hbook1(324,'dN/deta, K0',nymx,ymin,ymax,0.0)
	call hbook1(325,'dN/deta, K~0',nymx,ymin,ymax,0.0)
	call hbook1(326,'dN/deta, K_S0',nymx,ymin,ymax,0.0)
c	call hbook1(327,'dN/deta, K_L0',nymx,ymin,ymax,0.0)
	call hbook1(328,'dN/deta, K*0',nymx,ymin,ymax,0.0)
	call hbook1(340,'dN/deta, D(+,-)',nymx,ymin,ymax,0.0)

c       Create the pseudo-rapidity distributions of the baryons.  
	call hbook1(350,'dN/deta, p',nymx,ymin,ymax,0.0)
	call hbook1(351,'dN/deta, p~',nymx,ymin,ymax,0.0)
	call hbook1(352,'dN/deta, n',nymx,ymin,ymax,0.0)
	call hbook1(353,'dN/deta, n~',nymx,ymin,ymax,0.0)
	call hbook1(354,'dN/deta, Lambda',nymx,ymin,ymax,0.0)
	call hbook1(355,'dN/deta, Lambda~',nymx,ymin,ymax,0.0)
c
c       Create the det/dy distributions
	call hbook1(400,'det/dy',nymx,ymin,ymax,0.0)
	call hbook1(401,'det/dy char',nymx,ymin,ymax,0.0)
	call hbook1(402,'det/dy neut',nymx,ymin,ymax,0.0)

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
	tneut = 0.0
	tchrm = 0.0
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
	      ttot = ttot + w_event
	      if (luchge(katt(i,1)).eq.0) then
		 call hfill(402,rap,0.0,et*w_event/w_y)
		 tneut = tneut + w_event
	      else
		 if(abs(eta).le.eta_cut) then
		    call hfill(101,ptr,0.0,w_event/2.0/max(ptr,0.00001))
		    nch_cor=nch_cor+1
		    pt_cor=pt_cor+ptr
		 endif		 
		 call hfill(201,eta,0.0,w_event/w_y)
		 call hfill(301,rap,0.0,w_event/w_y)
		 call hfill(401,rap,0.0,et*w_event/w_y)
		 tchar = tchar + w_event
		 if (luchge(katt(i,1)).lt.0) tcharm = tcharm + w_event
	      end if	
	      ptot=ptot+ptr
c       
c       
c       record the gammas:
c 
	      if (katt(i,1).eq.22) then
c       gammas from all resonances except the pi0 decay !!
		 call hfill(115,eta,0.0,w_event/w_y)
		 call hfill(315,rap,0.0,w_event/w_y)
		 if (katt(i,2).ne.40) then       
		    if (abs(eta).le.eta_cut) call hfill(216,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		 else
c       gammas from direct QCD proc
		    if(abs(eta).le.eta_cut) call hfill(117,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		 endif
c	
c       record the mesons:
c       
c       pi0
	      else if (abs(katt(i,1)).lt.1000) then
		 if (katt(i,1).eq.111) then
		    call hfill(317,eta,0.0,w_event/w_y)
		    if (abs(eta).le.eta_cut) call hfill(117,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		    tpi0 = tpi0 + w_event
c       pi+
		 else if (katt(i,1).eq.211) then
		    call hfill(320,eta,0.0,w_event/w_y)
		    call hfill(220,rap,0.0,w_event/w_y)
		    if(abs(eta).le.eta_cut) call hfill(120,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		    tpip = tpip + w_event
c       pi-
		 else if(katt(i,1).eq.-211) THEN
		    call hfill(321,eta,0.0,w_event/w_y)
		    call hfill(221,rap,0.0,w_event/w_y)
		    if(abs(eta).le.eta_cut) call hfill(121,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))	
		    tpim = tpim + w_event
		    
c       K0
		 else if(katt(i,1).eq.311) then
		    call hfill(224,rap,0.0,w_event/w_y)
		    tK0 = tK0 + w_event
c       K~0
		 else if(katt(i,1).eq.-311) then
		    call hfill(326,rap,0.0,w_event/w_y)
		    taK0 = taK0 + w_event
c       K+
		 else if(katt(i,1).eq.321) then
		    call hfill(322,eta,0.0,w_event/w_y)
		    call hfill(222,rap,0.0,w_event/w_y)
		    if (abs(eta).le.eta_cut) call hfill(122,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		    tKp = tKp + w_event
c       K-
		 else if(katt(i,1).eq.-321) then
		    call hfill(323,eta,0.0,w_event/w_y)
		    call hfill(223,rap,0.0,w_event/w_y)
		    if (abs(eta).le.eta_cut) call hfill(123,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		    tKm = tKm + w_event
		    
c       K_S0
		 else if(katt(i,1).eq.310) then
		    call hfill(226,rap,0.0,w_event/w_y)
		    tKS0 = tKS0 + w_event
c       K*0
		 else if(katt(i,1).eq.313) then
		    call hfill(228,rap,0.0,w_event/w_y)
		    tKstar = tKstar + w_event
c       D+
		 else if(abs(katt(i,1)).eq.411) then
		    call hfill(340,eta,0.0,w_event/w_y)
		    call hfill(240,rap,0.0,w_event/w_y)
		    if(abs(eta).le.eta_cut) call hfill(140,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))		 
		    tD = tD + w_event
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
		    if(abs(eta).le.eta_cut) call hfill(150,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		    tp = tp + w_event
c       p~
		 else if(katt(i,1).eq.-2212) then
		    call hfill(351,eta,0.0,w_event/w_y)
		    call hfill(251,rap,0.0,w_event/w_y)
		    if(abs(eta).le.eta_cut) call hfill(151,ptr,0.0,
	1		 w_event/2.0/max(ptr,0.00001))
		    tap = tap + w_event
c       n
		 else if (katt(i,1).eq.2112) then
		    call hfill(252,rap,0.0,w_event/w_y)
		    tn = tn + w_event
c       n~
		 else if (katt(i,1).eq.-2112) then
		    call hfill(253,rap,0.0,w_event/w_y)
		    tan = tan + w_event
c       Lambda
		 else if (katt(i,1).eq.3122) then
		    call hfill(254,rap,0.0,w_event/w_y)
		    tlam = tlam + w_event
c       Lambda~
		 else if (katt(i,1).eq.-3122) then
		    call hfill(255,rap,0.0,w_event/w_y)
		    talam = talam + w_event
		    
c       Sigma-
		 else if (katt(i,1).eq.3112) then
		    call hfill(256,rap,0.0,w_event/w_y)
		    tsigm = tsigm + w_event
c       Sigma~+
		 else if (katt(i,1).eq.-3112) then
		    call hfill(257,rap,0.0,w_event/w_y)
		    tasigp = tasigp + w_event
c       Sigma0
		 else if (katt(i,1).eq.3212) then
		    call hfill(258,rap,0.0,w_event/w_y)
		    tsig0 = tsig0 + w_event
c       Sigma~0
		 else if (katt(i,1).eq.-3212) then
		    call hfill(259,rap,0.0,w_event/w_y)
		    tasig0 = tasig0 + w_event
c       Sigma+
		 else if (katt(i,1).eq.3222) then
		    call hfill(260,rap,0.0,w_event/w_y)
		    tsigp = tsigp + w_event
c       Sigma~-
		 else if (katt(i,1).eq.-3222) then
		    call hfill(261,rap,0.0,w_event/w_y)
		    tasigm = tasigm + w_event
c       Xi-
		 else if (katt(i,1).eq.3312) then
		    call hfill(262,rap,0.0,w_event/w_y)
		    txim = txim + w_event
c       Xi~+
		 else if (katt(i,1).eq.-3312) then
		    call hfill(263,rap,0.0,w_event/w_y)
		    taxip = taxip + w_event
c       Xi0
		 else if (katt(i,1).eq.3322) then
		    call hfill(264,rap,0.0,w_event/w_y)
		    txi0 = txi0 + w_event
c       Xi~0
		 else if (katt(i,1).eq.-3322) then
		    call hfill(265,rap,0.0,w_event/w_y)
		    taxi0 = taxi0 + w_event
c       Omega-
		 else if (katt(i,1).eq.3334) then
		    call hfill(266,rap,0.0,w_event/w_y)
		    tomem = tomem + w_event
c       Omega~+
		 else if (katt(i,1).eq.-3334) then
		    call hfill(267,rap,0.0,w_event/w_y)
		    taomep = taomep + w_event		    
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
	      if(jatt.eq.0) call hfill(161,xn,0.0,w_event/w_an)
	      if(jatt.eq.1) call hfill(162,xn,0.0,w_event/w_an)
	      if(jatt.gt.1) call hfill(163,xn,0.0,w_event/w_an)
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
		    call hfill(205,rap,0.0,w_event/w_y)
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
c       include the bin width of the pT histograms
	id(1) = 101
	id(2) = 117
	id(3) = 120
	id(4) = 121
	id(5) = 122
	id(6) = 123
	id(7) = 140
	id(8) = 150
	id(9) = 151
	jjmax = 9

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
	   
	
C*******************************************
 4001	FORMAT(1X,I10,3X,'Ntot=',I10,3X,'Nch=',I10,
     &	3X,'Njet=',I10/1X,'***************************'/)
	
	PTOT=PTOT*W_EVENT/ANTOT
	PT_CUT=PT_CUT_TOT*W_EVENT/AN_CUT
	
	R_RMS=SQRT(R_RMS/FLOAT(IAP+IAT)*W_EVENT)
	write(6,*) ' geometry rms radius av =',r_rms
	write(6,*) ' DETOT =',DETOT*w_event,' ANTOT =',ANTOT,' <PT> =',PTOT
	write(6,*) 'AN_CUT=', AN_CUT,'  PT_CUT=',PT_CUT
	write(6,*) ' E_CM=',HINT1(1),' I_DIPOLE=',IHPR2(1),' I_RAD=',IHPR2(2)
	write(6,*) 'Q_CUT=',HIPR1(5),' PTRAD_CUT=',HIDAT(2),' P_RAD=',HIDAT(3)
	write(6,*) 'MCUT=',HIPR1(1),' PT_INIRAD',HIDAT(1),' CRS_SD=',HIDAT(4)
	write(6,*) 'PTRMS=',HIPR1(2),' a=',PARJ(41),' b=',PARJ(42)
	write(6,*) ' J_JMAX=',IHPR2(8),' HINT1(11)',HINT1(11),'Stot',HINT1(13)
	write(6,*) 'HINT1(10) =',HINT1(10),'Sinel =',HINT1(12)
	write(6,*) 'S_Trig =',HINT1(60),' S_trig_t =',HINT1(59)
	
	write(6,*) ' '
	write(6,*) ' NHIT= ',NHIT,' OUT OF NEVENT =', nhit+nmiss
	write(6,*) ' '
	
c       reaction cross section in mb
	sigreac=3.1415*(bmax*bmax-bmin*bmin)*10.*float(nhit)/real(n_event)
	write(6,*) ' reaction cross section mb = ',sigreac
	write(6,*) ' '

	write(23,*) 'Total numbers of the various particles: '
	write(23,*) ' '
	write(23,*) 'total   =  ',ttot	
	write(23,*) 'charged =  ',tchar
	write(23,*) 'neutral =  ',tneut
	write(23,*) ' ' 
	write(23,*) 'h-      =  ',tcharm
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
	
	call hrput(0,hstfile,'N')
	
	close(23)
	
	END
	
	













