c     The subroutine bab models baryon-antibaryon creation according to
c     the Regge exchange of two Pomerons and a M^J_0 Reggeon.
c
c     Ref:  S.E. Vance and M. Gyulassy, Phys. Rev. Lett., submitted. 
c
c     by S.E. Vance, Sept 1, 1998.
c     
      SUBROUTINE HIJBAB(JTP,NTP,IERROR)

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/     
      COMMON/HISTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
      SAVE  /HISTRNG/
      COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &     PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &     PJPM(300,500),NTJ(300),KFTJ(300,500),
     &     PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &     PJTE(300,500),PJTM(300,500)
      SAVE  /HIJJET1/
      COMMON/HIJB/nb,ib(15),rb(15),sb(5),kb(10,900),kfb(10,15),
     &     pb(10,15),ncb(15)
      SAVE  /HIJB/
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)   
      SAVE  /LUJETS/
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE  /LUDAT1/
c     
      EXTERNAL ULMASS, RAPMAX, KFBARYON, LUKFDI
      
      DIMENSION e(5),px(5),py(5),pz(5),pm(5),y(5),yvq(5),ysq(5)
      DIMENSION kfvq(5),kfsq(5),emvq(5),epvq(5),emsq(5),epsq(5)
      DIMENSION pxvq(5),pyvq(5),pxsq(5),pysq(5),vqmt(5),sqmt(5)
      DIMENSION vqm(5),sqm(5),ihrd(5),kmes(5),dyms(5)
      DIMENSION kfbq(5),kfabq(5),embq(5),epbq(5),emabq(5),epabq(5)
      DIMENSION pxab(5),pyab(5),pxb(5),pyb(5),bqmt(5),abqmt(5)
      DIMENSION bqm(5),abqm(5),kfl(5)

      IERROR = 0
      nbhit = 0
      nb = 0

c     set parameters
      ptwid = sb(3)

c     determine if the string has sufficient energy to decay.
      if (ntp.eq.1) then
         smass = sqrt(max((pp(jtp,4))**2 - (pp(jtp,3))**2 
     1        - pp(jtp,2)**2 - pp(jtp,1)**2,0.0001))
      else
         smass = sqrt(max((pt(jtp,4))**2 - (pt(jtp,3))**2 
     1        - pt(jtp,2)**2 - pt(jtp,1)**2,0.0001))
      end if
      if (smass.lt.rb(7)) then
         ierror = 1
         nb = 1
         return
      end if 

c     Determine the rapidities of the quark and the diquark.
c     pmtzp:  mT of the +z moving particle 
c     pmtzn:  mT of the -z moving particle
c     
c     The constituent quark masses are used and are flagged by setting
c     mstj(93) = 1.
c
      if (ntp.eq.1) then
         kfq = min(iabs(nfp(jtp,1)),iabs(nfp(jtp,2)))
         kfqq = max(iabs(nfp(jtp,1)),iabs(nfp(jtp,2)))
         kfsign = nfp(jtp,1)/iabs(nfp(ntp,1))
         ep = pp(jtp,4) + pp(jtp,3)
         em = pp(jtp,4) - pp(jtp,3)
         pmtzp = sqrt(pp(jtp,15)**2 + pp(jtp,8)**2 + pp(jtp,9)**2)
         mstj(93) = 1
         pmtzn = sqrt(ulmass(kfq)**2 + pp(jtp,6)**2 + pp(jtp,7)**2)
         pxvq(1) = pp(jtp,8)
         pyvq(1) = pp(jtp,9)
         pxvq(2) = pp(jtp,6)
         pyvq(2) = pp(jtp,7)
         kfvq(1) = kfqq
         kfvq(2) = kfq
         vqm(1) = ulmass(kfqq)
         vqmt(1) = sqrt(vqm(1)**2 + pxvq(1)**2 + pyvq(1)**2)
         mstj(93) = 1
         vqm(2) = ulmass(kfq)
         vqmt(2) = sqrt(vqm(2)**2 + pxvq(2)**2 + pyvq(2)**2)

      else if (ntp.eq.2) then
         kfq = min(iabs(nft(jtp,1)),iabs(nft(jtp,2)))
         kfqq = max(iabs(nft(jtp,1)),iabs(nft(jtp,2)))
         kfsign = nft(jtp,1)/iabs(nft(jtp,1))
         ep = pt(jtp,4) + pt(jtp,3)
         em = pt(jtp,4) - pt(jtp,3)
         pmtzn = sqrt(pt(jtp,15)**2 + pt(jtp,8)**2 + pt(jtp,9)**2)
         mstj(93) = 1
         pmtzp = sqrt(ulmass(kfq)**2 + pt(jtp,6)**2 + pt(jtp,7)**2)
         pxvq(1) = pt(jtp,6)
         pyvq(1) = pt(jtp,7)
         pxvq(2) = pt(jtp,8)
         pyvq(2) = pt(jtp,9)
         kfvq(1) = kfq
         kfvq(2) = kfqq
         vqm(2) = ulmass(kfqq)
         vqmt(2) = sqrt(vqm(2)**2 + pxvq(2)**2 + pyvq(2)**2)
         mstj(93) = 1
         vqm(1) = ulmass(kfq)
         vqmt(1) = sqrt(vqm(1)**2 + pxvq(1)**2 + pyvq(1)**2)

      else 
         write(6,*) 'error in rapidity identification'
      end if

c   
c     solve for the rapidities of the quark and the diquark.
c     iend choses which particles rapidity is solved for the +z or -z moving
c     particle
      if (ntp.eq.1) then
         iend = 1
        else 
         iend = 0
      end if 
      CALL RAPMAX(iend,ep,em,pmtzp,pmtzn,ymax,ymin,irperr)
      if (ib(5).eq.1) write(6,*) 'ymin,ymax',ymin,ymax
      if (irperr.eq.1) then
         nbhit = 200
      end if 
       
      if (ymax.lt.ymin) then
         nbhit = 200
         write(6,*) 'BAB Failed, initial string is too short'
      end if

c     Determine the flavors of the three sea quarks and the baryon. 
      CALL KFBARYON(3,kfbq1,kfbq2,kfbq3,kbar)
      kfbq(1) = kfbq1
      kfbq(2) = kfbq2
      kfbq(3) = kfbq3
      kfbar   = kbar

      if (ntp.eq.1) then
         kfsq(2) = -kfbq(3)
      else
         kfsq(1) = -kfbq(3)
      end if 
      
c     Determine the flavors of the three anti-sea quarks and the anti-baryon.
      CALL KFBARYON(3,kfabq1,kfabq2,kfabq3,kabar)
      kfabq(1) = -kfabq1
      kfabq(2) = -kfabq2
      kfabq(3) = -kfabq3
      kfabar   = -kabar

      if (ntp.eq.1) then
         kfsq(1) = -kfabq(3)
      else
         kfsq(2) = -kfabq(3)
      end if 

c
 10   nbhit = nbhit + 1
      if (nbhit.gt.200) then
         ierror = 2
         nb = 1
c         write(6,*) 'kfbar,kfabar',kfbar,kfabar
         return
      endif 
c
c     set counting parameters 
      nstr = 0

c     statement so that for a given pT, several rapidities are chosen.
      if (mod(nbhit,5).ne.1) goto 35 

c     Determine the pT and the mT of the sea quarks, the valence quarks
c     and the baryon according to a gaussian pT distribution.      
      do 30, i = 1,3
         ptmag = ptwid*sqrt(-log(max(1E-10,(1-rlu(0)))))
         ptphi = 2.0*hipr1(40)*rlu(0)
         pxb(i) = ptmag*cos(ptphi)
         pyb(i) = ptmag*sin(ptphi)
 30   continue
      if (ntp.eq.1) then         
         pxsq(1) = -pxb(3)
         pysq(1) = -pyb(3)
      else
         pxsq(2) = -pxb(3)
         pysq(2) = -pyb(3)
      end if

       do 33, i = 1,3
            ptmag = ptwid*sqrt(-log(max(1E-10,(1-rlu(0)))))
            ptphi = 2.0*hipr1(40)*rlu(0)
            pxab(i) = ptmag*cos(ptphi)
            pyab(i) = ptmag*sin(ptphi)
 33      continue

      if (ntp.eq.1) then
         pxsq(2) = -pxab(3)
         pysq(2) = -pyab(3)
      else 
         pxsq(1) = -pxab(3)
         pysq(1) = -pyab(3)
      end if 

 35   continue

c     Determine the rapidity of the baryon and the antibaryon. 
      dystr = abs(ymax - ymin)
      yhalf = ymin + 0.5*dystr

c     The 2.0 in the statement below takes into account the rapidity 
c     intervals of the two outer strings. 
      ybab = yhalf + (rlu(0)-0.5)*(dystr - 2.0)

      ybabmin = 0.0
c
c     The 0.8 below takes into account the rapidity interval for 
c     the q-anti-q end string.  
      if (ntp.eq.1) then
         ybabmax = ybab - ymin - 0.8 
      else 
         ybabmax = ymax - ybab - 0.8
      end if 

      if (ybabmin.gt.ybabmax) goto 10

      abab = -0.5
      pnorm = exp(abab*ybabmax)-exp(abab*ybabmin)
      dy    = 1/abab*log(exp(abab*ybabmin) + pnorm*rlu(0))
      
      if (ntp.eq.1) then
         yb = ybab - dy
      else
         yb = ybab + dy
      end if 
      yab = ybab

c     allow for the difference in mass between a s and a u/d quark. 
      if (kfabq(1).eq.3.or.kfabq(2).eq.3) then
         dymax = 1.3
      else
         dymax = 0.85
      end if 

c     Allow mesons to be produced between the baryon - anti-baryon pair if 
c     possible. 
      if (dy.lt.dymax) then
         nstr = 0
         kfabq(1) = -kfbq(1)
         kfabq(2) = -kfbq(2)
         kfabq1 = -kfabq(1)
         kfabq2 = -kfabq(2)
         kfabq3 = -kfabq(3)
         CALL KFBARYON(0,kfabq1,kfabq2,kfabq3,kabar)
         kfab = -kabar

         pxab(1) = -pxb(1)
         pyab(1) = -pyb(1)
         pxab(2) = -pxb(2)
         pyab(2) = -pyb(2)

         epstr = 0.0
         emstr = 0.0

      else
         if (yb.gt.yab) then
            kfsq(3) = -kfbq(1)
            kfsq(4) = -kfbq(2)
            kfvq(3) = -kfabq(1)
            kfvq(4) = -kfabq(2)
            ysq(3) = yb
            yvq(3) = yab
            ysq(4) = yb
            yvq(4) = yab            

            pxsq(3) = -pxb(1)
            pysq(3) = -pyb(1)
            pxsq(4) = -pxb(2)
            pysq(4) = -pyb(2)
            pxvq(3) = -pxab(1)
            pyvq(3) = -pyab(1)
            pxvq(4) = -pxab(2)
            pyvq(4) = -pyab(2)

         else
            kfsq(3) = -kfabq(1)
            kfsq(4) = -kfabq(2)
            kfvq(3) = -kfbq(1)
            kfvq(4) = -kfbq(2)
            ysq(3) = yab
            yvq(3) = yb
            ysq(4) = yab
            yvq(4) = yb 

            pxvq(3) = -pxb(1)
            pyvq(3) = -pyb(1)
            pxvq(4) = -pxb(2)
            pyvq(4) = -pyb(2)
            pxsq(3) = -pxab(1)
            pysq(3) = -pyab(1)
            pxsq(4) = -pxab(2)
            pysq(4) = -pyab(2)

         end if 
         nstr = 2

         do 37, i = 3,4         
            mstj(93) = 1
            sqm(i)   = ulmass(kfsq(i))
            sqmt(i)  = sqrt(sqm(i)**2 + pxsq(i)**2 + pysq(i)**2)
            epsq(i)  = sqmt(i)*exp(ysq(i))
            emsq(i)  = sqmt(i)*exp(-ysq(i))
            mstj(93) = 1
            vqm(i)   = ulmass(kfvq(i))
            vqmt(i)  = sqrt(vqm(i)**2 + pxvq(i)**2 + pyvq(i)**2)
            epvq(i)  = vqmt(i)*exp(yvq(i))
            emvq(i)  = vqmt(i)*exp(-yvq(i))
 37      continue
      
         epstr = epsq(3) + epvq(3) + epsq(4) + epvq(4)
         emstr = emsq(3) + emvq(3) + emsq(4) + emvq(4)
      end if

c
      ptmag = 0.6*sqrt(-log(max(1E-10,(1-rlu(0))))) 
      ptphi = 2.0*hipr1(40)*rlu(0)
      pxkick = ptmag*cos(ptphi)
      pykick = ptmag*sin(ptphi)
c
c     Determine the kinematic properties of the baryon.
      pxbar  = pxb(1) + pxb(2) + pxb(3) + pxkick 
      pybar  = pyb(1) + pyb(2) + pyb(3) + pykick
      ptbar  = sqrt(pxbar**2 + pybar**2)
      mstj(24) = 0
      bmass  = ulmass(kfbar)     
      mstj(24) = 2
      barmt  = sqrt(bmass**2 + ptbar**2)
      pzb  = barmt*sinh(yb)
      eb   = barmt*cosh(yb) 
      epb  = eb + pzb
      emb  = eb - pzb

      if ((eb**2-pxbar**2-pybar**2-pzb**2).le.0) goto 10

c     Determine the kinematic properties of the anti-baryon.
      pxabar = pxab(1) + pxab(2) + pxab(3) - pxkick
      pyabar = pyab(1) + pyab(2) + pyab(3) - pykick
      mstj(24) = 0 
      abmass = ulmass(kfabar)
      mstj(24) = 2
      ptabar = sqrt(pxabar**2 + pyabar**2)
      abarmt = sqrt(abmass**2 + ptabar**2)
      pzb  = barmt*sinh(yb)
      eb   = barmt*cosh(yb) 
      epb  = eb + pzb
      emb  = eb - pzb
      pzab  = abarmt*sinh(yab)
      eab   = abarmt*cosh(yab) 
      epab  = eab + pzab
      emab  = eab - pzab
      
      if ((eab**2-pxabar**2-pyabar**2-pzab**2).le.0) goto 10
      
      if (ntp.eq.1) then
         ysq(1) = yab
         ysq(2) = yb
      else
         ysq(1) = yb
         ysq(2) = yab
      end if
         
      do 40, i = 1,2         
         mstj(93) = 1
         sqm(i)   = ulmass(kfsq(i))
         sqmt(i)  = sqrt(sqm(i)**2 + pxsq(i)**2 + pysq(i)**2)
         epsq(i)  = sqmt(i)*exp(ysq(i))
         emsq(i)  = sqmt(i)*exp(-ysq(i))
 40   continue

c     Determine the E+ and E- available for the creation of end strings.
      epr = ep - epb - epab - epsq(1) - epsq(2) - epstr
      emr = em - emb - emab - emsq(1) - emsq(2) - emstr
      if (epr.lt.0.0.or.emr.lt.0.0) then
c     if (ib(5).eq.1) write(6,*) 'ep,em,',ep,em,epb,epab,emb,emab
         if (ib(5).eq.1) write(6,*) 'epr,emr',epr,emr
         goto 10
      end if 
      
      if (ntp.eq.1) then
         CALL RAPMAX(0,epr,emr,vqmt(1),vqmt(2),yrmax,yrmin,irperr)
         if (irperr.eq.1) goto 10
         
         epvq(1) = pmtzp*exp(yrmax)
         emvq(1) = pmtzn*exp(-yrmax)
            
         epeff = ep - epb - epab - epsq(1) - epvq(1) - epsq(2) - epstr
         emeff = em - emb - emab - emsq(1) - emvq(1) - emsq(2) - emstr

         if (ib(5).eq.1) write(6,*) 'epeff,emeff',epeff,emeff

         if (epeff.le.0.0.or.emeff.le.0.0) goto 10
         
         epvq(2) = epeff
         emvq(2) = emeff
         
         yvq(2) = 0.5*log(epvq(2)/emvq(2))

      else
         CALL RAPMAX(1,epr,emr,vqmt(1),vqmt(2),yrmax,yrmin,irperr)
         if (irperr.eq.1) goto 10
         
         epvq(2) = pmtzp*exp(yrmin)
         emvq(2) = pmtzn*exp(-yrmin)
            
         epeff = ep - epb - epab - epsq(1) - epsq(2) - epvq(2) - epstr
         emeff = em - emb - emab - emsq(1) - emsq(2) - emvq(2) - emstr
         
c         if (ib(5).eq.1) write(6,*) 'epeff,emeff',epeff,emeff

         if (epeff.le.0.0.or.emeff.le.0.0) goto 10
         
         epvq(1) = epeff
         emvq(1) = emeff
         
         yvq(1) = 0.5*log(epvq(1)/emvq(1))  

      end if 

c     Determine the E, Pz, Px, Py, y and m of the strings and baryon.
      do 100, i = 1,2+nstr
         e(i) = ((epsq(i)+emsq(i))+(epvq(i)+emvq(i)))/2.0
         pz(i) = ((epsq(i)-emsq(i))+(epvq(i)-emvq(i)))/2.0
         if (abs(pz(i)).gt.e(i)) goto 10
         if (((epsq(i)+epvq(i))/(emsq(i)+emvq(i))).lt.0.0) goto 10
         y(i) = 1./2.*log((epsq(i)+epvq(i))/(emsq(i)
     1        +emvq(i)))
         px(i) = pxsq(i) + pxvq(i)
         py(i) = pysq(i) + pyvq(i) 
         if ((e(i)**2-px(i)**2-py(i)**2-pz(i)**2).le.0.0) goto 10
         pm(i) = sqrt(e(i)**2 - px(i)**2 - py(i)**2 - pz(i)**2)
c    allow for the difference in mass between a s and a u/d quark.
c     for the quark - diquark end string
         if (ntp.eq.1.and.i.eq.1.and.pm(i).lt.1.0+sqm(i)) goto 10
         if (ntp.eq.2.and.i.eq.2.and.pm(i).lt.1.0+sqm(i)) goto 10
c     for the quark - anti-quark string
         if (ntp.eq.1.and.i.eq.2.
     1        and.pm(i).lt.sqm(i)+vqm(i)-0.225) goto 10
         if (ntp.eq.2.and.i.eq.1.
     1        and.pm(i).lt.sqm(i)+vqm(i)-0.225) goto 10 
 100  continue
c
c     require that the string fragment correctly
c     
      do 110, i = 1,2 + nstr
         amt = pm(i)**2 + px(i)**2 + py(i)**2
         amt1 = vqmt(i)**2
         amt2 = sqmt(i)**2
         pzcm2 = (AMT**2+AMT1**2+AMT2**2-2.0*AMT*AMT1
     &        -2.0*AMT*AMT2-2.0*AMT1*AMT2)/4.0/AMT
         if (pzcm2.le.0) goto 10
 110  continue

      if (nstr.eq.0) then
         pxstr = 0.0
         pystr = 0.0
         pzstr = 0.0
         estr = 0.0
      else if (nstr.eq.2) then
         pxstr = px(3) + px(4)
         pystr = py(3) + py(4)
         pzstr = pz(3) + pz(4)
         estr = e(3) + e(4)
      end if

c     Check for energy conservation
      pxtot = px(1) + px(2) + pxbar + pxabar + pxstr 
      pytot = py(1) + py(2) + pybar + pyabar + pystr
      pztot = pz(1) + pz(2) + pzb + pzab + pzstr
      etot  = e(1)  + e(2)  + eab + eb + estr
      
      if (ntp.eq.1) then
         if (((etot-pp(jtp,4)).gt.0.0001).or.((pztot-pp(jtp,3))
     1        .gt.0.0001).or.((pytot-pp(jtp,2)).gt.0.0001)
     2        .or.((pxtot-pp(jtp,1)).gt.0.0001)) then
            write(6,*) 'energy not conserved'
         end if 
      else 
         if (((etot-pt(jtp,4)).gt.0.0001).or.((pztot-pt(jtp,3))
     1        .gt.0.0001).or.((pytot-pt(jtp,2)).gt.0.0001)
     2        .or.((pxtot-pt(jtp,1)).gt.0.0001)) then
            write(6,*) 'energy not conserved'   
         end if 
      end if 

      if (nstr.eq.0) then
         if (rlu(0).lt.0.5) then
            ihrd(1) = 0
            ihrd(2) = 1
         else
            ihrd(1) = 1
            ihrd(2) = 0
         end if
      else
         rdmn = 4.0*rlu(0)
         do 115, i = 1,4
            if (rdmn.gt.i-1.and.rdmn.lt.i) then
               ihrd(i) = 1
            else
               ihrd(i) = 0
            end if
 115     continue
      end if 
         
c     Record the flavor and 4-momentum of the strings and baryon.
      do 120, nb=1,2+nstr
         if (nb.eq.1) then
            kfb(nb,1) = kfvq(nb)
            kfb(nb,2) = kfsq(nb)
         else if (nb.eq.2) then
            kfb(nb,1) = kfsq(nb)
            kfb(nb,2) = kfvq(nb)
         else
            kfb(nb,1) = kfsq(nb)
            kfb(nb,2) = kfvq(nb)
         end if 
c     kfb(nb,3) states whether it is a baryon (0) or a string (1)
         kfb(nb,3) = 1

         kfb(nb,4) = ntp
         kfb(nb,5) = jtp
         kfb(nb,7) = ihrd(nb)
c      
         if (ihrd(nb).eq.1) then
            if (ntp.eq.1) then
               kfb(nb,8) = npj(jtp)
            else 
               kfb(nb,8) = ntj(jtp)
            end if 
         else
            kfb(nb,8) = 0
         endif
c     
         pb(nb,1) = px(nb)
         pb(nb,2) = py(nb)
         pb(nb,3) = pz(nb)
         pb(nb,4) = e(nb)
         pb(nb,5) = pm(nb)
         if (nb.eq.1) then
            pb(nb,6) = pxvq(nb)
            pb(nb,7) = pyvq(nb)
            pb(nb,8) = pxsq(nb)
            pb(nb,9) = pysq(nb)
            pb(nb,14) = vqm(nb)
            pb(nb,15) = sqm(nb)
         else if (nb.eq.2) then
            pb(nb,6) = pxsq(nb)
            pb(nb,7) = pysq(nb)
            pb(nb,8) = pxvq(nb)
            pb(nb,9) = pyvq(nb)
            pb(nb,14) = sqm(nb)
            pb(nb,15) = vqm(nb)
         else
            pb(nb,6) = pxsq(nb)
            pb(nb,7) = pysq(nb)
            pb(nb,8) = pxvq(nb)
            pb(nb,9) = pyvq(nb)
            pb(nb,14) = sqm(nb)
            pb(nb,15) = vqm(nb)
         end if 
 120  continue

      kfb(nb,1) = kfbar
      kfb(nb,2) = kfbar
c     kfb(nb,3) states whether the resonance is a baryon (0) or string (1)
      kfb(nb,3) = 0
      kfb(nb,4) = ntp
      kfb(nb,5) = jtp
      kfb(nb,7) = 0
      kfb(nb,8) = 0
c     
      pb(nb,1) = pxbar
      pb(nb,2) = pybar
      pb(nb,3) = pzb
      pb(nb,4) = eb
      pb(nb,5) = bmass
      
      nb = nb + 1
      kfb(nb,1) = kfabar
      kfb(nb,2) = kfabar
c     kfb(nb,3) states whether the resonance is a baryon (0) or string (1)
      kfb(nb,3) = 0
      kfb(nb,4) = ntp
      kfb(nb,5) = jtp
      kfb(nb,7) = 0
      kfb(nb,8) = 0
c     
      pb(nb,1) = pxabar
      pb(nb,2) = pyabar
      pb(nb,3) = pzab
      pb(nb,4) = eab
      pb(nb,5) = abmass

c     Write out the properties of the event.
      if (ib(5).eq.1) then
         if (ntp.eq.1) then 
            write (6,*) 'projectile string configuration  ---==--= '
         else
            write (6,*) 'target string configuration  =--==--- '
         end if 
         write(6,*) ' nb ', nb
         write(6,*) 'px(i),pxbar ',px(1),px(2),pxbar,pxabar
         write(6,*) 'py(i),pybar ',py(1),py(2),pybar,pyabar
         write(6,*) 'pz(i),pzbar ',pz(1),pz(2),pzb,pzab
         write(6,*) 'e(i),ebar   ',e(1),e(2),eb,eab
         write(6,*) 'y(i),ybar   ',y(1),y(2),yb,yab
         write(6,*) 'pm(i),bmass ',pm(1),pm(2),bmass,abmass 
         if (ntp.eq.1) then
            write(6,*) 'pp(jtp,4) ', pp(jtp,4),' etot ', etot
            write(6,*) 'pp(jtp,3) ', pp(jtp,3),' pztot ', pztot
            write(6,*) 'pp(jtp,2) ', pp(jtp,2),' pytot ', pytot
            write(6,*) 'pp(jtp,1) ', pp(jtp,1),' pxtot ', pxtot
         else 
            write(6,*) 'pt(jtp,4) ', pt(jtp,4),' etot ', etot
            write(6,*) 'pt(jtp,3) ', pt(jtp,3),' pztot ', pztot
            write(6,*) 'pt(jtp,2) ', pt(jtp,2),' pytot ', pytot
            write(6,*) 'pt(jtp,1) ', pt(jtp,1),' pxtot ', pxtot
         end if 
      endif
      
      return
      
      end
      
c     The subroutine HIJUNC models baryon stopping by creating the multi-
c     string configurations motivated by the Rossi-Veneziano baryon junction 
c     Regge exchange as calculated by Kharzeev.
c
c     Ref:  D. Kharzeev, Phys. Lett. B378 (1996) 238.
c
c     by S.E. Vance, 3/8/98
c     
c     
      SUBROUTINE HIJBJE(JTP,NTP,IERROR)

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/     
      COMMON/HISTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
      SAVE  /HISTRNG/
      COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &     PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &     PJPM(300,500),NTJ(300),KFTJ(300,500),
     &     PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &     PJTE(300,500),PJTM(300,500)
      SAVE  /HIJJET1/
      COMMON/HIJB/nb,ib(15),rb(15),sb(5),kb(10,900),kfb(10,15),
     &     pb(10,15),ncb(15)
      SAVE  /HIJB/

      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)   
      SAVE  /LUJETS/
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE  /LUDAT1/
c     
      EXTERNAL ULMASS, RAPMAX, KFBARYON, LUKFDI

      DIMENSION e(5),px(5),py(5),pz(5),pm(5),y(5)
c     The energy, p_x, p_y, p_z, mass and rapidity of the strings/baryon. 
      DIMENSION kfvq(5),emvq(5),epvq(5),pxvq(5),pyvq(5),vqmt(5),vqm(5)
c     The flavor, E+, E-, p_x, p_y, m_T and mass of the valence quarks.
      DIMENSION kfsq(5),emsq(5),epsq(5),pxsq(5),pysq(5),sqmt(5),sqm(5)
c     The flavor, E+, E-, p_x, p_y, m_T and mass of the sea quarks.
      DIMENSION yvq(5),ihrd(5),kmes(5),dyms(5)
c     Other various observables.
 
      IERROR = 0
      nbhit = 0
      nb = 0

c     set parameters
      ptwid = sb(3)
      itype = kb(ntp,jtp)
c
c     Determine if there exists enough phase space to create the 
c     configuration.
      
      if (ntp.eq.1) then
         smass = sqrt(max((pp(jtp,4))**2 - (pp(jtp,3))**2 
     1        - pp(jtp,2)**2 - pp(jtp,1)**2,0.0001))
      else
         smass = sqrt(max((pt(jtp,4))**2 - (pt(jtp,3))**2 
     1        - pt(jtp,2)**2 - pt(jtp,1)**2,0.0001))
      end if
      if (smass.lt.rb(4)) then
         ierror = 1
         nb = 1
         return
      end if 

c     Determine the rapidities of the quark and the diquark.
c     pmtzp:  mT of the +z moving particle 
c     pmtzn:  mT of the -z moving particle
c     
c     The constituent quark masses are used and are flagged by setting
c     mstj(93) = 1.
c

      if (ntp.eq.1) then
         kfq = min(iabs(nfp(jtp,1)),iabs(nfp(jtp,2)))
         kfqq = max(iabs(nfp(jtp,1)),iabs(nfp(jtp,2)))
         kfsign = nfp(jtp,1)/iabs(nfp(ntp,1))
         ep = pp(jtp,4) + pp(jtp,3)
         em = pp(jtp,4) - pp(jtp,3)
         pmtzp = sqrt(pp(jtp,15)**2 + pp(jtp,8)**2 + pp(jtp,9)**2)
         mstj(93) = 1
         pmtzn = sqrt(ulmass(kfq)**2 + pp(jtp,6)**2 + pp(jtp,7)**2)

      else if (ntp.eq.2) then
         kfq = min(iabs(nft(jtp,1)),iabs(nft(jtp,2)))
         kfqq = max(iabs(nft(jtp,1)),iabs(nft(jtp,2)))
         kfsign = nft(jtp,1)/iabs(nft(jtp,1))
         ep = pt(jtp,4) + pt(jtp,3)
         em = pt(jtp,4) - pt(jtp,3)
         pmtzn = sqrt(pt(jtp,15)**2 + pt(jtp,8)**2 + pt(jtp,9)**2)
         mstj(93) = 1
         pmtzp = sqrt(ulmass(kfq)**2 + pt(jtp,6)**2 + pt(jtp,7)**2)
      else 
         write(6,*) 'error in rapidity identification'
      end if

c     sample for the rapidity of the baryon 
c     iend choses which particles rapidity is solved for the +z or -z moving
c     particle
      if (ntp.eq.1) then
         iend = 1
      else 
         iend = 0
      end if 
      CALL RAPMAX(iend,ep,em,pmtzp,pmtzn,ymax,ymin,irperr)
      if (irperr.eq.1) then
c         write(6,*) 'ymin,ymax',ymin,ymax
c         write(6,*) 'mass',sqrt(ep*em)
         nbhit = 200
      end if 
       
      if (ymax.lt.ymin) then
         nbhit = 200
         write(6,*) 'BJE Failed, initial string is too short'
      end if

c     Determine the flavors of the valence quarks. 
      kfvq(1) = kfq
      if (kfqq.eq.1103) then
         kfvq(2) = 1
         kfvq(3) = 1
      else if (kfqq.eq.2203) then
         kfvq(2) = 2
         kfvq(3) = 2
      else if (kfqq.eq.2103.or.kfqq.eq.2101) then
         if (rlu(0).lt.0.5) then
            kfvq(2) = 1
            kfvq(3) = 2
         else
            kfvq(2) = 2
            kfvq(3) = 1
         endif 
      endif
      
c     Determine the flavors of the sea quarks and the resolved baryon. 
      CALL KFBARYON(3, kfbq1, kfbq2, kfbq3, kfbar)

      if (kfsign.eq.-1) then
         kfvq(1) = isign(kfvq(1),kfsign)            
         kfvq(2) = isign(kfvq(2),kfsign)
         kfvq(3) = isign(kfvq(3),kfsign)
c         
         kfbar   = isign(kfbar,kfsign)
         kfsq(1) = -isign(kfbq1,kfsign)
         kfsq(2) = -isign(kfbq2,kfsign)
         kfsq(3) = -isign(kfbq3,kfsign)
      else
         kfsq(1) = -kfbq1
         kfsq(2) = -kfbq2
         kfsq(3) = -kfbq3
      end if 
      
 10   nbhit = nbhit + 1
      if (nbhit.gt.70) then
         ierror = 2
         if (ib(5).eq.1) write(6,*) 'kfbar',kfbar
         nb = 1
         return
      endif
      if (mod(nbhit,5).ne.1) goto 33

c     Determine the pT and the mT of the sea quarks, the valence quarks
c     and the baryon according to a gaussian pT distribution.
      do 30, i = 1,3
         ptmag = ptwid*sqrt(-log(max(1E-10,(1-rlu(0)))))
         ptphi = 2.0*hipr1(40)*rlu(0)
         pxsq(i) = ptmag*cos(ptphi)
         pysq(i) = ptmag*sin(ptphi)
 30   continue
      if (ib(4).eq.1.and.ntp.eq.1) then
         pxbar = -1.0*(pxsq(1) + pxsq(2) + pxsq(3)) + PP(jtp,8)
         pybar = -1.0*(pysq(1) + pysq(2) + pysq(3)) + PP(jtp,9)
      else if (ib(4).eq.1.and.ntp.eq.2) then
         pxbar = -1.0*(pxsq(1) + pxsq(2) + pxsq(3)) + PT(jtp,8)
         pybar = -1.0*(pysq(1) + pysq(2) + pysq(3)) + PT(jtp,9)
      else
         pxbar = -1.0*(pxsq(1) + pxsq(2) + pxsq(3)) 
         pybar = -1.0*(pysq(1) + pysq(2) + pysq(3)) 
      end if 

      bmass  = ulmass(kfbar)
      ptbar  = sqrt(pxbar**2 + pybar**2) 
      barmt  = sqrt(bmass**2 + ptbar**2)

 33   continue
c     Determine the rapidity of the produced baryon.
c     itype = 2:   the cosh(y/2) dependence
c     itype = 3:   the constant rapidity dependence

      if (itype.eq.2) then
c     Sample from a Exp(abj*y) rapidity distribution.
         if (ntp.eq.1) then
            abj = 0.5
         else if (ntp.eq.2) then
            abj = -0.5
         end if 
         ybjmax = ymax - 1.5
         ybjmin = ymin + 1.5

         pnorm = exp(abj*ybjmax)-exp(abj*ybjmin)
         ybar  = 1/abj*log(exp(abj*ybjmin) + pnorm*rlu(0))

         if (ib(5).eq.1) then
            write(6,*) 'pnorm,ybar',pnorm,ybar
         end if 

         dystr = abs(ymax - ymin)
         yhalf = ymin + 0.5*dystr
   
      else if (itype.eq.3) then
c     Sample from a uniform rapidity distribution.
         dystr = abs(ymax - ymin)
         yhalf = ymin + 0.5*dystr

         ybjmax = yhalf + 0.25*dystr
         ybjmin = yhalf - 0.25*dystr
    
         ybar  = ybjmin + (ybjmax - ybjmin)*rlu(0)

      else
         IERROR = 1
         write(6,*) 'rapidity dependence of BJ not defined'
         write(6,*) 'itype',itype
	 nb = 1
         return
      end if 

      do 40, i = 1,3
         if (i.eq.1) then
            if (ntp.eq.1) then
               pxvq(i) = PP(jtp,6)
               pyvq(i) = PP(jtp,7)
            else 
               pxvq(i) = PT(jtp,6)
               pyvq(i) = PT(jtp,7)
            end if 
         else 
            if (ib(4).eq.0.and.ntp.eq.1) then
               pxvq(i) = PP(jtp,8)/2.0
               pyvq(i) = PP(jtp,9)/2.0
            else if (ib(4).eq.0.and.ntp.eq.2) then
               pxvq(i) = PT(jtp,8)/2.0
               pyvq(i) = PT(jtp,9)/2.0
            else
               pxvq(i) = 0.0
               pyvq(i) = 0.0
            end if 
         end if 
         
         mstj(93) = 1
         sqm(i)   = ulmass(kfsq(i))
         sqmt(i)  = sqrt(sqm(i)**2 + pxsq(i)**2 + pysq(i)**2)
         epsq(i)  = sqmt(i)*exp(ybar)
         emsq(i)  = sqmt(i)*exp(-ybar)

         mstj(93) = 1
         vqm(i)   = ulmass(kfvq(i))
         vqmt(i)  = sqrt(vqm(i)**2 + pxvq(i)**2 + pyvq(i)**2)

 40   continue

c     Determine the kinematic properties of the baryon.
      pzbar  = barmt*sinh(ybar)
      ebar   = barmt*cosh(ybar) 
      epbar  = ebar + pzbar
      embar  = ebar - pzbar
      if (ebar**2-pxbar**2-pybar**2-pzbar**2.le.0) goto 10
      
c     Determine the E+ and E- available for string creation.
      epr = ep - epbar - epsq(1) - epsq(2) - epsq(3)
      emr = em - embar - emsq(1) - emsq(2) - emsq(3)          
      if (epr.lt.0.0.or.emr.lt.0.0) goto 10

      vqmt23 = vqmt(2) + vqmt(3)
      if (ntp.eq.1) then
         iend = 1
         CALL RAPMAX(iend,epr,emr,vqmt23,vqmt(1),yrmax,yrmin,
     1        irperr)
         if (irperr.eq.1) goto 10
         yvq(1) = yrmin
      else 
         iend = 0
         CALL RAPMAX(iend,epr,emr,vqmt(1),vqmt23,yrmax,yrmin,
     1        irperr)
         if (irperr.eq.1) goto 10
         yvq(1) = yrmax
      end if 
      
      epvq(1) = vqmt(1)*exp(yvq(1))
      emvq(1) = vqmt(1)*exp(-yvq(1))

      epeff = ep - epbar - epsq(1) - epvq(1) - epsq(2) - epsq(3) 
      emeff = em - embar - emsq(1) - emvq(1) - emsq(2) - emsq(3)
      
      if (epeff.le.0.0.or.emeff.le.0.0) goto 10
      
      epvq(2) = epeff/2.0
      epvq(3) = epeff/2.0
      
      emvq(2) = emeff/2.0
      emvq(3) = emeff/2.0
      
      yvq(2) = 0.5*log(epvq(2)/emvq(2))
      yvq(3) = 0.5*log(epvq(3)/emvq(3))
      
      if ((ntp.eq.1.and.yvq(2).le.ybar).or.
     &     (ntp.eq.1.and.yvq(3).le.ybar).or.
     &     (ntp.eq.2.and.yvq(2).ge.ybar).or.
     &     (ntp.eq.2.and.yvq(3).ge.ybar)) goto 10
      
      epvq(1) = ep-epbar-epsq(1)-epsq(2)-epvq(2)-epsq(3)-epvq(3) 
      emvq(1) = em-embar-emsq(1)-emsq(2)-emvq(2)-emsq(3)-emvq(3)
      if (epvq(1).le.0.0.or.emvq(1).le.0.0) goto 10
      
      yvq(1) = 0.5*log(epvq(1)/emvq(1)) 
      if ((ntp.eq.1.and.yvq(1).ge.ybar).or.
     &     (ntp.eq.2.and.yvq(1).le.ybar)) goto 10

c     Determine the E, Pz, Px, Py, y and m of the strings and baryon.
      do 100, i = 1,3
         e(i) = ((epsq(i)+emsq(i))+(epvq(i)+emvq(i)))/2.0
         pz(i) = ((epsq(i)-emsq(i))+(epvq(i)-emvq(i)))/2.0
         if (pz(i).gt.e(i)) goto 10
         if (((epsq(i)+epvq(i))/(emsq(i)+emvq(i))).lt.0.0) goto 10
         y(i) = 1./2.*log((epsq(i)+epvq(i))/(emsq(i)
     &        +emvq(i)))
         px(i) = pxsq(i) + pxvq(i)
         py(i) = pysq(i) + pyvq(i) 
         if ((e(i)**2-px(i)**2-py(i)**2-pz(i)**2).le.0.0) goto 10
         pm(i) = sqrt(e(i)**2 - px(i)**2 - py(i)**2 - pz(i)**2)
 100  continue
c
c     require that the string fragment correctly
c     
      do 110, i = 1,3
         amt = pm(i)**2 + px(i)**2 + py(i)**2
         amt1 = vqmt(i)**2
         amt2 = sqmt(i)**2
         pzcm2 = (AMT**2+AMT1**2+AMT2**2-2.0*AMT*AMT1
     &        -2.0*AMT*AMT2-2.0*AMT1*AMT2)/4.0/AMT
         if (pzcm2.le.0) goto 10
 110  continue
      

c     Check for energy conservation
      pxtot = px(1) + px(2) + px(3) + pxbar 
      pytot = py(1) + py(2) + py(3) + pybar 
      pztot = pz(1) + pz(2) + pz(3) + pzbar 
      etot  = e(1)  + e(2)  + e(3)  + ebar 
      
      if (ntp.eq.1) then
         if (((etot-pp(jtp,4)).gt.0.0001).or.((pztot-pp(jtp,3))
     &        .gt.0.0001).or.((pytot-pp(jtp,2)).gt.0.0001)
     &        .or.((pxtot-pp(jtp,1)).gt.0.0001)) then
            write(6,*) 'energy not conserved'
         end if 
      else 
         if (((etot-pt(jtp,4)).gt.0.0001).or.((pztot-pt(jtp,3))
     &        .gt.0.0001).or.((pytot-pt(jtp,2)).gt.0.0001)
     &        .or.((pxtot-pt(jtp,1)).gt.0.0001)) then
            write(6,*) 'energy not conserved'
         end if 
      end if 

c     Determine which of the strings will have the hard interaction
      if (rlu(0).lt.0.33333) then
         ihrd(1) = 1
         ihrd(2) = 0
         ihrd(3) = 0
      else
         ihrd(1) = 0
         if (rlu(0).lt.0.5) then
            ihrd(2) = 0
            ihrd(3) = 1
         else
            ihrd(2) = 1
            ihrd(3) = 0
         end if 
      end if
      
c     Record the flavor and 4-momentum of the strings and baryon.
      do 120, nb=1,3
         if ((ntp.eq.1.and.nb.eq.1).or.(ntp.eq.2.and.nb.ne.1)) then
            kfb(nb,1) = kfsq(nb)
            kfb(nb,2) = kfvq(nb)
         else if ((ntp.eq.1.and.nb.ne.1).or.(ntp.eq.2.and.nb.eq.1)) then
            kfb(nb,1) = kfvq(nb)
            kfb(nb,2) = kfsq(nb)
         end if
c     kfb(nb,3) states whether it is a baryon (0) or a string (1)
         kfb(nb,3) = 1
         kfb(nb,4) = ntp
         kfb(nb,5) = jtp
         kfb(nb,7) = ihrd(nb)
c      
         if (ihrd(nb).eq.1) then
            if (ntp.eq.1) then
               kfb(nb,8) = npj(jtp)
            else 
               kfb(nb,8) = ntj(jtp)
            end if 
         else
            kfb(nb,8) = 0
         endif
c     
         pb(nb,1) = px(nb)
         pb(nb,2) = py(nb)
         pb(nb,3) = pz(nb)
         pb(nb,4) = e(nb)
         pb(nb,5) = pm(nb)
         if ((ntp.eq.1.and.nb.eq.1).or.(ntp.eq.2.and.nb.ne.1)) then
            pb(nb,6) = pxsq(nb)
            pb(nb,7) = pysq(nb)
            pb(nb,8) = pxvq(nb)
            pb(nb,9) = pyvq(nb)
            pb(nb,14) = sqm(nb)
            pb(nb,15) = vqm(nb)
         else if ((ntp.eq.1.and.nb.ne.1).or.(ntp.eq.2.and.nb.eq.1)) then
            pb(nb,6) = pxvq(nb)
            pb(nb,7) = pyvq(nb)
            pb(nb,8) = pxsq(nb)
            pb(nb,9) = pysq(nb)
            pb(nb,14) = vqm(nb)
            pb(nb,15) = sqm(nb)
         end if 
 120  continue

c     nb = 4 now
      kfb(nb,1) = kfbar
      kfb(nb,2) = kfbar
c     kfb(nb,3) states whether the resonance is a baryon (0) or string (1)
      kfb(nb,3) = 0
      kfb(nb,4) = ntp
      kfb(nb,5) = jtp
      kfb(nb,7) = 0
      kfb(nb,8) = 0
c     
      pb(nb,1) = pxbar
      pb(nb,2) = pybar
      pb(nb,3) = pzbar
      pb(nb,4) = ebar
      pb(nb,5) = sqrt(ebar**2-pxbar**2-pybar**2-pzbar**2)
      
c     Write out the properties of the event.
      if (ib(5).eq.1) then
         if (ntp.eq.1) then 
            write (6,*) 'projectile string configuration  ---=== '
         else
            write (6,*) 'target string configuration  ===--- '
         end if 
         write(6,*) ' nb ', nb
         write(6,*) 'px(i),pxbar ',px(1),px(2),px(3),pxbar
         write(6,*) 'py(i),pybar ',py(1),py(2),py(3),pybar
         write(6,*) 'pz(i),pzbar ',pz(1),pz(2),pz(3),pzbar
         write(6,*) 'e(i),ebar   ',e(1),e(2),e(3),ebar
         write(6,*) 'y(i),ybar   ',y(1),y(2),y(3),ybar
         write(6,*) 'pm(i),bmass ',pm(1),pm(2),pm(3),bmass 
         if (ntp.eq.1) then
            write(6,*) 'pp(jtp,4) ', pp(jtp,4),' etot ', etot
            write(6,*) 'pp(jtp,3) ', pp(jtp,3),' pztot ', pztot
            write(6,*) 'pp(jtp,2) ', pp(jtp,2),' pytot ', pytot
            write(6,*) 'pp(jtp,1) ', pp(jtp,1),' pxtot ', pxtot
         else 
            write(6,*) 'pt(jtp,4) ', pt(jtp,4),' etot ', etot
            write(6,*) 'pt(jtp,3) ', pt(jtp,3),' pztot ', pztot
            write(6,*) 'pt(jtp,2) ', pt(jtp,2),' pytot ', pytot
            write(6,*) 'pt(jtp,1) ', pt(jtp,1),' pxtot ', pxtot
         end if 
      endif
      
      return
      
      end
      

c     The subroutine HIJDQB models baryon number transport by the
c     breaking of a diquark as mediated by a gluon exchange. 
c
c     Ref:  B.Z. Kopeliovich and B.G. Zakharov, Z. Phys. C43 (1989) 241.
c
c     by S.E. Vance, Aug 1998
c     
c     
      SUBROUTINE HIJDQB(JTP,NTP,IERROR)

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/     
      COMMON/HISTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
      SAVE  /HISTRNG/
      COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &     PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &     PJPM(300,500),NTJ(300),KFTJ(300,500),
     &     PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &     PJTE(300,500),PJTM(300,500)
      SAVE  /HIJJET1/
      COMMON/HIJB/nb,ib(15),rb(15),sb(5),kb(10,900),kfb(10,15),
     &     pb(10,15),ncb(15)
      SAVE  /HIJB/

      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)   
      SAVE  /LUJETS/
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE  /LUDAT1/
c     
      EXTERNAL ULMASS, KFBARYON, LUKFDI
      
      DIMENSION e(5),px(5),py(5),pz(5),pm(5),y(5)
c     The energy, p_x, p_y, p_z, mass and rapidity of the strings/baryon. 
      DIMENSION kfvq(5),emvq(5),epvq(5),pxvq(5),pyvq(5),vqmt(5),vqm(5)
c     The flavor, E+, E-, p_x, p_y, m_T and mass of the valence quarks.
      DIMENSION kfsq(5),emsq(5),epsq(5),pxsq(5),pysq(5),sqmt(5),sqm(5)
c     The flavor, E+, E-, p_x, p_y, m_T and mass of the sea quarks.
      DIMENSION yvq(5),ihrd(5),kmes(5),dyms(5)
c     Other various observables.
 
      IERROR = 0
      nbhit = 0
      nb = 0

c     set parameters
      ptwid = sb(3)
c
c     Determine if there exists enough phase space to create the 
c     configuration.
c      
      if (ntp.eq.1) then
         smass = sqrt(max((pp(jtp,4))**2 - (pp(jtp,3))**2 
     1        - pp(jtp,2)**2 - pp(jtp,1)**2,0.0001))
      else
         smass = sqrt(max((pt(jtp,4))**2 - (pt(jtp,3))**2 
     1        - pt(jtp,2)**2 - pt(jtp,1)**2,0.0001))
      end if
      if (smass.lt.rb(12)) then
         ierror = 1
         nb = 1
         return
      end if 
         
c     Determine the rapidities of the quark and the diquark.
c     pmtzp:  mT of the +z moving particle 
c     pmtzn:  mT of the -z moving particle
c     
c     The constituent quark masses are used and are flagged by setting
c     mstj(93) = 1.
c
      if (ntp.eq.1) then
         kfq = min(iabs(nfp(jtp,1)),iabs(nfp(jtp,2)))
         kfqq = max(iabs(nfp(jtp,1)),iabs(nfp(jtp,2)))
         kfsign = nfp(jtp,1)/iabs(nfp(ntp,1))
         ep = pp(jtp,4) + pp(jtp,3)
         em = pp(jtp,4) - pp(jtp,3)
         pmtzp = sqrt(pp(jtp,15)**2 + pp(jtp,8)**2 + pp(jtp,9)**2)
         mstj(93) = 1
         pmtzn = sqrt(ulmass(kfq)**2 + pp(jtp,6)**2 + pp(jtp,7)**2)

      else if (ntp.eq.2) then
         kfq = min(iabs(nft(jtp,1)),iabs(nft(jtp,2)))
         kfqq = max(iabs(nft(jtp,1)),iabs(nft(jtp,2)))
         kfsign = nft(jtp,1)/iabs(nft(jtp,1))
         ep = pt(jtp,4) + pt(jtp,3)
         em = pt(jtp,4) - pt(jtp,3)
         pmtzn = sqrt(pt(jtp,15)**2 + pt(jtp,8)**2 + pt(jtp,9)**2)
         mstj(93) = 1
         pmtzp = sqrt(ulmass(kfq)**2 + pt(jtp,6)**2 + pt(jtp,7)**2)
      else 
         write(6,*) 'error in rapidity identification'
      end if

c     sample for the rapidity of the baryon 
c     iend choses which particle's rapidity is solved for the +z or -z 
c     moving particle
      if (ntp.eq.1) then
         iend = 1
      else 
         iend = 0
      end if 
      CALL RAPMAX(iend,ep,em,pmtzp,pmtzn,ymax,ymin,irperr)
      if (irperr.eq.1) then
         nbhit = 200
      end if 
       
      if (ymax.lt.ymin) then
         nbhit = 200
         write(6,*) 'BJ Failed, initial string is too short to break'
      end if

c     Determine the flavors of the valence quarks. 
c     The first valence quark is the trailing quark of the jet.
c     The second quark is the valence quark found in the baryon.
c     The third valence quark is the leading quark of the jet. 
      kfvq(1) = kfq
      if (kfqq.eq.1103) then
         kfvq(2) = 1
         kfvq(3) = 1
      else if (kfqq.eq.2203) then
         kfvq(2) = 2
         kfvq(3) = 2
      else if (kfqq.eq.2103.or.kfqq.eq.2101) then
         if (rlu(0).lt.0.5) then
            kfvq(2) = 1
            kfvq(3) = 2
         else
            kfvq(2) = 2
            kfvq(3) = 1
         endif 
      endif
      
      kfbq1 = kfvq(2)
c     Determine the flavors of the two sea quarks and the resolved baryon.  
      CALL KFBARYON(2, kfbq1, kfbq2, kfbq3, kfbar)

      if (kfsign.eq.-1) then
         kfvq(1) = isign(kfvq(1),kfsign)            
         kfvq(2) = isign(kfvq(2),kfsign)
         kfvq(3) = isign(kfvq(3),kfsign)
c         
         kfbar   = isign(kfbar,kfsign)
         kfsq(1) = -isign(kfbq2,kfsign)
         kfsq(3) = -isign(kfbq3,kfsign)
      else
         kfsq(1) = -kfbq2
         kfsq(3) = -kfbq3
      end if 

 10   nbhit = nbhit + 1
      if (nbhit.gt.70) then
         IERROR = 2
         nb = 1
         if (ib(5).eq.1) write(6,*) 'kfbar',kfbar
         return
      endif


c     Sample from a Exp(abj*y) rapidity distribution for the baryon.
      if (ntp.eq.1) then
         abj = 0.5
      else if (ntp.eq.2) then
         abj = -0.5
      end if 
      ybjmax = ymax - 1.5
      ybjmin = ymin + 1.5
      
      pnorm = exp(abj*ybjmax)-exp(abj*ybjmin)
      ybar  = 1/abj*log(exp(abj*ybjmin) + pnorm*rlu(0))
      
      if (ib(5).eq.1) then
         write(6,*) 'pnorm,ybar',pnorm,ybar
      end if 
      

c     Determine the pT and the mT of the sea quarks, the valence quarks
c     and the baryon according to a gaussian pT distribution.
      do 30, i = 1, 3, 2
         ptmag = ptwid*sqrt(-log(max(1E-10,(1-rlu(0)))))
         ptphi = 2.0*hipr1(40)*rlu(0)
         pxsq(i) = ptmag*cos(ptphi)
         pysq(i) = ptmag*sin(ptphi)
 30   continue
      if (ib(4).eq.1.and.ntp.eq.1) then
         pxbar = -1.0*(pxsq(1) + pxsq(3)) + PP(jtp,8)
         pybar = -1.0*(pysq(1) + pysq(3)) + PP(jtp,9)
      else if (ib(4).eq.1.and.ntp.eq.2) then
         pxbar = -1.0*(pxsq(1) + pxsq(3)) + PT(jtp,8)
         pybar = -1.0*(pysq(1) + pysq(3)) + PT(jtp,9)
      else
         pxbar = -1.0*(pxsq(1) + pxsq(3)) 
         pybar = -1.0*(pysq(1) + pysq(3)) 
      end if 

      bmass  = ulmass(kfbar)
      ptbar  = sqrt(pxbar**2 + pybar**2) 
      barmt  = sqrt(bmass**2 + ptbar**2)

      do 40, i = 1, 3, 2
         if (i.eq.1) then
            if (ntp.eq.1) then
               pxvq(i) = PP(jtp,6)
               pyvq(i) = PP(jtp,7)
            else 
               pxvq(i) = PT(jtp,6)
               pyvq(i) = PT(jtp,7)
            end if 
         else 
            if (ib(4).eq.0.and.ntp.eq.1) then
               pxvq(i) = PP(jtp,8)
               pyvq(i) = PP(jtp,9)
            else if (ib(4).eq.0.and.ntp.eq.2) then
               pxvq(i) = PT(jtp,8)
               pyvq(i) = PT(jtp,9)
            else
               pxvq(i) = 0.0
               pyvq(i) = 0.0
            end if 
         end if 
         
         mstj(93) = 1
         sqm(i)   = ulmass(kfsq(i))
         sqmt(i)  = sqrt(sqm(i)**2 + pxsq(i)**2 + pysq(i)**2)
         epsq(i)  = sqmt(i)*exp(ybar)
         emsq(i)  = sqmt(i)*exp(-ybar)

         mstj(93) = 1
         vqm(i)   = ulmass(kfvq(i))
         vqmt(i)  = sqrt(vqm(i)**2 + pxvq(i)**2 + pyvq(i)**2)

 40   continue

c     Determine the kinematic properties of the baryon.
      pzbar  = barmt*sinh(ybar)
      ebar   = barmt*cosh(ybar) 
      epbar  = ebar + pzbar
      embar  = ebar - pzbar
      if (ebar**2-pxbar**2-pybar**2-pzbar**2.le.0) goto 10
      
c     Determine the E+ and E- available for string creation.
      epr = ep - epbar - epsq(1) - epsq(3)
      emr = em - embar - emsq(1) - emsq(3)          
      if (epr.lt.0.0.or.emr.lt.0.0) goto 10

      if (ntp.eq.1) then
         iend = 1
         CALL RAPMAX(iend,epr,emr,vqmt(3),vqmt(1),yrmax,yrmin,irperr)
         if (irperr.eq.1) goto 10
         yvq(1) = yrmin
      else 
         iend = 0
         CALL RAPMAX(iend,epr,emr,vqmt(1),vqmt(3),yrmax,yrmin,irperr)
         if (irperr.eq.1) goto 10
         yvq(1) = yrmax
      end if 
      
      epvq(1) = vqmt(1)*exp(yvq(1))
      emvq(1) = vqmt(1)*exp(-yvq(1))

      epvq(3) = ep - epbar - epsq(1) - epvq(1) - epsq(3) 
      emvq(3) = em - embar - emsq(1) - emvq(1) - emsq(3)
      
      if (epvq(3).le.0.0.or.emvq(3).le.0.0) goto 10
      
      yvq(3) = 0.5*log(epvq(3)/emvq(3))
      
      if ((ntp.eq.1.and.yvq(3).le.ybar).or.
     &     (ntp.eq.2.and.yvq(3).ge.ybar)) goto 10
      
      yvq(1) = 0.5*log(epvq(1)/emvq(1)) 

      if ((ntp.eq.1.and.yvq(1).ge.ybar).or.
     &     (ntp.eq.2.and.yvq(1).le.ybar)) goto 10

c     Determine the E, Pz, Px, Py, y and m of the strings and baryon.
      do 100, i = 1,3,2
         e(i) = ((epsq(i)+emsq(i))+(epvq(i)+emvq(i)))/2.0
         pz(i) = ((epsq(i)-emsq(i))+(epvq(i)-emvq(i)))/2.0
         if (pz(i).gt.e(i)) goto 10
         if (((epsq(i)+epvq(i))/(emsq(i)+emvq(i))).lt.0.0) goto 10
         y(i) = 1./2.*log((epsq(i)+epvq(i))/(emsq(i)
     &        +emvq(i)))
         px(i) = pxsq(i) + pxvq(i)
         py(i) = pysq(i) + pyvq(i) 
         if ((e(i)**2-px(i)**2-py(i)**2-pz(i)**2).le.0.0) goto 10
         pm(i) = sqrt(e(i)**2 - px(i)**2 - py(i)**2 - pz(i)**2)
 100  continue
      
c     Check for energy conservation
      pxtot = px(1) + px(3) + pxbar 
      pytot = py(1) + py(3) + pybar 
      pztot = pz(1) + pz(3) + pzbar 
      etot  = e(1)  + e(3)  + ebar 
      
      if (ntp.eq.1) then
         if (((etot-pp(jtp,4)).gt.0.0001).or.((pztot-pp(jtp,3))
     &        .gt.0.0001).or.((pytot-pp(jtp,2)).gt.0.0001)
     &        .or.((pxtot-pp(jtp,1)).gt.0.0001)) then
            write(6,*) 'energy not conserved'
         end if 
      else 
         if (((etot-pt(jtp,4)).gt.0.0001).or.((pztot-pt(jtp,3))
     &        .gt.0.0001).or.((pytot-pt(jtp,2)).gt.0.0001)
     &        .or.((pxtot-pt(jtp,1)).gt.0.0001)) then
            write(6,*) 'energy not conserved'   
         end if 
      end if 

c     Determine which of the strings will have the hard interaction
      if (rlu(0).lt.0.5) then
         ihrd(1) = 0
         ihrd(3) = 1
      else
         ihrd(1) = 1
         ihrd(3) = 0
      end if 
      
c     Record the flavor and 4-momentum of the strings and baryon.
      do 120, i = 1,3,2
         nbary = nbary + 1
         if ((ntp.eq.1.and.i.eq.1).or.(ntp.eq.2.and.i.ne.1)) then
            kfb(i,1) = kfsq(i)
            kfb(i,2) = kfvq(i)
         else if ((ntp.eq.1.and.i.ne.1).or.(ntp.eq.2.and.i.eq.1)) then
            kfb(i,1) = kfvq(i)
            kfb(i,2) = kfsq(i)
         end if
c     kfb(i,3) states whether it is a baryon (0) or a string (1)
         kfb(i,3) = 1
         kfb(i,4) = ntp
         kfb(i,5) = jtp
         kfb(i,7) = ihrd(i)
c      
         if (ihrd(i).eq.1) then
            if (ntp.eq.1) then
               kfb(i,8) = npj(jtp)
            else 
               kfb(i,8) = ntj(jtp)
            end if 
         else
            kfb(i,8) = 0
         endif
c     
         pb(i,1) = px(i)
         pb(i,2) = py(i)
         pb(i,3) = pz(i)
         pb(i,4) = e(i)
         pb(i,5) = pm(i)
         if ((ntp.eq.1.and.i.eq.1).or.(ntp.eq.2.and.i.ne.1)) then
            pb(i,6) = pxsq(i)
            pb(i,7) = pysq(i)
            pb(i,8) = pxvq(i)
            pb(i,9) = pyvq(i)
            pb(i,14) = sqm(i)
            pb(i,15) = vqm(i)
         else if ((ntp.eq.1.and.i.ne.1).or.(ntp.eq.2.and.i.eq.1)) then
            pb(i,6) = pxvq(i)
            pb(i,7) = pyvq(i)
            pb(i,8) = pxsq(i)
            pb(i,9) = pysq(i)
            pb(i,14) = vqm(i)
            pb(i,15) = sqm(i)
         end if 
 120  continue

      nbary = nbary + 1
      kfb(2,1) = kfbar
      kfb(2,2) = kfbar
c     kfb(2,3) states whether the resonance is a baryon (0) or string (1)
      kfb(2,3) = 0
      kfb(2,4) = ntp
      kfb(2,5) = jtp
      kfb(2,7) = 0
      kfb(2,8) = 0
c     
      pb(2,1) = pxbar
      pb(2,2) = pybar
      pb(2,3) = pzbar
      pb(2,4) = ebar
      pb(2,5) = sqrt(ebar**2-pxbar**2-pybar**2-pzbar**2)
      
c     Write out the properties of the event.
      if (ib(5).eq.1) then
         if (ntp.eq.1) then 
            write (6,*) 'projectile string configuration  ---*-- '
         else
            write (6,*) 'target string configuration  --*--- '
         end if 
         write(6,*) ' nb ', nb
         write(6,*) 'px(i),pxbar ',px(1),px(3),pxbar
         write(6,*) 'py(i),pybar ',py(1),py(3),pybar
         write(6,*) 'pz(i),pzbar ',pz(1),pz(3),pzbar
         write(6,*) 'e(i),ebar   ',e(1),e(3),ebar
         write(6,*) 'y(i),ybar   ',y(1),y(3),ybar
         write(6,*) 'pm(i),bmass ',pm(1),pm(3),bmass 
         if (ntp.eq.1) then
            write(6,*) 'pp(jtp,4) ', pp(jtp,4),' etot ', etot
            write(6,*) 'pp(jtp,3) ', pp(jtp,3),' pztot ', pztot
            write(6,*) 'pp(jtp,2) ', pp(jtp,2),' pytot ', pytot
            write(6,*) 'pp(jtp,1) ', pp(jtp,1),' pxtot ', pxtot
         else 
            write(6,*) 'pt(jtp,4) ', pt(jtp,4),' etot ', etot
            write(6,*) 'pt(jtp,3) ', pt(jtp,3),' pztot ', pztot
            write(6,*) 'pt(jtp,2) ', pt(jtp,2),' pytot ', pytot
            write(6,*) 'pt(jtp,1) ', pt(jtp,1),' pxtot ', pxtot
         end if 
      endif
      
      return
      
      end
      

      SUBROUTINE RAPMAX(iend,ep,em,zpmt,znmt,ymax,ymin,irperr)
c
c     Using the total E+,E- and the mt of the end partons, this subroutine
c     determines the rapidities of the end partons of the string. 
c      
      double precision a,b,c,zpep,znep,zpy,zny
      
      irperr = 0
      if (ep.eq.0.0.or.em.eq.0.0.or.zpmt.eq.0.0.or.znmt.eq.0.0) then
         write(6,*) 'RAPMAX ERROR, unspecified values'
         return
      end if 

      if (iend.eq.0) then
         a = em
         b = znmt**2 - zpmt**2 - ep*em
         c = ep*(zpmt**2)
      else 
         a = em
         b = zpmt**2 - znmt**2 - ep*em
         c = ep*(znmt**2)
      end if 
         
      if ((b**2 - 4.0*a*c).lt.0.0) then
         irperr = 1
         return
      end if 

      if (iend.eq.0) then
         zpep = (-b + sqrt(b**2 - 4.0*a*c))/(2.0*a)
         znep = ep - zpep
      else
         znep = (-b - sqrt(b**2 - 4.0*a*c))/(2.0*a)
         zpep = ep - znep
      end if 
      
      if (zpep.lt.0.0.or.znep.lt.0.0) then
         irperr = 1
         return
      end if 

      zpy = dlog(zpep/zpmt)  
      zny = dlog(znep/znmt) 

      ymax = zpy
      ymin = zny

      return
      end

      SUBROUTINE KFBARYON(nq,kfq1,kfq2,kfq3,kfbar)
c
c     This subroutine is a modified version of the LUND subroutine
c     LUKFDI which produces the flavor contents of baryons and mesons
c     
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE  /LUDAT1/
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/
      COMMON/HIJB/nb,ib(15),rb(15),sb(5),kb(10,900),kfb(10,15),
     &     pb(10,15),ncb(15)
      SAVE  /HIJB/

      par2 = sb(1)
      par3 = 1.0
      par4 = 3.0*sb(2)
      kfl1 = 1
      kf2a = 0

c     note that nq=0 may cause problems with an infinite loop if 
c     the parameters are not chosen correctly.
      if (nq.eq.0) then
         kfla = kfq1
         kflb = kfq2
         kflc = kfq3
      else if (nq.eq.1) then
         kfla = kfq1
         kflb = kfq2
      else if (nq.eq.2) then
         kfla = kfq1
      else if (nq.eq.3) then
         KFLA = 1+INT((2.+PAR2*PAR3)*RLU(0))
      end if 
      
 130  continue
      if (nq.eq.1) then
         KFLC=1+INT((2.+PAR2*PAR3)*RLU(0))
      else if (nq.gt.2) then
         KFLB=1+INT((2.+PAR2*PAR3)*RLU(0)) 
         KFLC=1+INT((2.+PAR2*PAR3)*RLU(0)) 
      end if 
         
      KFLDS=1   
      IF(KFLB.GE.KFLC) KFLDS=3  
      IF(KFLDS.EQ.1.AND.PAR4*RLU(0).GT.1.) GOTO 130 
      IF(KFLDS.EQ.3.AND.PAR4.LT.RLU(0)) GOTO 130    
      KFL3=ISIGN(1000*MAX(KFLB,KFLC)+100*MIN(KFLB,KFLC)+KFLDS,KFL1) 
      
C...SU(6) factors for formation of baryon. Try again if fails.  
      KBARY=KFLDS 
      IF(KFLDS.EQ.3.AND.KFLB.NE.KFLC) KBARY=5 
      IF(KFLA.NE.KFLB.AND.KFLA.NE.KFLC) KBARY=KBARY+1 
      WT=PARF(60+KBARY)+PARJ(18)*PARF(70+KBARY)   
      IF(KF2A.EQ.0.AND.WT.LT.RLU(0)) GOTO 130 
    
C...Form baryon. Distinguish Lambda- and Sigmalike baryons. 
      KFLD=MAX(KFLA,KFLB,KFLC)    
      KFLF=MIN(KFLA,KFLB,KFLC)    
      KFLE=KFLA+KFLB+KFLC-KFLD-KFLF   
      KFLS=2  
      IF((PARF(60+KBARY)+PARJ(18)*PARF(70+KBARY))*RLU(0).GT.  
     &     PARF(60+KBARY)) KFLS=4  
      KFLL=0  
      IF(KFLS.EQ.2.AND.KFLD.GT.KFLE.AND.KFLE.GT.KFLF) THEN    
         IF(KFLDS.EQ.1.AND.KFLA.EQ.KFLD) KFLL=1    
         IF(KFLDS.EQ.1.AND.KFLA.NE.KFLD) KFLL=INT(0.25+RLU(0)) 
         IF(KFLDS.EQ.3.AND.KFLA.NE.KFLD) KFLL=INT(0.75+RLU(0)) 
      ENDIF   
      IF(KFLL.EQ.0) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+KFLS,KFL1)    
      IF(KFLL.EQ.1) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+KFLS,KFL1)    

      kfq1  = kfla
      kfq2  = kflb
      kfq3  = kflc
      kfbar = kf 

      RETURN    
      END
