      ! T - Temperature (input)
      ! xg - Gaussian points (output)
      ! wg - Gaussian quadrature weights (output)
      ! ne - number of terms (input)
 
      subroutine polyfermi(T,xg,wg,ne)   
      implicit real*8 (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      parameter (ntm=100000)
      complex*16 z,s,ephi,emphi,b,c
      complex*16 w
      parameter (nem=100)
      real*8 phi(ntm),x(ntm),wr(ntm),pl(ntm,3),sum(0:nem)
      real*8 al(ntm),bl(ntm)
      real*8 pe(0:nem,ntm),xe(ntm),we(ntm)
      integer ipvt(nem)
      complex*16 xg(ne),wg(ne)
      
      gamma=3.d0-sqrt(8.d0)
      if(ne.gt.nem) stop'polyfermi'
      n=2+int(0.5+0.1d0/T) 
      if(n.gt.ntm.or.ne.gt.n) stop'polyfermi'
      sum(0)=0.d0
      sum1=0.d0
      r=(1.d0+gamma)/(2.d0*n)
      q=(1.d0-gamma)/(2.d0*n)
      p=(1.d0-gamma*gamma)/(8.d0*n)
      a=r*r
      do i=1,n
        phi(i)=(i-0.5d0)*pi/n
        ephi=exp(dcmplx(0.d0,phi(i)))
        emphi=conjg(ephi)
        b=r*emphi+q*ephi
        c=emphi*emphi-1.d0
        s=sqrt(b*b-a*c)
        z=ephi*(-b+s)/(0.5d0*a)
        w=dcmplx(0.d0,-1.d0)*(1.d0+0.5d0*z*r)*(1.d0-z*q)/(1.d0-z*p)*T
        if(ne.eq.n) then
      	xg(i)=z*T
      	wg(i)=w
        else
        wr(i)=1.d0     ! changed
        x(i)=phi(i)    ! changed
      	pl(i,1)=1.d0
      	sum(0)=sum(0)+wr(i)
      	sum1=sum1+wr(i)*x(i)
        endif
      enddo ! i
      
      if(ne.lt.n) then
        sum(1)=0.d0
        sum2=0.d0
        do i=1,n
      	al(1)=sum1/sum(0)
      	bl(1)=0.d0
      	pl(i,2)=x(i)-al(1)
      	p2=pl(i,2)*pl(i,2)
      	sum(1)=sum(1)+wr(i)*p2
      	sum2=sum2+wr(i)*x(i)*p2
        enddo
      
        ilm1=1
        il=2
        do l=2,ne
          sum(l)=0.d0
          sum1=sum2
          sum2=0.d0
          ilm2=ilm1
          ilm1=il
          il=mod(il,3)+1
          do i=1,n
            al(l)=sum1/sum(l-1)
            bl(l)=sum(l-1)/sum(l-2)
            pl(i,il)=(x(i)-al(l))*pl(i,ilm1)-bl(l)*pl(i,ilm2)
            p2=pl(i,il)*pl(i,il)
            sum(l)=sum(l)+wr(i)*p2
            sum2=sum2+wr(i)*x(i)*p2
          enddo
        enddo ! l
      
        nz=0
        do i=1,n
        if(pl(i,il).eq.0.d0) then
          nz=nz+1
          pe(0,nz)=1.d0
          xe(nz)=x(i)
          pe(1,nz)=xe(nz)-al(1)
          do l=2,ne
            pe(l,nz)=(xe(nz)-al(l))*pe(l-1,nz)-bl(l)*pe(l-2,nz)
          enddo ! l
      	else if(i.gt.1.and.pl(i,il)*pl(i-1,il).lt.0.d0) then
          nz=nz+1
          pe(0,nz)=1.d0
          x1=x(i-1)
          x2=x(i)
          p1=pl(i-1,il)
          p2=pl(i,il)
      
201       continue
          xe(nz)=x1-p1*(x1-x2)/(p1-p2)
          pe(1,nz)=xe(nz)-al(1)
          do l=2,ne
            pe(l,nz)=(xe(nz)-al(l))*pe(l-1,nz)-bl(l)*pe(l-2,nz)
          enddo ! l
          if(abs(p1).lt.abs(p2)) then
            x2=xe(nz)
            p2=pe(ne,nz)
          else
            x1=xe(nz)
            p1=pe(ne,nz)
          endif
          if(abs(x1-x2).gt.1.d-14) goto 201
      	endif
        enddo  ! i
        we(nz)=0.d0
      
        we(1)=sum(0)
        call dgesv(ne,1,pe(0,1),nem+1,ipvt,we,nem,info)
      
        do nz=1,ne
          ephi=exp(dcmplx(0.d0,xe(nz)))   !changed
          emphi=conjg(ephi)
          b=r*emphi+q*ephi
          c=emphi*emphi-1.d0
          s=sqrt(b*b-a*c)
          z=ephi*(-b+s)/(0.5d0*a)
          w=dcmplx(0.d0,-1.d0)*(1.d0+0.5d0*z*r)*(1.d0-z*q)/(1.d0-z*p)
          xg(nz)=z*T
          wg(nz)=w*we(nz)*T    ! changed
        enddo ! nz
      endif ! ne
      return
      end !subroutine polyfermi
