program Matsu

      implicit real*8 (a-h,o-z)
      parameter (nem=100)
!      parameter (pi=3.141592653589793238d0)
      real*8, allocatable:: xg(:),wg(:)

      pi=4.d0*atan(1.d0)
      T=1.d0/256.d0
      eta=1.d-7
      ne=28
      allocate(xg(ne))
      allocate(wg(ne))
      ioffset=10
      offset=2.d0*pi*T*ioffset
      call polymatsu(T,xg,wg,ne,eta,offset)
      print *, 'pi*T=',pi*T
      print *, xg
      print *, '======'
      print *, wg
      sum1=0.d0
      sum2=0.d0
      do i=ne,1,-1
        sum1=sum1+wg(i)/(xg(i)*xg(i))
        sum2=sum2+wg(i)*exp(-xg(i)/T)
      enddo
      do i=ioffset,1,-1
        x=(2.d0*i-1.d0)*pi*T
        sum1=sum1+1.d0/(x*x)
        sum2=sum2+exp(-x/T)
      enddo
      sum1=sum1*(pi*T)**2
      print *, 'test 1:'
      print *, 'sum=',sum1,'  exact=',0.125d0*pi*pi
      print *, 'test 2:'
      print *, 'sum=',sum2,'  exact=',0.5d0/sinh(pi)
      end


      subroutine polymatsu(T,xg,wg,ne,eta,offset)
  ! xg-positions of Gauss points;
  ! wg-weight;
  ! ne-number of Gauss points
  ! eta - tolerance
  implicit real*8 (a-h,o-z)
  parameter (pi=3.141592653589793238d0)
  parameter (ntm=10000000)
integer ne
  Real*8 T ! temperature 
  parameter (nem=100)
  real*8 x(ntm),wr(ntm),pl(ntm,3),sum(0:nem)
  real*8 al(nem),bl(nem)
  real*8 xe(nem),we(nem)
  real*8 xg(ne),wg(ne)

real*8 ZZ(nem,nem),work2(nem)

      if(ne.gt.nem) stop 'polymatsu'
      n=1.d0/eta
      if(n.gt.ntm.or.ne.gt.n) stop 'polymatsu'
  do i=1,n
        omega=(2.d0*i)*pi*T+offset
        wr(i)=1.d0/(omega*omega)
        x(i)=1.d0/omega
        pl(i,1)=1.d0	!zeroth order polynomial
  enddo
  phi0=1.d0/(2.d0*pi*T)
  sum(0)=0.5d0*wr(n)+x(n)*phi0
  sum1=0.5d0*wr(n)*(x(n)+phi0)
  do i=n-1,1,-1
        sum(0)=sum(0)+wr(i)
        sum1=sum1+wr(i)*x(i)
  enddo			! i
 

  ! first order polynomial
     al(1)=sum1/sum(0)
     bl(1)=0.d0
     do i=1,n
        pl(i,2)=x(i)-al(1)
     enddo
     p2=pl(n,2)*pl(n,2)
     dp2=2.d0*pl(n,2)*phi0
     sum(1)=(0.5d0*wr(n)+x(n)*phi0)*p2+0.5d0*wr(n)*dp2
     sum2=0.5d0*wr(n)*((x(n)+phi0)*p2+x(n)*dp2/3.d0)
     do i=n-1,1,-1 ! reverse order to improve accuracy
        p2=pl(i,2)*pl(i,2)
        sum(1)=sum(1)+wr(i)*p2
        sum2=sum2+wr(i)*x(i)*p2
     enddo

     ! recursion formula of polynomials
     ilm1=1
     il=2
     do l=2,ne
        ilm2=ilm1
        ilm1=il
        il=mod(il,3)+1
        sum1=sum2
        al(l)=sum1/sum(l-1)
        bl(l)=sum(l-1)/sum(l-2)
        do i=1,n
           pl(i,il)=(x(i)-al(l))*pl(i,ilm1)-bl(l)*pl(i,ilm2)
        enddo
        p2=pl(n,il)*pl(n,il)
        dp2=(pl(n,il)+pl(n-1,il))*(pl(n,il)-pl(n-1,il))
     sum(l)=(0.5d0*wr(n)+x(n)*phi0)*p2+0.5d0*wr(n)*dp2
     sum2=0.5d0*wr(n)*((x(n)+phi0)*p2+x(n)*dp2/3.d0)
        do i=n-1,1,-1
           p2=pl(i,il)*pl(i,il)
           sum(l)=sum(l)+wr(i)*p2
           sum2=sum2+wr(i)*x(i)*p2
        enddo
	
     enddo		! l

!~ =======Golub-Welsch algorithm=====
do i=2,ne
bl(i)=sqrt(bl(i))
enddo
xe=al
call DSTEV("V", Ne, xe(1:ne), bl(2:ne), ZZ, nem, WORK2, INFO )

!~ ======Weight==========
do i=1,ne
we(i)=ZZ(1,i)*ZZ(1,i)*sum(0)
enddo
  
     do nz=1,ne
        xg(nz)=1.d0/xe(nz)-pi*T
        wg(nz)=we(nz)/(xe(nz)*xe(nz))
     enddo		! nz
      return
      end ! subroutine polymatsu

