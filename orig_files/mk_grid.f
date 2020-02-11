       subroutine mkgrid1(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8  delta(imax),Yplus,X,rmax,drr,Rei,fact
      real*8  pii,y,y1,y2,fA,fB,fC

      pii   = 4.0*atan(1.0)
      Rei   = 1.0/Re
      dphi  = 2*pii/(jmax_tot)
      dz    = 1.0*LoD/(kmax_tot)


      fA = 0.12
      fB = 2.4
       ru(0)=0.

!       do i = 1,imax
!       
!
!       ru(i) = ru(i-1)-.75*real(i**2)/imax**2.+1.
!       enddo



!     do i = 1,imax
!      ru(i) = ru(i-1)+0.5/real(imax)
!     enddo
      do i = 1,imax
        fact = (i-0.)/(imax-0.)
        ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB))
 
        ru(i) = ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
        delta(i)=0.5*(ru(i)-ru(i-1))
       enddo



      do i=0,imax
        ru(i)=ru(i)/(2.*ru(imax))
      enddo

      do i = 1,imax
        Rp(i) = (Ru(i)+Ru(i-1))/2.0
        dr(i) = (Ru(i)-Ru(i-1))
      enddo
      dr(i1) = dr(imax)
      Ru(i1) = Ru(imax) + dr(i1)
      Rp(i1) = Ru(imax) + dr(i1)/2.0
      dr(0)  = dr(1)
      Rp(0)  = Ru(0) - dr(0)/2.0

      if (rank.eq.0) then
         open(11,file = 'grid.txt')
         write(11,*) Re, imax
         do i=1,imax
            Yplus = (0.5-Rp(i))*Re
            write(11,'(i5,4F12.6)') i,yplus,Ru(i),Rp(i),dr(i)
         enddo
         close(11)
      endif
      end



      subroutine mkgrid(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer iimax,ii,imax0
      real*8  delta(imax),Yplus,rmax,drr,Rei,a,b,x,y,c,pi
      pi    = 4.0*atan(1.0)
      Rei   = 1.0/Re
      dphi  = 2*pi/(jmax_tot)
      dz    = 1.0*LoD/(kmax_tot)
      iimax=0!imax/8
      imax0=imax-iimax
      b=1.015
      a=(b+1.)/(b-1.)
      ru(0)=0.
      c=1.75
      DO i = 1,imax0
      x=1.*i/(imax0)
      y=1.*(i-1)/imax0
      delta(i)=0.24+1.2*(b-1.)*(((a**(c*(1.-y))-1.)/(a**(c*(-y))+1.))-((a**(c*(1.-x))-1.)/(a**(c*(-x))+1.)))
      ENDDO

       DO i =1,imax
       ru(i)=ru(i-1)+delta(i)
       ENDDO

      do i=0,imax_tot
        ru(i)=ru(i)/(2.*ru(imax_tot))
      enddo

      do i = 1,imax
        Rp(i) = (Ru(i)+Ru(i-1))/2.0
        dr(i) = (Ru(i)-Ru(i-1))
      enddo
      dr(i1) = dr(imax)
      Ru(i1) = Ru(imax) + dr(i1)
      Rp(i1) = Ru(imax) + dr(i1)/2.0
      dr(0)  = dr(1)
      Rp(0)  = Ru(0) - dr(0)/2.0

      if (rank.eq.0) then
         open(11,file = 'grid.txt')
         write(11,*) Re, imax
         do i=1,imax
            Yplus = (0.5-Rp(i))*Re
            write(11,'(i5,4F12.6)') i,yplus,Ru(i),Rp(i),dr(i)
         enddo
         close(11)
      endif
      end



