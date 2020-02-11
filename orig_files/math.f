
      subroutine TABLE_Data(rank)
      include 'param.txt'
      include 'common.txt'
       open(27,file='co2h_table.dat')
         do i=1,nTab
            read(27,*)tempTab(i),rhoTab(i),muTab(i),lamTab(i),cpTab(i),enthTab(i)
         enddo
         close(27)
       do i=1,nTab
         lamocpTab(i) = lamTab(i)/cpTab(i)
       enddo
      call spline(enthTab, rhoTab,    nTab, rho2Tab)
      call spline(enthTab, muTab,     nTab, mu2Tab)
      call spline(enthTab, lamTab,    nTab, lam2Tab)
      call spline(enthTab, cpTab,     nTab, cp2Tab)
      call spline(enthTab, lamocpTab, nTab, lamocp2Tab)
      call spline(enthTab, tempTab,   nTab, temp2Tab)
      end
      subroutine state_max(enth,rho,mu,lam,Con,CP)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer tabkhi,tabklo,l
      real*8 ,dimension(imax,jmax,kmax)::enth,rho,mu,lam,Con,CP

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               tabkhi = 0
               tabklo = 0
               call splint(enthTab,rhoTab,rho2Tab,nTab,enth(i,j,k),rho(i,j,k),tabkhi,tabklo)
               call splint(enthTab,muTab,mu2Tab,nTab,enth(i,j,k),mu(i,j,k), tabkhi,tabklo)
               call splint(enthTab,lamocpTab,lamocp2Tab,nTab,enth(i,j,k),lam(i,j,k),tabkhi,tabklo)
               call splint(enthTab,cpTab,cp2Tab,nTab,enth(i,j,k),Cp(i,j,k), tabkhi,tabklo)
               call splint(enthTab,lamTab,lam2Tab,nTab,enth(i,j,k),Con(i,j,k),tabkhi,tabklo)
               mu(i,j,k)  = mu(i,j,k)/Re
               lam(i,j,k) = lam(i,j,k)/(Re*Pr)
            enddo
         enddo
      enddo
      return
      end

      subroutine state_1(enth,rho,mu,lam,Con,CP)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer tabkhi,tabklo,l
      real*8 ,dimension(0:i1,0:j1,0:k1)::enth,rho,mu,lam,Con,CP

      do k=0,k1
         do j=0,j1
            do i=0,i1
               tabkhi = 0
               tabklo = 0
               call splint(enthTab,rhoTab,rho2Tab,nTab,enth(i,j,k),rho(i,j,k),tabkhi,tabklo)
               call splint(enthTab,muTab,mu2Tab,nTab,enth(i,j,k),mu(i,j,k), tabkhi,tabklo)
               call splint(enthTab,lamocpTab,lamocp2Tab,nTab,enth(i,j,k),lam(i,j,k),tabkhi,tabklo)
               call splint(enthTab,cpTab,cp2Tab,nTab,enth(i,j,k),Cp(i,j,k), tabkhi,tabklo)
               call splint(enthTab,lamTab,lam2Tab,nTab,enth(i,j,k),Con(i,j,k),tabkhi,tabklo)
               mu(i,j,k)  = mu(i,j,k)/Re
               lam(i,j,k) = lam(i,j,k)/(Re*Pr)
            enddo
         enddo
      enddo
      return
      end
!

      subroutine spline(x, y, n, y2)
      implicit none
      integer   i, k, n, nmax
      parameter  (nmax=5000)
      real*8    yp1, ypn, x(n), y(n), y2(n), p, qn, sig, un, u(nmax)

      y2(1) = 0.
      u(1)  = 0.
      do i=2, n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     &        (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

      do k=n-1, 1, -1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
      end

      subroutine splint(xa,ya,y2a,n,x,y,khi,klo)
      implicit none
      integer n,k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n), a,b,h

!     if ((khi.eq.0) .and. (klo.eq.0)) then
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif

!     endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

      subroutine funcNewtonSolveRG(enth_i1, enth_imax)
      implicit none
      integer nIterNewton,success
      real*8 enth_i1, enth_imax, fxValue, fxValue1,lamOcpinter
        success = 1 
        fxValue = 1000.0
        nIterNewton = 0

        enth_i1= enth_imax! +(Rp(i1)-Rp(imax))*Qwall/lamOcpinter 

        do while (abs(fxValue).gt.1.0e-10)
          call funcNewtonBC(enth_i1,        enth_imax, fxValue)
          call funcNewtonBC(enth_i1+1.0e-8, enth_imax, fxValue1)
          enth_i1 = enth_i1 - fxValue/((fxValue1-fxValue)/1.0e-8)
          if (nIterNewton.gt.200) then
            fxValue = 0.0
            success = 0
          endif
          nIterNewton = nIterNewton + 1
        enddo
        if (success.eq.0) then
          write (*,*) 'newton didnt converge, enthimax= ',enth_imax,',', nIterNewton, ', ', enth_i1
          stop
        endif
      end

      subroutine funcNewtonBC(enth, enthIMAX, fxValue)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer tabkhi,tabklo
      real*8 enth,lamOcpinter,enthIMAX,fxValue
        tabkhi = 0.0
        tabklo = 0.0
        call splint(enthTab,lamocpTab,lamocp2Tab,nTab,0.5*(enth+enthIMAX),lamOcpinter,tabkhi,tabklo)
        fxValue = enth - enthIMAX - (Rp(i1)-Rp(imax))*Qwall/lamOcpinter
      end


      subroutine shiftb(UT,UP,rank)
      implicit none
      include 'param.txt'
!      include 'common.txt'
      include 'mpif.h'
      integer ileng,rankb,rankf,ierr
      integer itag,status(MPI_STATUS_SIZE),l
      real*8 ut(imax,0:k1),up(imax),UTMP(imax)
      parameter (ileng= imax )
      itag = 11
      do i=1,imax
        utmp(i) =UT(i,1)
      enddo
      rankf=rank+1
      rankb=rank-1
      if(rank.eq.px-1)rankf=MPI_PROC_NULL
      if(rank.eq.   0)rankb=MPI_PROC_NULL
         call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag,
     &                     up,ileng,MPI_REAL8,rankf,itag,
     &                     MPI_COMM_WORLD,status,ierr)
      if (rank.eq.   0) then
        call MPI_Send(utmp,ileng,MPI_REAL8,px-1,itag, MPI_COMM_WORLD,ierr)
        endif
       if (rank.eq.px-1) then
           call MPI_RECV(up,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,status,ierr)
          endif
      end

      subroutine shiftf(UT,UP,rank)
      implicit none
      include 'param.txt'
!      include 'common.txt'
      include 'mpif.h'
      integer ileng,rankb,rankf
      integer  itag,status(MPI_STATUS_SIZE),l,ierr
      real*8 ut(imax,0:k1),up(imax),UTMP(imax)
      parameter (ileng= imax )
      itag = 10
        do i=1,imax
        UTMP(i) =UT(i,kmax)
        enddo
      rankf=rank+1
      rankb=rank-1
      if(rank.eq.px-1)rankf=MPI_PROC_NULL
      if(rank.eq.   0)rankb=MPI_PROC_NULL
       call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag, up,ileng,MPI_REAL8,rankb,itag,
     &                     MPI_COMM_WORLD,status,ierr)
      if (rank.eq.px-1) then
      call MPI_SEND(UTMP,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,ierr)
      endif
      if (rank.eq.0) then
      call MPI_RECV(UP,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,status,ierr)
      endif
      end


      subroutine shiftb1(UT,UP,rank)
      implicit none
      include 'param.txt'
!      include 'common.txt'
      include 'mpif.h'
      integer ileng,rankb,rankf,ierr
      integer itag,status(MPI_STATUS_SIZE),l
      real*8 ut(0:k1),up,UTMP
      parameter (ileng= 1 )
      itag = 11
        utmp =UT(1)
      rankf=rank+1
      rankb=rank-1
      if(rank.eq.px-1)rankf=MPI_PROC_NULL
      if(rank.eq.   0)rankb=MPI_PROC_NULL
         call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag,
     &                     up,ileng,MPI_REAL8,rankf,itag,
     &                     MPI_COMM_WORLD,status,ierr)
      if (rank.eq.   0) then
        call MPI_Send(utmp,ileng,MPI_REAL8,px-1,itag, MPI_COMM_WORLD,ierr)
        endif
       if (rank.eq.px-1) then
           call MPI_RECV(up,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,status,ierr)
          endif
      end

      subroutine shiftf1(UT,UP,rank)
      implicit none
      include 'param.txt'
!      include 'common.txt'
      include 'mpif.h'
      integer ileng,rankb,rankf
      integer  itag,status(MPI_STATUS_SIZE),l,ierr
      real*8 ut(0:k1),up,UTMP
      parameter (ileng= 1 )
      itag = 10
        UTMP =UT(kmax)
      rankf=rank+1
      rankb=rank-1
      if(rank.eq.px-1)rankf=MPI_PROC_NULL
      if(rank.eq.   0)rankb=MPI_PROC_NULL
       call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag, up,ileng,MPI_REAL8,rankb,itag,
     &                     MPI_COMM_WORLD,status,ierr)
      if (rank.eq.px-1) then
      call MPI_SEND(UTMP,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,ierr)
      endif
      if (rank.eq.0) then
      call MPI_RECV(UP,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,status,ierr)
      endif
      end

