      subroutine Write_MOM(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr,tabkhi,tabklo,l
      real*8 DivU(imax,kmax),Norm(imax,kmax),cons,RR,ekmm,lamocp

       character*4 cha
      integer ss1,ssmax
      real*8 Ub1,dUbdx(kmax),Ubb(0:k1),BB(1:8)

       do j=1,6



       if (rank.eq.pos(1,j))then
       write(cha,'(I4.4)')pos(1,j)
       k=pos(2,j)
       open (100,file='data/KIN.'//cha)
       do i=1,imax
        write(100,'(7E18.11)')(0.5-rp(i))*2.,KIN(i,k),KIN(i,k)/Rnnff(imax,k,1),Ruwpp(i,k),RU2pp(i,k),RV2pp(i,k),RW2pp(i,k)
       enddo
       close(100)

       open (100,file='data/RMS.'//cha)
       do i=1,imax
        write(100,'(10E18.11)')(0.5-rp(i))*2.,sqrt(Ru2pp(i,k)/Rnnff(i,k,2)),sqrt(Rv2pp(i,k)/Rnnff(i,k,3))
     1                                      ,sqrt(Rw2pp(i,k)*Rnnff(i,k,4)),sqrt(R_rms(i,k)),sqrt(ekm_rms(i,k)),sqrt(ekh_rms(i,k))
     1                                      ,sqrt(Cp_rms(i,k)),sqrt(P_rms(i,k))

       enddo
       close(100)

       open (12,file='data/PD.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,PD(i,k,1:3)
        enddo
       close(12)

       open (12,file='data/PS.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,PS(i,k,1:3)
        enddo
       close(12)


       open (12,file='data/TD.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,TD(i,k,1:3)
        enddo
       close(12)
       open (12,file='data/DS.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,DS(i,k,1:3)
        enddo
       close(12)

       open (12,file='data/VD.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,VD(i,k,1:3)
        enddo
       close(12)

       open (12,file='data/PRO.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,PRO(i,k,1:3)
        enddo
       close(12)

       open (12,file='data/Budgets_M.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,sum(TD(i,k,1:3))/2.,sum(PRO(i,k,1:3))/2.,
     1                                      sum(PD(i,k,1:3))/2.,sum(PS(i,k,1:3))/2.,
     1                                      sum(DS(i,k,1:3))/2.,sum(VD(i,k,1:3))/2.,
     1                                      sum(EntU1(i,k,1:3))/2.,sum(EntU2(i,k,1:3))/2.,
     1                                      sum(BY(i,k,1:3))/2.
        enddo
       close(12)
       open (12,file='data/Budgets_test.'//cha)
       do i=2,imax

       BB(1)=(TD(i,k,1)+TD(i-1,k,1)+TD(i,k,2)+TD(i,k,2)+TD(i,k,3)+TD(i,k-1,3))/4.
       BB(2)=(PRO(i,k,1)+PRO(i-1,k,1)+PRO(i,k,2)+PRO(i,k,2)+PRO(i,k,3)+PRO(i,k-1,3))/4.
       BB(3)=(PD(i,k,1)+PD(i-1,k,1)+PD(i,k,2)+PD(i,k,2)+PD(i,k,3)+PD(i,k-1,3))/4.
       BB(4)=(PS(i,k,1)+PS(i-1,k,1)+PS(i,k,2)+PS(i,k,2)+PS(i,k,3)+PS(i,k-1,3))/4.
       BB(5)=(DS(i,k,1)+DS(i-1,k,1)+DS(i,k,2)+DS(i,k,2)+DS(i,k,3)+DS(i,k-1,3))/4.
       BB(6)=(VD(i,k,1)+VD(i-1,k,1)+VD(i,k,2)+VD(i,k,2)+VD(i,k,3)+VD(i,k-1,3))/4.
       BB(7)=(EntU1(i,k,1)+EntU1(i-1,k,1)+EntU1(i,k,2)+EntU1(i,k,2)+EntU1(i,k,3)+EntU1(i,k-1,3))/4.
       BB(8)=(EntU2(i,k,1)+EntU2(i-1,k,1)+EntU2(i,k,2)+EntU2(i,k,2)+EntU2(i,k,3)+EntU2(i,k-1,3))/4.
       write(12,'(16E18.10)')(0.5-rp(i))*2.,BB(1:8)
        enddo
       close(12)

       open (12,file='data/Addi_M.'//cha)
       do i=2,imax
       write(12,'(16E18.10)')(0.5-rp(i))*2.,EntU1(i,k,1:3),EntU2(i,k,1:3)
        enddo
       close(12)


        open (10,file='data/Xpp.'//cha)
        do i=1,imax
        write(10,'(10E18.10)')(0.5-rp(i))*2.,Upp(i,k),Wpp(i,k),Cpp(i,k)
        enddo
        close(10)

        open (10,file='data/Utz.'//cha)
        do i=1,imax
        write(10,'(12E18.9)')(0.5-rp(i))*2.,Vel_t(i,k,4,3),Velnnff(i,k,4,3)
       enddo
        close(10)


        open (12,file='data/Budgets_Uh.'//cha)
        do i=2,imax
        write(12,'(15E18.10)')(0.5-rp(i))*2.,-TDVT(i,k,1),-PROVT(i,k,1)
     1                                      ,-PDVT(i,k,1)-PSVT(i,k,1)
     1                                      ,-DSVT(i,k,1),-VDVT(i,k,1)
     1                                      ,-EntCU1(i,k,1),-ENTCU2(i,k,1)
        enddo
       close(12)

!       open (12,file='data/Budgets_Uh.'//cha)
!      do i=2,imax
!      BB(1)=(TDVT(i,k,1)+TDVT(i-1,k,1))/2.
!      BB(2)=(PROVT(i,k,1)+PROVT(i-1,k,1))/2.
!      BB(3)=(PDVT(i,k,1)+PDVT(i-1,k,1))/2.
!      BB(4)=(PSVT(i,k,1)+PSVT(i-1,k,1))/2.
!      BB(5)=(DSVT(i,k,1)+DSVT(i-1,k,1))/2.
!      BB(6)=(VDVT(i,k,1)+VDVT(i-1,k,1))/2.
!      BB(7)=(EntCU1(i,k,1)+EntCU1(i-1,k,1))/2.
!      BB(8)=(EntCU2(i,k,1)+EntCU2(i-1,k,1))/2.
!      write(12,'(16E18.10)')(0.5-rp(i))*2.,BB(1:8)
!       enddo
!      close(12)


        open (100,file='data/SKE_FLA.'//cha)
        do i=2,imax
        write(100,'(6E18.9)')(0.5-rp(i))*2.,Ru3pp(i,k)/Ru2pp(i,k)**(1.5),
     1                                      RW3pp(i,k)/RW2pp(i,k)**(1.5),
     1                                      Ru4pp(i,k)/Ru2pp(i,k)**2.,
     1                                      Rw4pp(i,k)/Rw2pp(i,k)**2.
        enddo
        close(100)

       endif

        enddo                  

        end

      subroutine Write_Energy(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
       character*4 cha
       real*8 qrri,qrro,cons,ff,PRB,REB,NUD,n
       integer istap,jstart,ss1,ssmax
       do j=1,6
       if (rank.eq.pos(1,j))then
       write(cha,'(I4.4)')pos(1,j)
       k=pos(2,j)

       open (12,file='data/heat_con.'//cha)
        do i=2,imax
            write(12,'(13E18.10)')(1.-2*ru(i)),Qff(i,k,1),-RCUpp(i,k)
        enddo
        close(12)

        open (12,file='data/Hr.'//cha)
        do i=2,imax
           write(12,'(13E18.10)')(0.5-ru(i))*2.,Cnnff(i,k,2)/Cnnff(imax,k,2),C_t(i,k,2)/C_t(imax,k,2)
        enddo
       close(12)



        open (12,file='data/Budgets_H.'//cha)
        do i=2,imax
        write(12,'(16E18.10)')(0.5-rp(i))*2.,TDT(i,k),PROT(i,k),CTT(i,k)
     1                                      ,DST(i,k) ,VDT(i,k),EntC(i,k)
        enddo
        close(12)
       endif
       enddo




       end



      subroutine load_3D(counter_post,rank)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr,counter_post
      character*5 cha
      write(cha,'(I5.5)')counter_post
      call decomp_2d_read_one(1,unew,'../DATA/U.'//cha)
      call decomp_2d_read_one(1,vnew,'../DATA/V.'//cha)
      call decomp_2d_read_one(1,wnew,'../DATA/W.'//cha)
      call decomp_2d_read_one(1,Cnew,'../DATA/C.'//cha)
      call decomp_2d_read_one(1,P_2d,'../DATA/P.'//cha)
      return
      end
      subroutine save_avg(ini,counter_post,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer counter_post,ini
      character*5 cha
      write(cha,'(I5.5)')rank
      if (ini.eq.1)then
      open (100,file='Save_Avg/save_avg.'//cha,form='unformatted')

      write(100) counter_post
      write(100) PRO,TD,PS,PD,VD,DS,CTm,PROT,TDT,VDT,DST,CTT,PROVT,TDVT,PSVT,PDVT,VDVT,DSVT,CTVT,KIN,RU2pp,RV2pp,RW2pp,
     1 RC2pp,RU3pp,RV3pp,RW3pp,RC3pp,RU4pp,RV4pp,RW4pp,
     1 RC4pp,upp,vpp,wpp,cpp,P_rms,R_rms,ekh_rms,ekm_rms,Cp_rms,RCUpp,RCWpp,RUWpp,
     1 pdfuw,pdfuc,pdfwc,pdfru,n_uw,n_uc,n_wc,n_ru,Quw,Quh,pdfu,pdfw,pdfc
      close(100)
      endif
      if (ini.eq.0)then
      open (100,file='Save_Avg/save_avg.'//cha,form='unformatted')
      read(100) counter_post
      read(100) PRO,TD,PS,PD,VD,DS,CTm,PROT,TDT,VDT,DST,CTT,PROVT,TDVT,PSVT,PDVT,VDVT,DSVT,CTVT,KIN,RU2pp,RV2pp,RW2pp,
     1 RC2pp,RU3pp,RV3pp,RW3pp,RC3pp,RU4pp,RV4pp,RW4pp,
     1 RC4pp,upp,vpp,wpp,cpp,P_rms,R_rms,ekh_rms,ekm_rms,Cp_rms,RCUpp,RCWpp,RUWpp,
     1 pdfuw,pdfuc,pdfwc,pdfru,n_uw,n_uc,n_wc,n_ru,Quw,Quh,pdfu,pdfw,pdfc
      close(100)
      endif
      return
      end
      subroutine load_avg(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      write(cha,'(I5.5)')rank
      open (100,file='Load_Avg/save_avg.'//cha,form='unformatted')
      read(100)Pnnff,Velnnff,RVelnnff,Taunnff,Rnnff,ekmnnff,ekhnnff,Knnff,Cpnnff,Cnnff,RCnnff,Qnn,Qff  
      close(100)
      return
      end

      subroutine Q_P(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
       character*4 cha
       real*8 qrri,qrro,cons,ff,PRB,REB,NUD,n,n1,n2
       integer istap,jstart,ss1,ssmax
       do j=1,6
       if (rank.eq.pos(1,j))then
       write(cha,'(I4.4)')pos(1,j)
       k=pos(2,j)

        open (10,file='data/quad.'//cha)
        do i=1,imax
        write(10,'(9E18.9)')(0.5-Rp(i))*2 ,Quw(i,k,1),Quw(i,k,2),
     &                                     Quw(i,k,3),Quw(i,k,4),
     &                                     Quh(i,k,1)*Re*Pr/Qwall,Quh(i,k,2)*Re*Pr/Qwall,
     &                                     Quh(i,k,3)*Re*Pr/Qwall,Quh(i,k,4)*Re*Pr/Qwall
        enddo
        close(10)



        open (15,file='data/PDFu.'//cha)
        do n=-conter1,conter1
         write(15,'(6E18.9)')1.*n*bu1,pdfu(imax,k,n)/(sum(pdfu(imax,k,-conter1:conter1))*bu1)
     1   ,pdfu(imax-3,k,n)/(sum(pdfu(imax-3,k,-conter1:conter1)) *bu1)
     1   ,pdfu(imax-4,k,n)/(sum(pdfu(imax-4,k,-conter1:conter1)) *bu1)
     1   ,pdfu(imax-10,k,n)/(sum(pdfu(imax-10,k,-conter1:conter1))*bu1)
     1   ,pdfu(imax-20,k,n)/(sum(pdfu(imax-20,k,-conter1:conter1))*bu1)
         enddo
         close(15)

        open (15,file='data/PDFw.'//cha)
        do n=-conter1,conter1
         write(15,'(6E18.9)')1.*n*bw1,pdfw(imax,k,n)/(sum(pdfw(imax,k,-conter1:conter1))*bw1)
     1   ,pdfw(imax-3,k,n)/(sum(pdfw(imax-3,k,-conter1:conter1)) *bw1)
     1   ,pdfw(imax-4,k,n)/(sum(pdfw(imax-4,k,-conter1:conter1)) *bw1)
     1   ,pdfw(imax-10,k,n)/(sum(pdfw(imax-10,k,-conter1:conter1))*bw1)
     1   ,pdfw(imax-20,k,n)/(sum(pdfw(imax-20,k,-conter1:conter1))*bw1)
         enddo
         close(15)

        open (15,file='data/PDFc.'//cha)
        do n=-conter1,conter1
         write(15,'(6E18.9)')1.*n*bc1,pdfc(imax,k,n)/(sum(pdfc(imax,k,-conter1:conter1))*bc1)
     1    ,pdfc(imax-3,k,n)/(sum(pdfc(imax-3,k,-conter1:conter1)) *bc1)
     1    ,pdfc(imax-4,k,n)/(sum(pdfc(imax-4,k,-conter1:conter1)) *bc1)
     1    ,pdfc(imax-10,k,n)/(sum(pdfc(imax-10,k,-conter1:conter1))*bc1)
     1    ,pdfc(imax-20,k,n)/(sum(pdfc(imax-20,k,-conter1:conter1))*bc1)
         enddo
         close(15)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


      open(15,file='2D/PDFuw.'//cha)
       write(15,*)'VARIABLES="X","Y","uw1","uw2","uw3","uw4"'
       write(15,*) 'ZONE I=  ', 2*conter2+1,' J=  ',2*conter2+1,' F=POINT '
      do n1=-conter2,conter2
      do n2=-conter2,conter2
       write(15,'(6E18.9)')real(n1*bu2),real(n2*bw2),pdfuw(imax-2,n1,n2),pdfuw(imax-4,n1,n2)
     1                     ,pdfuw(imax-8,n1,n2),pdfuw(imax-12,n1,n2)!/(sum(pdfu(i,k,-conter2:conter2))*bu2)+1e-10)
         enddo
      enddo
      close(15)
      open(15,file='2D/PDFuc.'//cha)
       write(15,*)'VARIABLES="X","Y","uc1","uc2","uc3","uc4"'
       write(15,*) 'ZONE I=  ', 2*conter2+1,' J=  ',2*conter2+1,' F=POINT '
      do n1=-conter2,conter2
      do n2=-conter2,conter2
       write(15,'(6E18.9)')real(n1*bu2),real(n2*bc2),pdfuc(imax-2,n1,n2),pdfuc(imax-4,n1,n2)
     1     ,pdfuc(imax-8,n1,n2),pdfuc(imax-12,n1,n2)!/(sum(pdfu(i,k,-conter2:conter2))*bu2)+1e-10)
         enddo
      enddo
      close(15)
      open(15,file='2D/PDFwc.'//cha)
       write(15,*)'VARIABLES="X","Y","wc1","wc2","wc3","wc4"'
       write(15,*) 'ZONE I=  ', 2*conter2+1,' J=  ',2*conter2+1,' F=POINT '
      do n1=-conter2,conter2
      do n2=-conter2,conter2
       write(15,'(6E18.9)')real(n1*bw2),real(n2*bc2),pdfwc(imax-2,n1,n2),pdfwc(imax-4,n1,n2)
     1    ,pdfwc(imax-8,n1,n2),pdfwc(imax-12,n1,n2)!/(sum(pdfu(i,k,-conter2:conter2))*bu2)+1e-10)
         enddo
      enddo
      close(15)
      open(15,file='2D/PDFru.'//cha)
       write(15,*)'VARIABLES="X","Y","ru1","ru2","ru3","ru4"'
       write(15,*) 'ZONE I=  ', 2*conter2+1,' J=  ',2*conter2+1,' F=POINT '
      do n1=-conter2,conter2
      do n2=-conter2,conter2
       write(15,'(6E18.9)')real(n1*br2),real(n2*bu2),pdfru(imax-2,n1,n2)*n1*n2,pdfru(imax-4,n1,n2)*n1*n2
     1   ,pdfru(imax-8,n1,n2)*n1*n2,pdfru(imax-12,n1,n2)*n1*n2!/(sum(pdfu(i,k,-conter2:conter2))*bu2)+1e-10)
         enddo
      enddo
      close(15)
      endif
      enddo
       end

       subroutine pack_cbuffer_name(chars,buffer,buffer_count,reset)
       implicit none
       integer, dimension(*), intent(inout) :: buffer
       character(len=*), intent(in) :: chars
       integer, intent(inout) :: buffer_count
       integer :: i, n, reset
         if (reset.eq.0) buffer_count = 0
         n = len_trim(chars)
         do i = 1,n
           buffer(buffer_count+i) = ICHAR(chars(i:i))
         end do
         buffer(buffer_count+n+1) = 0
         buffer_count = buffer_count+n+1
       end 

       subroutine output3d(U1,V1,W1,C1,rank)
       implicit none
       include 'param.txt'
       include 'common.txt'
       include 'mpif.h'
       integer*4 buffer(1000), tec_version, tec_byte,tec_prec,tec_nvalues_3D,tec_nvalues_2D,tec_format,tec_color,tec_repeat
       real*4 tec_marker1, tec_marker2
       integer bufferCount
       real*8 U1(0:i1,0:j1,0:k1),V1(0:i1,0:j1,0:k1),W1(0:i1,0:j1,0:k1),C1(0:i1,0:j1,0:k1)
       integer istap
       real*8 theta
       real*8 Unode,Vnode,Wnode,Cnode,Ux,Uy
       character*5 cha
       tec_byte = 1             ! byte order integer value
       tec_marker1 = 299.0
       tec_marker2 = 357.0
       tec_format = 1           ! 1..structured data
       tec_color = -1
       tec_repeat = 0
       tec_prec = 2             ! 1..single, 2..double precission
       tec_nvalues_3D = 4          ! number of variables to be written
       if (rank<48)then
       write(cha,'(I5.5)')rank
       open(45,file='3D/zf.'//cha,form='unformatted',access='stream')

       if (rank.eq.0) then
         write(45) '#!TDV75 '
         write(45) tec_byte
         call pack_cbuffer_name('DNS pipe',buffer,bufferCount, 0)
         write(45) buffer(1:bufferCount)
         write(45) tec_nvalues_3D
           call pack_cbuffer_name('x',buffer,bufferCount, 0)
           call pack_cbuffer_name('y',buffer,bufferCount, 1)
           call pack_cbuffer_name('z',buffer,bufferCount, 1)
           call pack_cbuffer_name('Vel x',buffer,bufferCount, 1)
!           call pack_cbuffer_name('Vel y',buffer,bufferCount, 1)
!           call pack_cbuffer_name('Vel z',buffer,bufferCount, 1)
!           call pack_cbuffer_name('C',buffer,bufferCount, 1)
           write(45) buffer(1:bufferCount)

         ! begin tecplot binary header section
         write(45) tec_marker1
           call pack_cbuffer_name('flowfield',buffer,bufferCount, 0)
           write(45) buffer(1:bufferCount)
           write(45) tec_format
         write(45) tec_color
         write(45) imax-41, jmax/2, kmax_tot/2
         write(45) tec_marker2
!        end tecplot binary header section
!        begin tecplot data
         write(45) tec_marker1
           write(45) tec_repeat
           do i=1,tec_nvalues_3D
             write(45) tec_prec   ! precision of the data
           enddo
       endif
         do k=1,kmax
         do j=1,jmax/2
           theta = j*dphi
           do i=42,imax

             Unode = 0.5*(U1(i,j,k)+U1(i-1,j,k))/Utau(i,k)
             Vnode = 0.5*(V1(i,j,k)+V1(i,j-1,k))/Utau(i,k)
!             Wnode = 0.5*(W1(i,j,k)+W1(i,j,k-1))/Utau(i,k)
             Wnode = (ru(i)*(V1(i+1,j,k)+V1(i,j,k)+V1(i+1,j-1,k)+V1(i,j-1,k))-
     1               ru(i-1)*(V1(i,j,k)+V1(i-1,j,k)+V1(i,j-1,k)+V1(i-1,j-1,k)))/(4.*rp(i)*dr(i))-
     1               (U1(i,j+1,k)+U1(i-1,j+1,k)-U1(i,j-1,k)-V1(i-1,j-1,k))/(4.*rp(i)*dphi)



             Ux = Unode*cos(theta) - Vnode*sin(theta)
             Uy = Unode*sin(theta) + Vnode*cos(theta)

                Cnode = C1(i,j,k)-Cnnff(i,k,1)
                write(45)rp(i)*cos(theta),rp(i)*sin(theta),(k+rank*kmax)*dz,Wnode!,Cnode
             enddo
          enddo
       enddo
       close(45)
       endif

       if (rank>47)then
       write(cha,'(I5.5)')rank-48
       open(45,file='3D/ft.'//cha,form='unformatted',access='stream')

       if (rank-48.eq.0) then
         write(45) '#!TDV75 '
         write(45) tec_byte
         call pack_cbuffer_name('DNS pipe',buffer,bufferCount, 0)
         write(45) buffer(1:bufferCount)
         write(45) tec_nvalues_3D
           call pack_cbuffer_name('x',buffer,bufferCount, 0)
           call pack_cbuffer_name('y',buffer,bufferCount, 1)
           call pack_cbuffer_name('z',buffer,bufferCount, 1)
           call pack_cbuffer_name('Vel x',buffer,bufferCount, 1)
!           call pack_cbuffer_name('Vel y',buffer,bufferCount, 1)
!           call pack_cbuffer_name('Vel z',buffer,bufferCount, 1)
!           call pack_cbuffer_name('C',buffer,bufferCount, 1)
           write(45) buffer(1:bufferCount)

         ! begin tecplot binary header section
         write(45) tec_marker1
           call pack_cbuffer_name('flowfield',buffer,bufferCount, 0)
           write(45) buffer(1:bufferCount)
           write(45) tec_format
         write(45) tec_color
         write(45) imax-41, jmax/2, kmax_tot/2
         write(45) tec_marker2
!        end tecplot binary header section
!        begin tecplot data
         write(45) tec_marker1
           write(45) tec_repeat
           do i=1,tec_nvalues_3D
             write(45) tec_prec   ! precision of the data
           enddo
       endif
         do k=1,kmax
         do j=1,jmax/2
           theta = j*dphi
           do i=42,imax

             Unode = 0.5*(U1(i,j,k)+U1(i-1,j,k))/Utau(i,k)
             Vnode = 0.5*(V1(i,j,k)+V1(i,j-1,k))/Utau(i,k)
!             Wnode = 0.5*(W1(i,j,k)+W1(i,j,k-1))/Utau(i,k)
             Wnode = (ru(i)*(V1(i+1,j,k)+V1(i,j,k)+V1(i+1,j-1,k)+V1(i,j-1,k))-
     1               ru(i-1)*(V1(i,j,k)+V1(i-1,j,k)+V1(i,j-1,k)+V1(i-1,j-1,k)))/(4.*rp(i)*dr(i))-
     1               (U1(i,j+1,k)+U1(i-1,j+1,k)-U1(i,j-1,k)-V1(i-1,j-1,k))/(4.*rp(i)*dphi)



             Ux = Unode*cos(theta) - Vnode*sin(theta)
             Uy = Unode*sin(theta) + Vnode*cos(theta)

                Cnode = C1(i,j,k)-Cnnff(i,k,1)
                write(45)rp(i)*cos(theta),rp(i)*sin(theta),(k+(rank-48)*kmax)*dz,Wnode!,Cnode
             enddo
          enddo
       enddo
       close(45)
       endif
        end


