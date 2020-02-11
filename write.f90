subroutine load_avg(rank)
  implicit none
  include 'param.txt'
  include 'common.txt'
  include 'mpif.h'
  character*5 cha
  write(cha,'(I5.5)')rank
  open (100,file='Load_Avg/save_avg.'//cha,form='unformatted')
  read (100)Pnnff,Velnnff,RVelnnff,Taunnff,Rnnff,ekmnnff,ekhnnff,Knnff,Cpnnff,Cnnff,RCnnff,Qnn,Qff  
  close(100)
  return
end
subroutine load_avg_post(ini,counter_post,rank)
  implicit none
  include 'param.txt'
  include 'common.txt'
  include 'mpif.h'
  integer counter_post,ini
  character*5 cha
  write(cha,'(I5.5)')rank
  open (100,file='Save_Avg/save_avg.'//cha,form='unformatted',  action='read' )
  read(100) counter_post
  read(100) PRO,TD,PS,PD,VD,DS,CTm,PROT,TDT,VDT,DST,CTT,PROVT,TDVT,PSVT,PDVT,VDVT,DSVT,CTVT,KIN,RU2pp,RV2pp,RW2pp, &
   RC2pp,RU3pp,RV3pp,RW3pp,RC3pp,RU4pp,RV4pp,RW4pp, &
   RC4pp,upp,vpp,wpp,cpp,P_rms,R_rms,ekh_rms,ekm_rms,Cp_rms,RCUpp,RCWpp,RUWpp,&
   pdfuw,pdfuc,pdfwc,pdfru,n_uw,n_uc,n_wc,n_ru,Quw,Quh,pdfu,pdfw,pdfc
  close(100)
end subroutine
