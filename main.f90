program post_hassan
  implicit none
  include 'param.txt'
  include 'common.txt'
  integer counter_post, f
  !make the grid
  call mkgrid(0)

   

  OPEN(101, file= 'output.dat', status='replace')
  !load the data
  do f=0,95
     call load_avg_post(0,counter_post, f)
     do k =1,kmax
       do i= 1,imax
         write(101,"(3E20.10)") k*dz+kmax*f*dz,rp(i),upp(i,k),vpp(i,k), wpp(i,k),RUWpp(i,k)
       enddo
     enddo
  enddo
  close(101)


end program post_hassan
