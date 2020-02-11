program post_hassan
  implicit none
  include 'param.txt'
  include 'common.txt'
  integer counter_post, f
  !make the grid
  call mkgrid(0)


  OPEN(100, file='output.dat', status='new')
  !load the data
  do f=0,95
     call load_avg_post(0,counter_post, f)

     do k =1,kmax
       do i= 1,imax
         write(100,*) k*dz+(kmax*dz)*f,rp(i),upp(i,k)
       enddo
     enddo
  enddo


end program post_hassan
