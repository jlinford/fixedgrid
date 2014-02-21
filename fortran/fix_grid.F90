#define RUN_CASE 291+5
!#define REFINE_ALL 1
#define REFINE_NONE 1
program  fix_grid

  use global_vars

  implicit none

  ! local amr variables

  integer           :: NX,NY,i_source,j_source
  integer           :: loop_count, i, j!, k, kk
  integer           :: minstp, maxstp, istep, L
  integer           :: shift

  real              :: dt,dx,dy,FinalTime,time
  character*4	    :: chProc,chCase
  character*6	    :: chIter
  character*256     :: chTemp

#ifndef REFINE_ALL
#ifndef REFINE_NONE
#error ERROR
#endif
#endif

#ifdef REFINE_ALL
  real,dimension(320,96) :: conc
  real,dimension(320) :: conc_x
  real,dimension(96) :: conc_y
#endif
  
#ifdef REFINE_NONE
  real,dimension(40,12) :: conc 
  real,dimension(40) :: conc_x
  real,dimension(12) :: conc_y
#endif

  real              :: domain_size_x,domain_size_y,cell_size
  integer::n_var
  real::   cell_x,cell_y,bsize_x,bsize_y

  write(*,*) 'Running on ',0,' processors'
  write(chProc,"(I4)") 1000+0
  write(chCase,"(I4)") 1000+RUN_CASE


  ! ################################
  u=5.0    ! [m/s]
  v=0.0    ! [m/s]
  kh=100.0 ! [m^2/s]

  ! ################################

  i_MaxSpecies = N_VAR 
  conc(:,:)=O3_intial;

  minstp = 1
  maxstp = 500 ![Final time = maxstp*dt]
  time=0.0
  !! 12.0
  FinalTime=12.0*3600.0
  dt=50 ![sec]
  maxstp=INT((FinalTime-time)/dt)
  !maxstp=100
  ! ################################

  ! 4x4
  domain_size_x=400.0*1E+03 ![M]
  domain_size_y=120.0*1E+03 ![M]


  ! ################################
  EMMFlag=1
  ! 4x4
  point_source_x=62.5*1E+03 ![M]
  point_source_y=62.5*1E+03 ![M]


  ! ################################
#ifdef REFINE_ALL
  dx = 1250.0
  dy = 1250.0
  NX = 320
  NY = 96
  cell_size=1.25*1E+03      ![M]
  i_source=50
  j_source=50
#endif

#ifdef REFINE_NONE
  dx = 10000.0
  dy = 10000.0
  NX = 40
  NY = 12
  cell_size=10.0*1E+03      ![M]
  i_source=7
  j_source=7
#endif

  print*,'Running for ',(FinalTime-time)/3600.0,' [hours]'
  print*,'Steps=',maxstp-minstp+1,', each with ',dt,' [sec]'
  
  open(25,file='OUT_initial_solution'//&
       '_'//chCase(2:4)//'.'//&
       chProc(2:4),POSITION='REWIND')
  do j=1,NY
     do i=1,NX
        cell_x=i*dx-dx*0.5
        cell_y=j*dy-dy*0.5
        write(25,2020) &
             cell_x,cell_y,conc(i,j)
     enddo
  enddo
  close(25)
  
  open(25,file='OUT_initial_grid'//'_'//chCase(2:4)//'.'&
       //chProc(2:4),POSITION='REWIND')
  cell_x=domain_size_x/2
  cell_y=domain_size_y/2
  bsize_x=domain_size_x
  bsize_y=domain_size_y
  write(25,4040) &
       0,1,cell_x,cell_y,bsize_x,bsize_y 
  close(25)


  do istep = minstp, maxstp

     if(time .lt. 6*3600.0) then
        EMMFlag=1
     else
        EMMFlag=0
     endif

    
 
 

     do j=1,NY
         call DISCRETIZE(NX, reshape(conc(1:NX,j),(/NX/)), &
              u, kh, dx, dt/2, conc_x)
        conc(1:NX,j)=conc_x
     enddo
     
 
     do i=1,NX
        call DISCRETIZE(NY, reshape(conc(i,1:NY),(/NY/)), &
             v, kh, dy, dt, conc_y)
       conc(i,1:NY)=conc_y
     enddo


     do j=1,NY
        call DISCRETIZE(NX, reshape(conc(1:NX,j),(/NX/)), &
             u, kh, dx, dt/2, conc_x)
        conc(1:NX,j)=conc_x
     enddo

     if(EMMFlag.eq.1) then
        conc(i_source,j_source)=conc(i_source,j_source)+&
             dt*(EmmSource)/(dx*dy*1E+03)
!        print*,dt*(EmmSource)/(dx*dy*1E+03)
!        STOP
     endif
!     print*,'conc',istep,'---',conc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$     ! Save solution
     write(chIter,"(I6)") 100000+istep   
     open(25,file='OUT_grid_'//chCase(2:4)//'_'//chIter(2:6)//'.'//chProc(2:4),POSITION='REWIND')
     cell_x=domain_size_x/2
     cell_y=domain_size_y/2
     bsize_x=domain_size_x
     bsize_y=domain_size_y
     write(25,4040) &
          0,1,cell_x,cell_y,bsize_x,bsize_y 
     close(25)

    
     open(25,file='OUT_solution_'//chCase(2:4)//'_'//chIter(2:6)//'.'//chProc(2:4),POSITION='REWIND')
     
     do j=1,NY
        do i=1,NX
           cell_x=i*dx-dx*0.5
           cell_y=j*dy-dy*0.5
           write(25,2020) &
                cell_x,cell_y,conc(i,j)
        enddo
     enddo
     close(25)

     time=time+dt

     if(MOD(istep,10).eq. 0) then
        write(*,"(A,I5,A,I5,A,F9.2,A,F7.4,A)") &
             'It ',istep,' out of ',maxstp,'; Time=',&
             time,'[sec] (',time/3600.0,'[hours])'

     endif
     
  end do

  write(*,"(A,F9.2,A,F7.4,A)") 'Final time=',&
       time,'[sec]; ',time/3600.0,'[hours]'


  open(25,file='OUT_config_'//chCase(2:4)//'.txt',POSITION='REWIND')
  write(25,"(I4,2X,I4,2X,I4,2X,I4,2X,I4,2X,I4,2X,E22.16,2X,E22.16,2X,E22.16,2X,I6)") 1,&
       1, 4, NX, NY, i_NoSpecies, &
       domain_size_x, domain_size_y, cell_size, maxstp
  close(25)
  
  stop

2020 FORMAT(E22.16,2X,E22.16,2X,100(E22.16,2X,))
4040 FORMAT(I4,2X,I4,2X,E22.16,2X,E22.16,2X,E22.16,2X,E22.16)
end program fix_grid

subroutine DISCRETIZE(N, CIN, uin, kin, dx, dt, COUT)
  implicit none
  integer, intent(in)           :: N
  double precision, intent(in)  :: CIN(N),uin,kin,dx,dt
  double precision, intent(out) :: COUT(N)
  
  double precision ::  C(N+4), U(N+4), K(N+4)
  
  !! Local Variables
  double precision :: Wind,Diff,DiffusionTerm
  double precision :: AdvectionTerm,AdvectionTermL,AdvectionTermR
  integer          :: shift, i


  shift=2
!  print*,N,uin, kin, dx, dt
  C(1)=CIN(N)
  C(2)=CIN(N-1)
  C(1+shift:N+shift)=CIN(1:N)
  C(N+shift+1)=CIN(1)
  C(N+shift+2)=CIN(2)
  U(1:N+4)=uin
  K(1:N+4)=kin

!  print*,'IN C = ',C
  call ADVDIFF_MIHAI(dt,N,1,dx,U,K,C)
!  print*,'OUT C = ',C

  COUT(1:N)=C(1+shift:N+shift)

end subroutine DISCRETIZE


subroutine ADVDIFF_MIHAI(DT, N, Nspec, DX, U, K,&
     C)
 
  implicit none
  integer, intent(in) :: N,Nspec
  double precision, intent(in)  :: DT, DX, U(1:N+4), K(1:N+4)
  double precision, intent(inout) :: C(1:N+4,Nspec)


  !  Local Variables
  double precision :: C1(1:N+4), DCDX(N)
  integer          :: shift, i

  shift=2
  do i=1,Nspec
     C1=C(:,i)

     call SPACE_ADVDIFF_MIHAI(N, C(:,i), U, K, dx, DCDX)

     C1(1+shift:N+shift)=C(1+shift:N+shift,i)+dt*DCDX(1:N)

     call SPACE_ADVDIFF_MIHAI(N, C1, U, K, dx, DCDX)
     C1(1+shift:N+shift)=C1(1+shift:N+shift)+dt*DCDX(1:N)
     C(1+shift:N+shift,i)=&!C1(1+shift:N+shift)!&
          0.5*(C(1+shift:N+shift,i)+C1(1+shift:N+shift))
     C(1+shift:N+shift,i)=max(C(1+shift:N+shift,i),0.0)
  enddo
  
end subroutine ADVDIFF_MIHAI


subroutine SPACE_ADVDIFF_MIHAI(N, C, U, K, dx, DCDX)
 
  implicit none
  integer, intent(in)           :: N
  double precision, intent(in)  :: C(N+4), U(N+4), K(N+4),dx
  double precision, intent(out) :: DCDX(N)

  
  !! Local Variables
  double precision :: Wind,Diff,DiffusionTerm
  double precision :: AdvectionTerm,AdvectionTermL,AdvectionTermR
  integer          :: shift, i

  shift=2
!!$  if(do_print.eq.1) then
!!$     print*,'dx=',dx
!!$  endif
  do i=1+shift,N+shift
     Wind=(U(i-1)+U(i))/2.0
     if(Wind .ge. 0.0) then
        AdvectionTermL=(1.0/6.0)*((-C(i-2))+5.0*C(i-1)+2.0*C(i))
     else
        AdvectionTermL=(1.0/6.0)*(2.0*C(i-1)+5.0*C(i)+(-C(i+1)))
     endif
     AdvectionTermL=AdvectionTermL*Wind
    
     Wind=(U(i+1)+U(i))/2.0
     if(Wind .ge. 0.0) then
        AdvectionTermR=(1.0/6.0)*((-C(i-1))+5.0*C(i)+2.0*C(i+1))
     else
        AdvectionTermR=(1.0/6.0)*(2.0*C(i)+5.0*C(i+1)+(-C(i+2)))
     endif
     AdvectionTermR=AdvectionTermR*Wind
     
     AdvectionTerm=(AdvectionTermL-AdvectionTermR)/dx
!!$     if(do_print.eq.1) then
!!$        print*,'Advection Term=',AdvectionTerm,AdvectionTermL,AdvectionTermR
!!$     endif
     DiffusionTerm=((K(i-1)+K(i))/2)*(C(i-1)-C(i))-((K(i)+K(i+1))/2)*(C(i)-C(i+1))
     DiffusionTerm=DiffusionTerm/(dx**2)
     
     DCDX(i-shift)=AdvectionTerm+DiffusionTerm
     
  end do
  
end subroutine SPACE_ADVDIFF_MIHAI


