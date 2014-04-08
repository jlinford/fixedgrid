
 module global_vars
   ! ################################################
   real::pi
   parameter(pi=3.141592653589793)
   integer::selectedlb
   integer::do_print
   ! ################################################
   
   ! ################################################
   integer:: i_O3,i_O,i_NO,i_NO2,i_NoSpecies,i_MaxSpecies
   real   :: O3_intial,EmmSource,O_intial,NO_initial,NO2_initial

   integer:: EMMFlag

   ! Species number
   parameter(i_NoSpecies=1)

   !  Species IDs
   parameter(i_O3=1,i_O=2,i_NO=3,i_NO2=4)

   !  Species initial values
   parameter(O3_intial = 8.61E+09)!15)
   parameter(O_intial  = 1.00E-03)
   parameter(NO_intial = 4.47E+08)
   parameter(NO2_intial= 4.47E+09)

   parameter(EmmSource = 4.67E+23)
   ! ################################################
   !  
   real   :: point_source_x,point_source_y
   real   :: max_initial_distance,x_initial_coord,y_initial_coord
   integer:: advtest,MihaiTimestepping
   real   :: u,v,kh

!   parameter(u=5.0)    ! [m/s]
!   parameter(v=0.0)    ! [m/s]
!   parameter(kh=100.0) ! [m/s^2]
   
 end module global_vars
