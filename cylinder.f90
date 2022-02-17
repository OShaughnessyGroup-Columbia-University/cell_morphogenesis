!**********************************************************************
! Program: Axisymmetric tubular shape. Finite volumes
! Serial code



!************************************************************************
!   MODULE CONSTANTS
!************************************************************************
module constants

double precision, parameter  :: pi = 3.1415926535897d0 ! Number pi

double precision, parameter  ::	R0 = 4.d0  ! Initial cylinder radius
double precision, parameter  ::	L0 = 6.d0  ! Initial arc-length domain

double precision, parameter  ::	cstar = 3.d0    ! Density at which stress saturates
double precision, parameter  ::	c0 = 1.d0    ! Preferred concentration of the stress regulator
double precision, parameter  ::	beta = 500.d0  ! Contractility
double precision, parameter  ::	s_tens = 1.d0  ! Constant surface tension
double precision, parameter  ::	diff = 1.d-1  ! Translational diffusion

double precision, parameter  :: bend_mod = 0.d0 ! Bending modulus\
double precision, parameter  :: drag = 10.d0 ! Drag coefficient
double precision, parameter  :: eta = 20.d0 ! 2-D viscosity
double precision, parameter  :: kturn = 2.d-1 ! Inverse of turnover time

integer, parameter           :: Nu = 140  ! Number of equally spaced collocation points
double precision, parameter  :: du = dble(L0/Nu) ! Distance between neighboring collocation points

end module 


!************************************************************************
!   MAIN PROGRAM
!************************************************************************

program tubular
 use constants
 implicit none

 double precision                            :: dt,dd,pressure,exponent
 integer                                     :: cycles,i,j
 double precision,dimension(:),allocatable   :: u,c,r,alpha,h,css,ctt,gamma,f,du_f,du_h,du_css,duu_css
 double precision,dimension(:),allocatable   :: a1,a2,a3,a4,a5,b1,b2,b3,RHS,vs,vn,ss,z
 double precision,dimension(:),allocatable   :: r_p,h_p,css_p,ctt_p,gamma_p,c_p,du_vs,du_vn,duu_vn,du_c,duu_c
 double precision,dimension(:,:),allocatable :: matrix
 integer,dimension(:),allocatable            :: indice
 integer*8                                   :: itime,icycles,ifile
 character*32 filename1

 allocate(u(Nu),c(Nu),r(Nu),alpha(Nu),h(Nu),css(Nu),ctt(Nu),gamma(Nu),f(Nu),du_f(Nu),du_h(Nu),du_css(Nu),duu_css(Nu))
 allocate(a1(Nu),a2(Nu),a3(Nu),a4(Nu),a5(Nu),b1(Nu),b2(Nu),b3(Nu),RHS(2*Nu+1),vs(Nu),vn(Nu),ss(Nu),z(Nu))
 allocate(r_p(Nu),h_p(Nu),css_p(Nu),ctt_p(Nu),gamma_p(Nu),c_p(Nu),du_vs(Nu),du_vn(Nu),duu_vn(Nu),du_c(Nu),duu_c(Nu))
 allocate(matrix(2*Nu+1,2*Nu+1),indice(2*Nu+1))

! Time step
dt=1.d-6
cycles=2

! Coordinates of collocation points: located in the center of the finite volumes
do i=1,Nu
 u(i)=(i-5.d-1)*du
enddo

exponent=4.d0

! Initial values: Concentration field c(u,t=0), radial coordinate r(u,t=0), tangent angle alpha(u,t=0)
! parameter transformation h(u,t=0), curvature css(u,t=0), curvature ctt(u,t=0), christoffel symbol gamma(u,t=0)

do i=1,Nu
!  c(i)=c0
 c(i)=c0+1.d-3*sin(2.d0*pi*u(i)/L0)
! c(i)=c0+1.d-3*exp(-(u(i)-L0/2.d0)**2.d0/4.d-2)
 r(i)=R0
 alpha(i)=pi/2.d0
 h(i)=1.d0
 css(i)=0.d0
 ctt(i)=1.d0/r(i)
 gamma(i)=0.d0
enddo

! Linear function for the active stress
do i=1,Nu
 f(i)=erf(c(i)/cstar)
! c(i)**2.d0/(c(i)**2.d0+c0**2.d0)
enddo
! Derivatives of f(u,t),h(u,t),css(u,t)
call first_deriv(f,du_f)
call first_deriv(h,du_h)
call first_deriv(css,du_css)
call second_deriv(css,duu_css)

! Parameters used to construct the hydrodynamic matrix
do i=1,Nu
 a1(i)=2.d0*eta/(h(i)**2.d0)
 a2(i)=(2.d0*eta/h(i))*(gamma(i)-(du_h(i)/(h(i)**2.d0)))
 a3(i)=-(2.d0*eta*gamma(i)**2.d0)-drag
 a4(i)=2.d0*eta*(du_css(i)/h(i)+gamma(i)*css(i)-ctt(i)*gamma(i))
 a5(i)=2.d0*eta*css(i)/h(i)
 b1(i)=2.d0*eta*css(i)/h(i)
 b2(i)=2.d0*eta*ctt(i)*gamma(i)
 b3(i)=2.d0*eta*(css(i)**2.d0+ctt(i)**2.d0) 
enddo

do i=1,2*Nu+1
 RHS(i)=0
 do j=1,2*Nu+1
  matrix(i,j)=0
 enddo
enddo

! Coefficients of the hydrodynamic matrix. Linear problem solving tangential and normal force balance equations
do i=1,Nu
  if(i.eq.1) then 
   matrix(i,Nu)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
   matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
   matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
   matrix(i,Nu+Nu)=-(a5(i)/(2.d0*du))
   matrix(i,Nu+i)=a4(i)
   matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
   RHS(i)=-(beta/h(i))*du_f(i)   
   matrix(Nu+i,Nu)=-b1(i)/(2.d0*du)   
   matrix(Nu+i,i)=b2(i)   
   matrix(Nu+i,i+1)=b1(i)/(2.d0*du)   
   matrix(Nu+i,Nu+i)=b3(i)   
   RHS(Nu+i)=-(s_tens+beta*f(i))*(css(i)+ctt(i))+bend_mod*(ctt(i)+css(i))*(ctt(i)-css(i))**2.d0+&
   2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-(du_css(i)*du_h(i)/(h(i)**3.d0))+&
   (2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
  elseif(i.eq.Nu) then
   matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
   matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
   matrix(i,1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
   matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
   matrix(i,Nu+i)=a4(i)
   matrix(i,Nu+1)=(a5(i)/(2.d0*du))
   RHS(i)=-(beta/h(i))*du_f(i) 
   matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)   
   matrix(Nu+i,i)=b2(i)
   matrix(Nu+i,1)=b1(i)/(2.d0*du)   
   matrix(Nu+i,Nu+i)=b3(i)
   RHS(Nu+i)=-(s_tens+beta*f(i))*(css(i)+ctt(i))+bend_mod*(ctt(i)+css(i))*(ctt(i)-css(i))**2.d0+&
   2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-(du_css(i)*du_h(i)/(h(i)**3.d0))+&
   (2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
  else     
   matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
   matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
   matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
   matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
   matrix(i,Nu+i)=a4(i)
   matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
   RHS(i)=-(beta/h(i))*du_f(i)
   matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)  
   matrix(Nu+i,i)=b2(i)
   matrix(Nu+i,i+1)=b1(i)/(2.d0*du)  
   matrix(Nu+i,Nu+i)=b3(i)
   RHS(Nu+i)=-(s_tens+beta*f(i))*(css(i)+ctt(i))+bend_mod*(ctt(i)+css(i))*(ctt(i)-css(i))**2.d0+&
   2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-(du_css(i)*du_h(i)/(h(i)**3.d0))+&
   (2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
  endif
enddo
RHS(2*Nu+1)=0.d0
do i=1,Nu
 matrix(Nu+i,2*Nu+1)=-1.d0   
 matrix(2*Nu+1,Nu+i)=r(i)*h(i)*du    
enddo

! Solution of the system of equations. Matrix is destroyed. Solution is stored in RHS vector
 call ludcmp(matrix,2*Nu+1,2*Nu+1,indice,dd)
 call lubksb(matrix,2*Nu+1,2*Nu+1,indice,RHS)

! We store the values of vs, vn and p.
 do i=1,Nu
  vs(i)=RHS(i)
  vn(i)=RHS(Nu+i)
 enddo
 pressure=RHS(2*Nu+1)

! Obtain the arc-length 
ss(1)=h(1)*du/2.d0
do i=2,Nu
 ss(i)=ss(i-1)+(h(i-1)+h(i))*du/2.d0
enddo
  
! Obtain z coordinate
z(1)=r(1)*ctt(1)*h(1)*du/2.d0
do i=2,Nu
 z(i)=z(i-1)+(r(i-1)*ctt(i-1)*h(i-1)+r(i)*ctt(i)*h(i))*du/2.d0
enddo  

! Time marching
do icycles=1,cycles
 do itime=1,200000000

  ! Some derivatives
  call first_deriv(vs,du_vs)
  call first_deriv(vn,du_vn)
  call second_deriv(vn,duu_vn)
  call first_deriv(h,du_h)
  call first_deriv(c,du_c)
  call first_deriv(css,du_css)
  call second_deriv(c,duu_c)

  ! Obtain values for the next time step
  do i=1,Nu
   r_p(i)=r(i)+dt*vn(i)*r(i)*ctt(i)
   h_p(i)=h(i)+dt*h(i)*css(i)*vn(i)
   css_p(i)=css(i)+dt*(-(css(i)**2.d0)*vn(i)-(duu_vn(i)/(h(i)**2.d0))+(du_h(i)*du_vn(i)/h(i)**3.d0))
   ctt_p(i)=ctt(i)+dt*(-(ctt(i)**2.d0)*vn(i)-gamma(i)*du_vn(i)/h(i))
   gamma_p(i)=gamma(i)+dt*(ctt(i)*du_vn(i)/h(i)-ctt(i)*gamma(i)*vn(i))
   c_p(i)=c(i)+dt*(-(c(i)*du_vs(i)/h(i))-(vs(i)*du_c(i)/h(i))-gamma(i)*c(i)*vs(i)-c(i)*vn(i)*(css(i)+ctt(i))+&
   diff*((duu_c(i)/h(i)**2.d0)-du_h(i)*du_c(i)/h(i)**3.d0+gamma(i)*du_c(i)/h(i))-(c(i)-c0)*kturn)
  enddo

  ! Update values
  do i=1,Nu
   r(i)=r_p(i)
   h(i)=h_p(i)
   css(i)=css_p(i)
   ctt(i)=ctt_p(i)
   gamma(i)=gamma_p(i) 
   c(i)=c_p(i)
   f(i)=erf(c(i)/cstar)
! c(i)**2.d0/(c(i)**2.d0+c0**2.d0)
  enddo

  ! Derivatives of f(u,t),h(u,t),css(u,t)
  call first_deriv(f,du_f)
  call first_deriv(h,du_h)
  call first_deriv(css,du_css)
  call second_deriv(css,duu_css)
  
  ! Parameters used to construct the hydrodynamic matrix
  do i=1,Nu
   a1(i)=2.d0*eta/(h(i)**2.d0)
   a2(i)=(2.d0*eta/h(i))*(gamma(i)-(du_h(i)/(h(i)**2.d0)))
   a3(i)=-(2.d0*eta*gamma(i)**2.d0)-drag
   a4(i)=2.d0*eta*(du_css(i)/h(i)+gamma(i)*css(i)-ctt(i)*gamma(i))
   a5(i)=2.d0*eta*css(i)/h(i)
   b1(i)=2.d0*eta*css(i)/h(i)
   b2(i)=2.d0*eta*ctt(i)*gamma(i)
   b3(i)=2.d0*eta*(css(i)**2.d0+ctt(i)**2.d0) 
  enddo

  do i=1,2*Nu+1
   RHS(i)=0
   do j=1,2*Nu+1
    matrix(i,j)=0
   enddo
  enddo

  ! Coefficients of the hydrodynamic matrix. Linear problem solving tangential and normal force balance equations
  do i=1,Nu
   if(i.eq.1) then 
    matrix(i,Nu)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
    matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
    matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
    matrix(i,Nu+Nu)=-(a5(i)/(2.d0*du))
    matrix(i,Nu+i)=a4(i)
    matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
    RHS(i)=-(beta/h(i))*du_f(i)   
    matrix(Nu+i,Nu)=-b1(i)/(2.d0*du)   
    matrix(Nu+i,i)=b2(i)   
    matrix(Nu+i,i+1)=b1(i)/(2.d0*du)   
    matrix(Nu+i,Nu+i)=b3(i)   
    RHS(Nu+i)=-(s_tens+beta*f(i))*(css(i)+ctt(i))+bend_mod*(ctt(i)+css(i))*(ctt(i)-css(i))**2.d0+&
    2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-(du_css(i)*du_h(i)/(h(i)**3.d0))+&
    (2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
   elseif(i.eq.Nu) then
    matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
    matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
    matrix(i,1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
    matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
    matrix(i,Nu+i)=a4(i)
    matrix(i,Nu+1)=(a5(i)/(2.d0*du))
    RHS(i)=-(beta/h(i))*du_f(i) 
    matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)   
    matrix(Nu+i,i)=b2(i)
    matrix(Nu+i,1)=b1(i)/(2.d0*du)   
    matrix(Nu+i,Nu+i)=b3(i)
    RHS(Nu+i)=-(s_tens+beta*f(i))*(css(i)+ctt(i))+bend_mod*(ctt(i)+css(i))*(ctt(i)-css(i))**2.d0+&
    2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-(du_css(i)*du_h(i)/(h(i)**3.d0))+&
    (2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
   else     
    matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
    matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
    matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
    matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
    matrix(i,Nu+i)=a4(i)
    matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
    RHS(i)=-(beta/h(i))*du_f(i)
    matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)  
    matrix(Nu+i,i)=b2(i)
    matrix(Nu+i,i+1)=b1(i)/(2.d0*du)  
    matrix(Nu+i,Nu+i)=b3(i)
    RHS(Nu+i)=-(s_tens+beta*f(i))*(css(i)+ctt(i))+bend_mod*(ctt(i)+css(i))*(ctt(i)-css(i))**2.d0+&
    2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-(du_css(i)*du_h(i)/(h(i)**3.d0))+&
    (2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
   endif
  enddo
  RHS(2*Nu+1)=0.d0
  do i=1,Nu
   matrix(Nu+i,2*Nu+1)=-1.d0   
   matrix(2*Nu+1,Nu+i)=r(i)*h(i)*du    
  enddo

  ! Solution of the system of equations. Matrix is destroyed. Solution is stored in RHS vector
  call ludcmp(matrix,2*Nu+1,2*Nu+1,indice,dd)
  call lubksb(matrix,2*Nu+1,2*Nu+1,indice,RHS)

  ! We store the values of vs, vn and p.
  do i=1,Nu
   vs(i)=RHS(i)
   vn(i)=RHS(Nu+i)
  enddo
  pressure=RHS(2*Nu+1)

  ! Obtain the arc-length 
  ss(1)=h(1)*du/2.d0
  do i=2,Nu
   ss(i)=ss(i-1)+(h(i-1)+h(i))*du/2.d0
  enddo
  
  ! Obtain z coordinate
  z(1)=r(1)*ctt(1)*h(1)*du/2.d0
  do i=2,Nu
   z(i)=z(i-1)+(r(i-1)*ctt(i-1)*h(i-1)+r(i)*ctt(i)*h(i))*du/2.d0
  enddo  

  ! Save variables from time to time
  if (10000.0*(itime/10000).eq.1.0*itime) then
   ifile = itime/10000
   write(filename1,'(A4,I6,A3,I2,A4)') 'data',ifile,'cyc',icycles,'.txt'
   open(unit=5,file=filename1,status='unknown')
   write(5,*) 'Variables = "vs","vn","h","r","css","ctt","gamma","ss","c","z" '
   do i=1,Nu
    write(5,*) vs(i),vn(i),h(i),r(i),css(i),ctt(i),gamma(i),ss(i),c(i),z(i)
   enddo
   close(5)
  endif
    
 enddo
enddo

end

!-------------------------------------------------------------------------------------------------------------------------------------------------
! Subroutines that solve the hydrodynamic linear system of equations. 
! This routine is used in combination with lubksb to solve linear equations or invert a matrix. We solve the linear set of equations A · x = b
! call ludcmp(a,n,np,indx,d)
! call lubksb(a,n,np,indx,b)
! The answer x will be returned in b. Your original matrix A will be destroyed.
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine ludcmp(a,n,np,indx,d)

 integer n,np,indx(n),NMAX
 double precision d,a(np,np),TINY
 PARAMETER (NMAX=40000,TINY=1.0e-20) ! Largest expected n, and a small number.
 integer i,imax,j,k
 double precision aamax,dum,sum,vv(NMAX)

 d=1.d0

 do i=1,n
  aamax=0.d0
  do j=1,n
   if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
  enddo
  if (aamax.eq.0.)  then
!   write (*,*) 333	
  endif		
  vv(i)=1.d0/aamax
 enddo
 do j=1,n
  do i=1,j-1
   sum=a(i,j)
   do k=1,i-1
    sum=sum-a(i,k)*a(k,j)
   enddo
   a(i,j)=sum
  enddo
  aamax=0.d0
  do i=j,n
   sum=a(i,j)
   do k=1,j-1
    sum=sum-a(i,k)*a(k,j)
   enddo
   a(i,j)=sum
   dum=vv(i)*abs(sum)
   if (dum.ge.aamax) then
    imax=i
    aamax=dum
   endif
  enddo
  if (j.ne.imax)then
   do k=1,n
    dum=a(imax,k)
    a(imax,k)=a(j,k)
    a(j,k)=dum
   enddo
   d=-d
   vv(imax)=vv(j)
  endif
  indx(j)=imax
  if(a(j,j).eq.0.) a(j,j)=TINY
  if(j.ne.n)then
   dum=1./a(j,j)
   do i=j+1,n
    a(i,j)=a(i,j)*dum
   enddo
  endif
 enddo

 return

end subroutine ludcmp

subroutine lubksb(a,n,np,indx,b)

 integer n,np,indx(n)
 double precision a(np,np),b(n)
 integer i,ii,j,ll
 double precision sum

 ii=0

 do i=1,n
  ll=indx(i)
  sum=b(ll)
  b(ll)=b(i)
  if (ii.ne.0)then
   do j=ii,i-1
    sum=sum-a(i,j)*b(j)
   enddo 
  else if (sum.ne.0.) then
   ii=i
  endif
  b(i)=sum
 enddo
 do i=n,1,-1
  sum=b(i)
  do j=i+1,n
   sum=sum-a(i,j)*b(j)
  enddo
  b(i)=sum/a(i,i)
 enddo

 return

end subroutine lubksb


!-------------------------------------------------------------------------------------------------------------------------------------------------
! It calculates the first derivative of vect and store the results in dx_vect. 2nd order accurate. Periodic BCs. 
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine first_deriv(vect,dx_vect)
 use constants
 implicit none

 integer i
 double precision vect(Nu), dx_vect(Nu) 

 do i=1,Nu
  if(i.eq.1) then
   dx_vect(i)=(vect(i+1)-vect(Nu))/(2.d0*du)
  elseif(i.eq.Nu) then
   dx_vect(i)=(vect(1)-vect(i-1))/(2.d0*du)
  else
   dx_vect(i)=(vect(i+1)-vect(i-1))/(2.d0*du)
  endif
 enddo

end subroutine first_deriv


!-------------------------------------------------------------------------------------------------------------------------------------------------
! It calculates the second derivative of vect and store the results in dxx_vect. 2nd order accurate. Periodic BCs. 
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine second_deriv(vect,dxx_vect)
 use constants
 implicit none

 integer i
 double precision vect(Nu), dxx_vect(Nu) 

 do i=1,Nu
  if(i.eq.1) then
   dxx_vect(i)=(vect(i+1)-2.d0*vect(i)+vect(Nu))/(du**2.d0);
  elseif(i.eq.Nu) then
   dxx_vect(i)=(vect(1)-2.d0*vect(i)+vect(i-1))/(du**2.d0);
  else
   dxx_vect(i)=(vect(i+1)-2.d0*vect(i)+vect(i-1))/(du**2.d0);
  endif
 enddo

end subroutine second_deriv








