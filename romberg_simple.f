      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
c
c     Robert Thacker
c     Romberg Integration
c     integrate a function - in this case natural log of 10.
c
      dimension T(0:50,51)
c
      open(unit=10,file='natural_log.out')
      open(unit=11,file='T_1.tex')
      open(unit=13,file='T_12.tex')
      open(unit=14,file='T_8901234.tex')
      open(unit=15,file='T_15_to_21.tex')
c
c
c==========================================================
c     Preamble    equation 4
c     set limits of integration and convergence value
c
      a         =   1.q0
      b         =  10.q0
      exact     = log(b)
c
      h         = b - a
c
      tol       = 1.d-15
      niter_max = 26  
c
      write(10,*) 'a   = ',a
      write(10,*) 'b   = ',b
      write(10,*) 'h   = ',h
      write(10,*) 'tol = ',tol
c
      write(10,800)
 800  format(/,
     &'   Best solution values for each n',/
     &'  (i - 1)/i are the rows and j is the oolumn of T',/
     &'  n','  i','  j','   T(i,j)/T(i,j) ','      difference  ' )
c==========================================================
c
c
c 
      do 1 i = 0,50
      do 2 j = 1,51
      T(i,j) = 0.d0
 2    continue
 1    continue
c
c=====================================================
c
c     equation 10 
      n = 0

      begin  = y_nat_log(a)

      end    = y_nat_log(b)

      T(n,1) = h*(( begin + end  )/2.q0)
      write(6,*) 'n = ',n    ! write to screen

c-------------------------------------------------
c     begin romberg iter n = 1,niter_max first column of T
c
c     Trapezoidal Rule with interval halving

      do 10 n =1,niter_max
      write(6,*) 'n = ',n    ! write to screen
   
c

      first   = h/( 2.d0**( n - 1 ) )  ! see equation 11
      iend    = 2**( n ) - 1
      denom   = 2.d0**n
c
      value   = 0.d0
c
      do 20 i = 1,iend,2

      x = (a + ( h*dfloat(i) )/denom )

      value = value + y_nat_log(x)
 20   continue
c
      T( n,1 ) = .5d0*( T(n - 1,1 ) + first*value )  ! equation 11
c
c-----------------------------------------------------
c-    perform richardson extrapolation
c     equation 12
c
      nrow     = n           ! start with this row of T
      jcol_end = nrow + 1    ! this is the last column of T
c
      do 30 j  = 2,jcol_end  ! j is which column - 
      nrow = nrow - 1
      T(nrow,j) = ( 4.d0**(j-1)*T(nrow+1,j-1) - T(nrow,j-1) ) /  
     &             (4.d0**(j-1) -1.d0)
 30   continue

c------------------------------------------------------
c
c     check for convergence    
c
      check_converge = 1.0d10
      irow     = n


      if( n .le. 3 )  go to 10
c
      do 35 j = 1, n - 1    ! go across columns 1 through n - 1
c
c 
      if( dabs( T(irow,j) - T(irow-1,j) ) .le. tol )    then   ! we have converged
c
c
      write(6,*) 
     &'solution has converged: see natural_log.out for details'
      write(10,810) n,
     &              irow - 1,j,irow - 1,j,T(irow-1,j),
     &              irow    ,j,irow,    j,T(irow  ,j),
     &              irow,j,irow-1,j,
     &             ( T(irow,j)- T(irow-1,j) ),( T(irow,j)- T(irow-1,j) )
 810  format(/,
     &'the run has converged:',/
     &       'n        = ',i5,/
     &       'irow-1,j = ',2i5,' T(',i2,',',i2,') = ',f20.15,/
     &       'irow  ,j = ',2i5,' T(',i2,',',i2,') = ',f20.15,/
     &       '           T(',i2,',',i2,') - T(',i2,',',i2,') = ',f20.15,
     &              ' = ',1pd22.15)
c 
      write(10,820)   T(irow,j),irow,j

 820  format(36x,'1 234567890123456789012345678901234567890',/
     &       17x,'solution      = ',f20.15,' =   T(',i2,',',i2,')' )
c
c
      go to 799

                                                        endif !if( dabs( T(irow,j) - T(irow-1,j) ) 
c==============================================================
c
c     look for best solution in T fr this value of n
c
      check = T(irow,j) - T(irow-1,j) 

      if( dabs(check) .lt. dabs(check_converge)   )     then
      check_converge = check
      keep_irow      = irow
      keep_j         =    j
                                                        endif

c==============================================================
c
      irow = irow - 1
 35   continue
c
c===============================================================
c
c     print out best solution for this n
c
      write(10,805) n,keep_irow - 1,keep_j,T(keep_irow - 1,keep_j),
     &                check,
     &                keep_irow,keep_j,T(keep_irow,keep_j)
 805  format(3i3,f20.15,f20.15,/
     &    3x,2i3,f20.15)
c

c      write(10,910) n,first,iend,denom,keep_irow,keep_j,
c     &check_converge
 910  format('n cpu_time first iend denom i j = ',i5,f20.10,i15,
     &     f10.0,i3,i3,1pd20.12)
      nn = n
 10   continue
c
c-----------------------------------------------------------------------
c
 799  continue

      end
c================================
c
c     natural log
c
      real*8 function y_nat_log(x)
      implicit real*8(a-h,o-z)
c
      y_nat_log = 1.d0/x
      return
      end
c================================
      real*8 function xi_s(s_val,u)
      implicit real*8(a-h,o-z)
c     page 59   equation 32
c
      xi_s = u**s_val/(dsinh(2.*u) + 2.*u)
      return
      end
      real*8 function xj_s(s_val,u)
      implicit real*8(a-h,o-z)
c     page 59   equation 32
c
      xj_s = ( (u**s_val)*exp(-2.*u) )/(dsinh(2.*u) + 2.*u)
      return
      end
