! KDM ken.mankoff@giss.nasa.gov
! Helper routines

      module lime
      implicit none      
  
      contains

      subroutine dbgcsv(msg,val,i,j,root)
! Print "msg,val" statements with ISO date string

! Example usage:
! In rundeck:
! use LIME
! call dbgcsv(msg="float value", val=42.0)
! call dbgcsv(val=42 * 1.0, msg="convert int to float")
! call dbgcsv(msg="IJ supported", val=v i=II, j=JJ)
! call dbgcsv(msg="limited MPI support", val=v, root=.true.)
!
! To turn off, either 
! #define dbgcsv_silent ! rundeck
! EXTRA_FFLAGS="-Ddbgcsv_silent" ! make setup

      use model_com, only: modelEclock
      use BaseTime_mod, only: BaseTime
      use Rational_mod, only : Real
      Use Domain_decomp_atm, only: AM_I_ROOT
    
      implicit none
    
      character(len=*), intent(in) :: msg
      real*8, optional, intent(in) :: val
      integer, optional, intent(in) :: i,j
      logical, optional, intent(in) :: root

      character(len=127) :: tmp_str
      character(len=127) :: out_str

      integer :: y, m, d, h, mm, s
      real*8 rs                 ! seconds
      type (BaseTime) :: t      ! time

#ifdef dbgscv_silent
      return
#endif
      
      y = modelEclock%getYear()
      m = modelEclock%getMonth()
      d = modelEclock%getDate()
      t = modelEclock%getTimeInSecondsFromDate(y,m,d,0)
      rs = Real(t)

      h = int(int(rs)/3600)
      mm = int((int(rs)-h*3600)/60)
      s = int(rs)-(h*3600)-(mm*60)

! iso8601
      write(tmp_str,
     *     '(I0.4,"-",I0.2,"-",I0.2,"T",I0.2,":",I0.2,":",I0.2)')
     *     y,m,d,h,mm,s
      out_str = tmp_str
      
      if (present(i).and.present(j)) then
         write(tmp_str, '("i,", I0.3, ",j,", I0.3)') i,j
         out_str = trim(out_str)//","//tmp_str
      end if

      out_str = trim(out_str)//","//msg

      write(tmp_str, '(F0.5)') val
      out_str = trim(out_str)//","//tmp_str

      out_str = "DBGCSV,"//trim(out_str)
      out_str = trim(out_str)

      if ((present(root) .and. root) .or. am_i_root()) then
         print '(a)', trim(out_str)
      else
         out_str = trim(out_str)//",not-root"
         print '(a)', trim(out_str)
      endif
      
      end subroutine dbgcsv
      
      end module lime
