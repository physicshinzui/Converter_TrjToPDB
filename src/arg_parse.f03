module arg_parse 
!program main
implicit none

  character(len=*), parameter :: version = '1.0'
  character(len=32) :: arg
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  logical :: do_time = .false.
  logical :: IsRef   = .false.
  logical :: IsTrj   = .false.
  character(len=100) :: file_ref, file_trj
  integer :: i, length
contains

  subroutine argument_parse()
      do i = 1, command_argument_count()
         call get_command_argument(i, arg, length=length)

         select case (arg)
         case ('-v', '--version')
            print '(2a)', 'version ', version
            stop
         case ('-h', '--help')
            call print_help()
            stop
         case ('-r', '--reference')
            IsRef = .true.
         case ('-t', '--trajectory')
            IsTrj = .true.
         case default
            if (IsRef) then 
              file_ref=arg 
              print*, "Reference : ", file_ref
              IsRef = .false.
            elseif (IsTrj) then
              file_trj=arg
              print*, "Trajectroy: ", file_trj
              IsTrj = .false.
            endif
            !print '(a,a,/)', 'Unrecognized command-line option: ', arg
            !call print_help()
            !stop
         end select
      end do
      print*, file_ref
  end subroutine

  subroutine print_help()
    print '(a)', 'usage: cmdline [OPTIONS]'
    print '(a)', ''
    print '(a)', 'Without further options, cmdline prints the date and exits.'
    print '(a)', ''
    print '(a)', 'options:'
    print '(a)', ''
    print '(a)', '  -v, --version     print version information and exit'
    print '(a)', '  -h, --help        print usage information and exit'
    print '(a)', '  -t, --time        print time'
  end subroutine print_help

end module
