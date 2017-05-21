module restart_to_pdb
  implicit none
contains

  subroutine read_restart_file(unit, file_name, atom_no, cord)
    integer  , intent(in)  :: unit  
    character, intent(in)  :: file_name
    integer  , intent(in)  :: atom_no
    real(4)  , intent(out), allocatable :: cord(:,:)   

    integer :: i

    open(unit, file = file_name, form = "unformatted", status="old")

    read(unit)
    read(unit)
    read(unit)
    read(unit) !(cord(1:3,i),i=1,atom_no)
    read(unit)
     
    !write(unit) (cord(1:3,i),i=1,atom_no)

  end subroutine

end module

program main
  use restart_to_pdb
  implicit none 

  call read_restart_file(11, & 
                         "/Users/siida/tbm2_mountpoint/02_s100b_ctd_work/md0_2_generating_init/n1/md.restart",&
                         atom_no,&
                         cord)


end program
