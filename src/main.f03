program main
    use Convert_trj_PDB !,only : read_pdb, readtrj 
    implicit none
    integer :: icrd
    real(8) :: t1, t2
!    integer :: n_atoms, n_residues, n_chains
!    real(8) :: cellsize(3)
    character(len=100) :: file_ref 
    character(len=100) :: file_trj

    call cpu_time(t1)
    !????????  Caution  ????????  
    !This main and Convert_trj_PDB share the follow variables!!!
    print*, "Input the number of atoms"
    read(*,*) n_atoms
    print*, "Input the number of residues"
    read(*,*) n_residues 
    print*, "Input the number of chains"
    read(*,*) n_chains 
    print*, "Input cell size" 
    read(*,*) cellsize(1:3) 
    !??????????????????????????????????????????????????

    print*, "Input a file name of a reference PDB"
    read(*,"(a120)") file_ref
    print*, "Input a file name of a trajectory"
    read(*,"(a120)") file_trj
 
    call read_pdb(file_ref)
    call readtrj(file_trj)
!    do iatom = 1, n_atoms 
!      write(10,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
!        "ATOM  ", &
!        AtomNum(iatom), &
!        AtomName(iatom),&
!        ResName(iatom),&
!        ChainId(iatom), &
!        ResNum(iatom), &
!        x(iatom), y(iatom), z(iatom),&
!        blnk
!    enddo
!    write(10,"(a6)") "ENDMDL"

    call cpu_time(t2)
    print*, "CPU TIME: ", t2 - t1
end program
