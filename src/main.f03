program main
    use variables
    use Convert_trj_PDB !,only : read_pdb, readtrj 
    implicit none
    real(8) :: t1, t2
    type(var_PDB) :: Ref
    type(var_snapshot) :: Trj
    character(len=100) :: file_ref
    character(len=100) :: file_trj
    character(len=100) :: fnameOut
    logical :: IsOutput

    call cpu_time(t1)

    print*, "*Output?"
    read(*,*) IsOutput

    print*, "*No. of atoms    :"
    read(*,*) Ref%n_atoms

    print*, "*No. of residues :"
    read(*,*) Ref%n_residues 

    print*, "*No. of chains   :"
    read(*,*) Ref%n_chains 

!add 2017-2-22
!    print*, "*The ranges of atom no.: "
!    do i = 1, nchains
!      read(*,*) Ref%atom_ranges(,)
!    enddo

    print*, "*Cell size       :" 
    read(*,*) Trj%cellsize(1:3) 

    print*, "*Reference PDB   :"
    read(*,"(a120)") file_ref

    print*, "*Trajectory file :"
    read(*,"(a120)") file_trj

    if (IsOutput) then 
        print*,"*Trajectory output:", IsOutput
        read(*,*) fnameOut
        print*,"*Output file Name :", fnameOut
    endif

    call read_pdb(file_ref, Ref)

    call analyze_trj(file_trj, fnameOut, Ref, Trj, IsOutput)

    call cpu_time(t2)
    print*, "CPU TIME (sec): ", t2 - t1

end program
