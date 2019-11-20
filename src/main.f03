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

    write(*, fmt="(a)", advance="no") "*No. of atoms     :"
    read(*,*) Ref%n_atoms

    write(*, fmt="(a)", advance="no") "*Reference PDB    :"
    read(*,"(a120)") file_ref

    write(*, fmt="(a)", advance="no") "*Trajectory file  :"
    read(*,"(a120)") file_trj

    write(*, fmt="(a)", advance="no") "*Output file Name :"
    read(*,*) fnameOut

    call read_pdb(file_ref, Ref)

    call analyze_trj(file_trj, fnameOut, Ref, Trj)

    call cpu_time(t2)
    print*, "CPU TIME (sec): ", t2 - t1

end program
