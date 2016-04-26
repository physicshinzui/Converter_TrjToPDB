program main
    use variables
    !use arg_parse
    use Convert_trj_PDB !,only : read_pdb, readtrj 
    implicit none
    real(8) :: t1, t2
    type(var_PDB) :: Ref
    type(var_Trj) :: Trj
    character(len=100) :: file_ref 
    character(len=100) :: file_trj

    call cpu_time(t1)

    !call argument_parse()
    print*, "*No. of atoms    :"
    read(*,*) Ref%n_atoms
    print*, "*No. of residues :"
    read(*,*) Ref%n_residues 
    print*, "*No. of chains   :"
    read(*,*) Ref%n_chains 
    print*, "*Cell size       :" 
    read(*,*) Trj%cellsize(1:3) 
    print*, "*Reference PDB   :"
    read(*,"(a120)") file_ref
    print*, "*Trajectory file :"
    read(*,"(a120)") file_trj

    call read_pdb(file_ref, Ref)
    call outputPDB(12,Ref)
    call analyze_trj(file_trj, Ref, Trj, .false.)

    call cpu_time(t2)
    print*, "CPU TIME: ", t2 - t1
end program
