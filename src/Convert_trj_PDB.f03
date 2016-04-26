module Convert_trj_PDB 
  use variables 
  implicit none

  public  read_pdb
  public  analyze_trj 
!  private outputPDB

contains

  subroutine read_pdb(FilName, PDB)
    type(var_PDB), intent(inout) :: PDB 
    character(len=*), intent(in) :: FilName
    character(len=6)   :: atom
    integer, parameter :: FilUnit = 11
    integer, parameter :: OutPDBUnit  = 12
    integer :: i, ios, ter 
    integer :: NotAtomLine,NumAtom 
  
    open(FilUnit, file=FilName, status="old")
    print*,"  File Name: ",FilName
    allocate(PDB%AtomNum(PDB%n_atoms),PDB%AtomName(PDB%n_atoms),PDB%ResName(PDB%n_atoms), &
             PDB%ChainId(PDB%n_atoms),PDB%ResNum(PDB%n_atoms),PDB%x(PDB%n_atoms),PDB%y(PDB%n_atoms),PDB%z(PDB%n_atoms))
    PDB%AtomNum(:) = 0
    PDB%ResNum(:)  = 0
    PDB%x(:) = 0; PDB%y(:) = 0; PDB%z(:) = 0
  
  !***Detect first atom line
    do
      read(FilUnit,*) atom
      if(atom /= "ATOM" .and. atom /= "HETATM") then
        NotAtomLine = NotAtomLine + 1  
        cycle
      elseif(atom == "ATOM" .or. atom == "HETATM") then
        exit
      endif
    enddo
    backspace(FilUnit)
  !***start reading
    do i = 1, PDB%n_atoms 
      read(FilUnit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") &
           atom       , &
           PDB%AtomNum(i) , &
           PDB%AtomName(i), &
           PDB%ResName(i) , &
           PDB%ChainId(i) , &
           PDB%ResNum(i)  , &
           PDB%x(i),        &
           PDB%y(i),        &
           PDB%z(i)  
           if (atom == "TER") stop "TER lines are included, remove these."
    enddo
    close(FilUnit)
  end subroutine

!------------------------------------------------------
  subroutine analyze_trj(filename, PDB, A_snapshot, isOutputTrj)
      use apply_PBC 
      implicit none
      type(var_PDB), intent(inout)   :: PDB 
      type(var_Trj), intent(inout)   :: A_snapshot
      logical, optional :: isOutputTrj
      character(len=*) :: filename
      integer    :: unit = 11, unit_outPDB = 10
      integer    :: ios, iatom, iconf
      integer(4) :: istp,iyn15v,iyn15h
      real(4)    :: sitime, sec, et, kinetic, temperature, rmsf, rmsd 
      if (.not. present(isOutputTrj)) isOutputTrj = .false.

  
      open(unit, file = filename, form="unformatted", status="old")
      allocate(A_snapshot%x(PDB%n_atoms),A_snapshot%y(PDB%n_atoms),A_snapshot%z(PDB%n_atoms))
      
      !***Initialization
      iconf = 0; A_snapshot%x(:) = 0; A_snapshot%y(:) = 0; A_snapshot%z(:) = 0
      
      call prepare_apply_PBC(PDB%n_atoms,PDB%ResNum,PDB%n_residues)
      !***Reading trajectory
      do 
          read(unit, iostat=ios) istp,sitime,sec,et,kinetic,temperature,&
                              A_snapshot%potential,rmsf,iyn15v,iyn15h,rmsd
          if (ios /= 0) exit
          read(unit) (A_snapshot%x(iatom), A_snapshot%y(iatom), A_snapshot%z(iatom), iatom = 1, PDB%n_atoms)
          iconf = iconf + 1
          call ReturnAtom(PDB%n_atoms,PDB%n_chains,PDB%n_residues,A_snapshot%x, A_snapshot%y, A_snapshot%z, A_snapshot%cellsize)
          if (isOutputTrj) call outputPDB(unit_outPDB, PDB, A_snapshot, iconf)
          print*,"Conf NO:", iconf 
      enddo
      111 close(unit_outPDB)
      A_snapshot%n_confs = iconf
      write(*,'("#Number of conformation ",i8)') A_snapshot%n_confs 
  end subroutine
!
  subroutine outputPDB(unit, PDB, a_snapshot, iconf)
      integer      , intent(in) :: unit 
      integer :: i
      type(var_PDB), intent(in) :: PDB 
      type(var_Trj), optional, intent(in) :: a_snapshot 
      integer, optional, intent(in) :: iconf

      if (present(a_snapshot)) then 
          open(unit, file ="frames.pdb")!, status="replace")
          write(unit,"('MODEL', i7)") iconf
          write(unit,"('#Potential, kcal/mol ', f10.3)") a_snapshot%potential 
          do i = 1, PDB%n_atoms 
            write(unit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
              "ATOM  ", &
              PDB%AtomNum(i), &
              PDB%AtomName(i),&
              PDB%ResName(i),&
              PDB%ChainId(i), &
              PDB%ResNum(i), &
              a_snapshot%x(i), a_snapshot%y(i), a_snapshot%z(i),&
              PDB%blnk
          enddo
          write(unit,"(a6)") "ENDMDL"
      else
          open(unit, file ="Reference.pdb")!, status="replace")
          do i = 1, PDB%n_atoms 
            write(unit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
              "ATOM  ", &
              PDB%AtomNum(i), &
              PDB%AtomName(i),&
              PDB%ResName(i),&
              PDB%ChainId(i), &
              PDB%ResNum(i), &
              PDB%x(i),        &
              PDB%y(i),        &
              PDB%z(i),        & 
              PDB%blnk
          enddo
      endif
      !close(unit)
  end subroutine

  
end module
