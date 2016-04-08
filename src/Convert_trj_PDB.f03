module Convert_trj_PDB 
  implicit none

  public  read_pdb
  public  readtrj
  private outputPDB

  !***Should I convert the follow variables to type ones?
  integer                       :: iatom, n_atoms, n_residues, n_chains
  integer         , allocatable :: AtomNum(:)
  character(len=4), allocatable :: AtomName(:)
  character(len=3), allocatable :: ResName(:)
  character(len=1), allocatable :: ChainId(:)
  integer         , allocatable :: ResNum(:)
  double precision, allocatable :: x(:),y(:),z(:)
  character(len=126) :: blnk=""

  real(4) :: potential
  real(4), allocatable   :: trj_x(:), trj_y(:), trj_z(:)
  integer :: iconf, n_confs
  real(8) :: cellsize(3)
contains

  subroutine read_pdb(FilName)
    character(len=*), intent(in) :: FilName
    character(len=6)   :: atom
    integer, parameter :: FilUnit = 11
    integer, parameter :: OutPDBUnit  = 12
    integer :: i, ios, ter 
    integer :: NotAtomLine,NumAtom 
  
    open(FilUnit, file=FilName, status="old")
    print*,"  File Name: ",FilName
    allocate(AtomNum(n_atoms),AtomName(n_atoms),ResName(n_atoms), &
             ChainId(n_atoms),ResNum(n_atoms),x(n_atoms),y(n_atoms),z(n_atoms))
    AtomNum(:) = 0
    ResNum(:)  = 0
    x(:) = 0; y(:) = 0; z(:) = 0
  
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
    do i = 1, n_atoms 
      read(FilUnit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") &
           atom       , &
           AtomNum(i) , &
           AtomName(i), &
           ResName(i) , &
           ChainId(i) , &
           ResNum(i)  , &
           x(i),        &
           y(i),        &
           z(i)  
    enddo
    close(FilUnit)
  end subroutine

!------------------------------------------------------
  subroutine readtrj(filename)
      use apply_PBC 
      implicit none
      character(len=*) :: filename
      integer    :: unit = 11, unit_outPDB = 10
      integer(4) :: istp,iyn15v,iyn15h
      real(4)    :: sitime, sec, et, kinetic, temperature, rmsf, rmsd 

  
      open(unit, file = filename, form="unformatted", status="old")
      allocate(trj_x(n_atoms), trj_y(n_atoms), trj_z(n_atoms))
      
      !***Initialization
      iconf = 0; trj_x(:) = 0; trj_y(:) = 0; trj_z(:) = 0
      
      !***Reading trajectory
      call prepare_apply_PBC(n_atoms,ResNum,n_residues)
      do 
          read(unit, end=111) istp,sitime,sec,et,kinetic,temperature,&
                              potential,rmsf,iyn15v,iyn15h,rmsd
          read(unit) (trj_x(iatom), trj_y(iatom), trj_z(iatom), iatom = 1, n_atoms)
          iconf = iconf + 1
          call ReturnAtom(n_atoms,n_chains,n_residues,trj_x, trj_y, trj_z, cellsize)

          call outputPDB(unit_outPDB)
          print*,"Conf NO:", iconf 
      enddo
      111 close(unit_outPDB)
      n_confs = iconf
      write(*,'("#Number of conformation ",i8)') n_confs 
  end subroutine

  subroutine outputPDB(unit)
      integer, intent(in) :: unit 
      open(unit, file ="frames.pdb")!, status="replace")
      write(unit,"('MODEL', i7)") iconf
      write(unit,"('#Potential, kcal/mol ', f10.3)") potential 
      do iatom = 1, n_atoms 
        write(unit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
          "ATOM  ", &
          AtomNum(iatom), &
          AtomName(iatom),&
          ResName(iatom),&
          ChainId(iatom), &
          ResNum(iatom), &
          trj_x(iatom), trj_y(iatom), trj_z(iatom),&
          blnk
      enddo
      write(unit,"(a6)") "ENDMDL"
      !close(unit)
  end subroutine

  
end module
