module variables

    type var_PDB
        integer                       :: iatom, n_atoms, n_residues, n_chains
        integer         , allocatable :: AtomNum(:)
        character(len=4), allocatable :: AtomName(:)
        character(len=3), allocatable :: ResName(:)
        character(len=1), allocatable :: ChainId(:)
        integer         , allocatable :: ResNum(:)
        double precision, allocatable :: x(:),y(:),z(:)
        character(len=126) :: blnk=""
    end type

    type var_trj
        real(4) :: potential
        real(4), allocatable   :: x(:), y(:), z(:)
        integer :: n_confs
        real(8) :: cellsize(3)
    end type


end module
