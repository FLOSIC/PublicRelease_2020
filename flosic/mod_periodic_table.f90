!> A Fortran Periodic Table Database
!> @author Sam Miller (https://github.com/smillerc)
!> @file mod_periodic_table.f90
!> @note Source code is taken from  https://github.com/smillerc/fortran_periodic_table
module mod_periodic_table
  use iso_fortran_env, only: ik => int32, rk => real64

  private
  public :: get_element, make_compound, element_t, compound_t

  type :: element_t
    !< A simple type to represent an element from the periodic table
    character(len=:), allocatable :: name
    character(len=:), allocatable :: symbol
    real(rk) :: mass = 0.0_rk !< Mass of the element in amu
    integer(ik) :: n_protons = 0
  end type

  type :: compound_t
    !< Simple type to represent a compound of multiple elements
    type(element_t), dimension(:), allocatable :: elements
    real(rk), dimension(:), allocatable :: element_fractions
    character(len=:), allocatable :: name
    real(rk) :: average_z = 0.0_rk
    real(rk) :: mass = 0.0_rk !< Mass of the element in amu
    integer(ik) :: n_elements = 0
  contains
    procedure, pass :: write => write_compound
    generic, public :: write(formatted) => write
  end type
contains

  type(element_t) function get_element(symbol_or_name) result(elem)
    !< Factory for creating element types
    character(len=*) :: symbol_or_name

    select case(symbol_or_name)
    case("hydrogen", "protium", "hydrogen-1", "H")
      elem = element_t(name="hydrogen", symbol="H", mass=1.00794_rk, n_protons=1)
    case("deuterium", "hydrogen-2", "2H", "D")
      elem = element_t(name="deuterium", symbol="D", mass=2.014101778_rk, n_protons=1)
    case("tritium", "hydrogen-3", "3H", "T")
      elem = element_t(name="tritium", symbol="T", mass=3.0160492675_rk, n_protons=1)
    case("helium", "He")
      elem = element_t(name="helium", symbol="He", mass=4.002602_rk, n_protons=2)
    case("helium-3", "3He", "He-3")
      elem = element_t(name="helium", symbol="He", mass=3.0160293097_rk, n_protons=2)
    case("lithium", "Li")
      elem = element_t(name="lithium", symbol="Li", mass=6.941_rk, n_protons=3)
    case("beryllium", "Be")
      elem = element_t(name="beryllium", symbol="Be", mass=9.012182_rk, n_protons=4)
    case("boron", "B")
      elem = element_t(name="boron", symbol="B", mass=10.811_rk, n_protons=5)
    case("carbon", "C")
      elem = element_t(name="carbon", symbol="C", mass=12.0107_rk, n_protons=6)
    case("nitrogen", "N")
      elem = element_t(name="nitrogen", symbol="N", mass=14.0067_rk, n_protons=7)
    case("oxygen", "O")
      elem = element_t(name="oxygen", symbol="O", mass=15.9994_rk, n_protons=8)
    case("fluorine", "F")
      elem = element_t(name="fluorine", symbol="F", mass=18.9984032_rk, n_protons=9)
    case("neon", "Ne")
      elem = element_t(name="neon", symbol="Ne", mass=20.1797_rk, n_protons=10)
    case("sodium", "Na")
      elem = element_t(name="sodium", symbol="Na", mass=22.98977_rk, n_protons=11)
    case("magnesium", "Mg")
      elem = element_t(name="magnesium", symbol="Mg", mass=24.305_rk, n_protons=12)
    case("aluminum", "Al")
      elem = element_t(name="aluminum", symbol="Al", mass=26.981538_rk, n_protons=13)
    case("silicon", "Si")
      elem = element_t(name="silicon", symbol="Si", mass=28.0855_rk, n_protons=14)
    case("phosphorus", "P")
      elem = element_t(name="phosphorus", symbol="P", mass=30.973761_rk, n_protons=15)
    case("sulfur", "S")
      elem = element_t(name="sulfur", symbol="S", mass=32.065_rk, n_protons=16)
    case("chlorine", "Cl")
      elem = element_t(name="chlorine", symbol="Cl", mass=35.453_rk, n_protons=17)
    case("argon", "Ar")
      elem = element_t(name="argon", symbol="Ar", mass=39.948_rk, n_protons=18)
    case("potassium", "K")
      elem = element_t(name="potassium", symbol="K", mass=39.0983_rk, n_protons=19)
    case("calcium", "Ca")
      elem = element_t(name="calcium", symbol="Ca", mass=40.078_rk, n_protons=20)
    case("scandium", "Sc")
      elem = element_t(name="scandium", symbol="Sc", mass=44.95591_rk, n_protons=21)
    case("titanium", "Ti")
      elem = element_t(name="titanium", symbol="Ti", mass=47.867_rk, n_protons=22)
    case("vanadium", "V")
      elem = element_t(name="vanadium", symbol="V", mass=50.9415_rk, n_protons=23)
    case("chromium", "Cr")
      elem = element_t(name="chromium", symbol="Cr", mass=51.9961_rk, n_protons=24)
    case("manganese", "Mn")
      elem = element_t(name="manganese", symbol="Mn", mass=54.938049_rk, n_protons=25)
    case("iron", "Fe")
      elem = element_t(name="iron", symbol="Fe", mass=55.845_rk, n_protons=26)
    case("cobalt", "Co")
      elem = element_t(name="cobalt", symbol="Co", mass=58.9332_rk, n_protons=27)
    case("nickel", "Ni")
      elem = element_t(name="nickel", symbol="Ni", mass=58.6934_rk, n_protons=28)
    case("copper", "Cu")
      elem = element_t(name="copper", symbol="Cu", mass=63.546_rk, n_protons=29)
    case("zinc", "Zn")
      elem = element_t(name="zinc", symbol="Zn", mass=65.409_rk, n_protons=30)
    case("gallium", "Ga")
      elem = element_t(name="gallium", symbol="Ga", mass=69.723_rk, n_protons=31)
    case("germanium", "Ge")
      elem = element_t(name="germanium", symbol="Ge", mass=72.64_rk, n_protons=32)
    case("arsenic", "As")
      elem = element_t(name="arsenic", symbol="As", mass=74.9216_rk, n_protons=33)
    case("selenium", "Se")
      elem = element_t(name="selenium", symbol="Se", mass=78.96_rk, n_protons=34)
    case("bromine", "Br")
      elem = element_t(name="bromine", symbol="Br", mass=79.904_rk, n_protons=35)
    case("krypton", "Kr")
      elem = element_t(name="krypton", symbol="Kr", mass=83.798_rk, n_protons=36)
    case("rubidium", "Rb")
      elem = element_t(name="rubidium", symbol="Rb", mass=85.4678_rk, n_protons=37)
    case("strontium", "Sr")
      elem = element_t(name="strontium", symbol="Sr", mass=87.62_rk, n_protons=38)
    case("yttrium", "Y")
      elem = element_t(name="yttrium", symbol="Y", mass=88.90585_rk, n_protons=39)
    case("zirconium", "Zr")
      elem = element_t(name="zirconium", symbol="Zr", mass=91.224_rk, n_protons=40)
    case("niobium", "Nb")
      elem = element_t(name="niobium", symbol="Nb", mass=92.90638_rk, n_protons=41)
    case("molybdenum", "Mo")
      elem = element_t(name="molybdenum", symbol="Mo", mass=95.94_rk, n_protons=42)
    case("technetium", "Tc")
      elem = element_t(name="technetium", symbol="Tc", mass=98_rk, n_protons=43)
    case("ruthenium", "Ru")
      elem = element_t(name="ruthenium", symbol="Ru", mass=101.07_rk, n_protons=44)
    case("rhodium", "Rh")
      elem = element_t(name="rhodium", symbol="Rh", mass=102.9055_rk, n_protons=45)
    case("palladium", "Pd")
      elem = element_t(name="palladium", symbol="Pd", mass=106.42_rk, n_protons=46)
    case("silver", "Ag")
      elem = element_t(name="silver", symbol="Ag", mass=107.8682_rk, n_protons=47)
    case("cadmium", "Cd")
      elem = element_t(name="cadmium", symbol="Cd", mass=112.411_rk, n_protons=48)
    case("indium", "In")
      elem = element_t(name="indium", symbol="In", mass=114.818_rk, n_protons=49)
    case("tin", "Sn")
      elem = element_t(name="tin", symbol="Sn", mass=118.71_rk, n_protons=50)
    case("antimony", "Sb")
      elem = element_t(name="antimony", symbol="Sb", mass=121.76_rk, n_protons=51)
    case("tellurium", "Te")
      elem = element_t(name="tellurium", symbol="Te", mass=127.6_rk, n_protons=52)
    case("iodine", "I")
      elem = element_t(name="iodine", symbol="I", mass=126.90447_rk, n_protons=53)
    case("xenon", "Xe")
      elem = element_t(name="xenon", symbol="Xe", mass=131.293_rk, n_protons=54)
    case("cesium", "Cs")
      elem = element_t(name="cesium", symbol="Cs", mass=132.90545_rk, n_protons=55)
    case("barium", "Ba")
      elem = element_t(name="barium", symbol="Ba", mass=137.327_rk, n_protons=56)
    case("lanthanum", "La")
      elem = element_t(name="lanthanum", symbol="La", mass=138.9055_rk, n_protons=57)
    case("cerium", "Ce")
      elem = element_t(name="cerium", symbol="Ce", mass=140.116_rk, n_protons=58)
    case("praseodymium", "Pr")
      elem = element_t(name="praseodymium", symbol="Pr", mass=140.90765_rk, n_protons=59)
    case("neodymium", "Nd")
      elem = element_t(name="neodymium", symbol="Nd", mass=144.24_rk, n_protons=60)
    case("promethium", "Pm")
      elem = element_t(name="promethium", symbol="Pm", mass=145_rk, n_protons=61)
    case("samarium", "Sm")
      elem = element_t(name="samarium", symbol="Sm", mass=150.36_rk, n_protons=62)
    case("europium", "Eu")
      elem = element_t(name="europium", symbol="Eu", mass=151.964_rk, n_protons=63)
    case("gadolinium", "Gd")
      elem = element_t(name="gadolinium", symbol="Gd", mass=157.25_rk, n_protons=64)
    case("terbium", "Tb")
      elem = element_t(name="terbium", symbol="Tb", mass=158.92534_rk, n_protons=65)
    case("dysprosium", "Dy")
      elem = element_t(name="dysprosium", symbol="Dy", mass=162.5_rk, n_protons=66)
    case("holmium", "Ho")
      elem = element_t(name="holmium", symbol="Ho", mass=164.93032_rk, n_protons=67)
    case("erbium", "Er")
      elem = element_t(name="erbium", symbol="Er", mass=167.259_rk, n_protons=68)
    case("thulium", "Tm")
      elem = element_t(name="thulium", symbol="Tm", mass=168.93421_rk, n_protons=69)
    case("ytterbium", "Yb")
      elem = element_t(name="ytterbium", symbol="Yb", mass=173.04_rk, n_protons=70)
    case("lutetium", "Lu")
      elem = element_t(name="lutetium", symbol="Lu", mass=174.967_rk, n_protons=71)
    case("hafnium", "Hf")
      elem = element_t(name="hafnium", symbol="Hf", mass=178.49_rk, n_protons=72)
    case("tantalum", "Ta")
      elem = element_t(name="tantalum", symbol="Ta", mass=180.9479_rk, n_protons=73)
    case("tungsten", "W")
      elem = element_t(name="tungsten", symbol="W", mass=183.84_rk, n_protons=74)
    case("rhenium", "Re")
      elem = element_t(name="rhenium", symbol="Re", mass=186.207_rk, n_protons=75)
    case("osmium", "Os")
      elem = element_t(name="osmium", symbol="Os", mass=190.23_rk, n_protons=76)
    case("iridium", "Ir")
      elem = element_t(name="iridium", symbol="Ir", mass=192.217_rk, n_protons=77)
    case("platinum", "Pt")
      elem = element_t(name="platinum", symbol="Pt", mass=195.078_rk, n_protons=78)
    case("gold", "Au")
      elem = element_t(name="gold", symbol="Au", mass=196.96655_rk, n_protons=79)
    case("mercury", "Hg")
      elem = element_t(name="mercury", symbol="Hg", mass=200.59_rk, n_protons=80)
    case("thallium", "Tl")
      elem = element_t(name="thallium", symbol="Tl", mass=204.3833_rk, n_protons=81)
    case("lead", "Pb")
      elem = element_t(name="lead", symbol="Pb", mass=207.2_rk, n_protons=82)
    case("bismuth", "Bi")
      elem = element_t(name="bismuth", symbol="Bi", mass=208.98038_rk, n_protons=83)
    case("polonium", "Po")
      elem = element_t(name="polonium", symbol="Po", mass=209_rk, n_protons=84)
    case("astatine", "At")
      elem = element_t(name="astatine", symbol="At", mass=210_rk, n_protons=85)
    case("radon", "Rn")
      elem = element_t(name="radon", symbol="Rn", mass=222_rk, n_protons=86)
    case("francium", "Fr")
      elem = element_t(name="francium", symbol="Fr", mass=223_rk, n_protons=87)
    case("radium", "Ra")
      elem = element_t(name="radium", symbol="Ra", mass=226_rk, n_protons=88)
    case("actinium", "Ac")
      elem = element_t(name="actinium", symbol="Ac", mass=227_rk, n_protons=89)
    case("thorium", "Th")
      elem = element_t(name="thorium", symbol="Th", mass=232.0381_rk, n_protons=90)
    case("protactinium", "Pa")
      elem = element_t(name="protactinium", symbol="Pa", mass=231.03588_rk, n_protons=91)
    case("uranium", "U")
      elem = element_t(name="uranium", symbol="U", mass=238.02891_rk, n_protons=92)
    case("neptunium", "Np")
      elem = element_t(name="neptunium", symbol="Np", mass=237_rk, n_protons=93)
    case("plutonium", "Pu")
      elem = element_t(name="plutonium", symbol="Pu", mass=244_rk, n_protons=94)
    case("americium", "Am")
      elem = element_t(name="americium", symbol="Am", mass=243_rk, n_protons=95)
    case("curium", "Cm")
      elem = element_t(name="curium", symbol="Cm", mass=247_rk, n_protons=96)
    case("berkelium", "Bk")
      elem = element_t(name="berkelium", symbol="Bk", mass=247_rk, n_protons=97)
    case("californium", "Cf")
      elem = element_t(name="californium", symbol="Cf", mass=251_rk, n_protons=98)
    case("einsteinium", "Es")
      elem = element_t(name="einsteinium", symbol="Es", mass=252_rk, n_protons=99)
    case("fermium", "Fm")
      elem = element_t(name="fermium", symbol="Fm", mass=257_rk, n_protons=100)
    case("mendelevium", "Md")
      elem = element_t(name="mendelevium", symbol="Md", mass=258_rk, n_protons=101)
    case("nobelium", "No")
      elem = element_t(name="nobelium", symbol="No", mass=259_rk, n_protons=102)
    case("lawrencium", "Lr")
      elem = element_t(name="lawrencium", symbol="Lr", mass=262_rk, n_protons=103)
    case("rutherfordium", "Rf")
      elem = element_t(name="rutherfordium", symbol="Rf", mass=261_rk, n_protons=104)
    case("dubnium", "Db")
      elem = element_t(name="dubnium", symbol="Db", mass=262_rk, n_protons=105)
    case("seaborgium", "Sg")
      elem = element_t(name="seaborgium", symbol="Sg", mass=266_rk, n_protons=106)
    case("bohrium", "Bh")
      elem = element_t(name="bohrium", symbol="Bh", mass=264_rk, n_protons=107)
    case("hassium", "Hs")
      elem = element_t(name="hassium", symbol="Hs", mass=277_rk, n_protons=108)
    case("meitnerium", "Mt")
      elem = element_t(name="meitnerium", symbol="Mt", mass=268_rk, n_protons=109)
    case("darmstadtium", "Ds")
      elem = element_t(name="darmstadtium", symbol="Ds", mass=281_rk, n_protons=110)
    case("roentgenium", "Rg")
      elem = element_t(name="roentgenium", symbol="Rg", mass=272_rk, n_protons=111)
    case("copernicium", "Cn")
      elem = element_t(name="copernicium", symbol="Cn", mass=285_rk, n_protons=112)
    case("nihonium", "Nh")
      elem = element_t(name="nihonium", symbol="Nh", mass=286_rk, n_protons=113)
    case("flerovium", "Fl")
      elem = element_t(name="flerovium", symbol="Fl", mass=289_rk, n_protons=114)
    case("moscovium", "Mc")
      elem = element_t(name="moscovium", symbol="Mc", mass=289_rk, n_protons=115)
    case("livermorium", "Lv")
      elem = element_t(name="livermorium", symbol="Lv", mass=293_rk, n_protons=116)
    case("tennessine", "Ts")
      elem = element_t(name="tennessine", symbol="Ts", mass=294_rk, n_protons=117)
    case("oganesson", "Og")
      elem = element_t(name="oganesson", symbol="Og", mass=294_rk, n_protons=118)
    case default
      error stop "Unknown element or symbol"
    end select

  end function

  type(compound_t) function make_compound(element_symbols, element_fractions) result(compound)
    !< Create a compound made of multiple elements
    character(len=*), dimension(:) :: element_symbols
    real(rk), dimension(:) :: element_fractions

    type(element_t), dimension(size(element_symbols)) :: elements
    character(len=50) :: compound_name
    real(rk) :: ave_z, mass
    integer(ik) :: i

    compound_name = ''
    if(size(element_symbols) /= size(element_fractions)) then
      error stop 'size(element_symbols) /= size(element_fractions)'
    end if

    if(abs(sum(element_fractions) - 1.0_rk) > epsilon(1.0_rk)) then
      error stop "sum(element_fractions) /= 1.0_rk"
    end if

    mass = 0.0_rk
    ave_z = 0.0_rk
    do i = 1, size(element_symbols)
      elements(i) = get_element(trim(element_symbols(i)))
      compound_name(i:i + len(elements(i)%symbol)) = elements(i)%symbol
      mass = elements(i)%mass * element_fractions(i)
      ave_z = ave_z + elements(i)%n_protons
    end do
    ave_z = ave_z / size(element_symbols)
    compound = compound_t(average_z=ave_z, mass=mass, name=compound_name, &
                          element_fractions=element_fractions, &
                          elements=elements, n_elements=size(elements))

  end function

  subroutine write_compound(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of `write(*,*) compound_t`

    class(compound_t), intent(in) :: self
    integer, intent(in) :: unit           !< input/output unit
    character(*), intent(in) :: iotype    !< LISTDIRECTED or DTxxx
    integer, intent(in) :: v_list(:)      !< parameters from fmt spec.
    integer, intent(out) :: iostat        !< non zero on error, etc.
    character(*), intent(inout) :: iomsg  !< define if iostat non zero.
    integer(ik) :: i

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "Name: "//self%name//new_line('a')

    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "Total mass [amu]: ", self%mass, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "Average Z: ", self%average_z, new_line('a')

    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) new_line('a')//"Elements:"//new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "========="//new_line('a')
    do i = 1, self%n_elements
      write(unit, '(a, f6.4, a, f6.4, a)', iostat=iostat, iomsg=iomsg) self%elements(i)%symbol// &
        ": fraction = ", self%element_fractions(i), &
        ", mass [amu] = ", self%element_fractions(i) * self%elements(i)%mass, new_line('a')
    end do
  end subroutine

  function lowercase(s1) result(s2)
    !< Convert a string to lowercase
    character(len=*):: s1
    character(len(s1)) :: s2
    character(len=1) :: ch
    integer(ik), parameter :: duc = ichar('a') - ichar('a')
    integer(ik) :: i

    do i = 1, len(s1)
      ch = s1(i:i)
      if(ch >= 'a' .and. ch <= 'z') ch = char(ichar(ch) - duc)
      s2(i:i) = ch
    end do
  end function lowercase

end module
