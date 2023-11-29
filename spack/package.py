from spack import *

class Flosic(MakefilePackage):
    """FLOSIC is a software package for performing Fermi-Löwdin orbital self-interaction corrections."""

    homepage = "https://github.com/FLOSIC/PublicRelease_2020"
    url      = "https://github.com/FLOSIC/PublicRelease_2020/archive/refs/tags/v1.0.tar.gz"
    git      = "https://github.com/FLOSIC/PublicRelease_2020.git"

    version('1.0', sha256='b98973d29a556ab1c91e4ee1e83f11c465d5d8d57cbd37f445620e00fe0bef5a')

    variant('mpi', default=True, description='Enable MPI support')
    variant('group', default=False, description='Enable group calculation')
    variant('atomforce', default=False, description='Enable atomic force in FLOSIC calculation')

    depends_on('openblas')
    depends_on('mpi', when='+mpi')

    def edit(self, spec, prefix):
       # Add linker flags to the Makefile before compiling
        makefile_path = join_path('flosic', 'Makefile.fedora')
        openblas_flags = self.spec['openblas'].libs
        filter_file('-llapack -lblas', ' '.join(openblas_flags), makefile_path)

        makefile = FileFilter('flosic/Makefile.fedora')

        # Adjust include and library directories
        include_path = spec['mpi'].prefix.include
        lib_path = spec['openblas'].prefix.lib

        # Adjust the Makefile to include the correct paths for MPI and OpenBLAS
        makefile.filter('^IDIR = .*', 'IDIR = -I' + include_path)
        makefile.filter('^LDIR = .*', 'LDIR = -L' + lib_path)

        # Link against the OpenBLAS library
        openblas_flags = self.spec['openblas'].libs
        makefile.filter('LIBS = -llapack -lblas', 'LIBS = {}'.format(openblas_flags))

        with working_dir('flosic'):
            makefile = FileFilter('Makefile.fedora')

            # Set MPI option
            if '+mpi' in spec:
                makefile.filter('MPI=N', 'MPI=Y')
            else:
                makefile.filter('MPI=Y', 'MPI=N')

            # Add any additional filters or changes needed here

            # Set group option
            if '+group' in spec:
                makefile.filter('GROUP=N', 'GROUP=Y')
            else:
                makefile.filter('GROUP=Y', 'GROUP=N')

            # Set atomforce option
            if '+atomforce' in spec:
                makefile.filter('ATOMFORCE=N', 'ATOMFORCE=Y')
            else:
                makefile.filter('ATOMFORCE=Y', 'ATOMFORCE=N')

            # Rename Makefile.fedora to Makefile
            move('Makefile.fedora', 'Makefile')

            # Adjust compilers and flags if needed
            # Example:
            # makefile.filter('CC = gcc', 'CC = {}'.format(spack_cc))
            # makefile.filter('FC = mpif90', 'FC = {}'.format(spack_fc))

    def build(self, spec, prefix):
        with working_dir('flosic'):
            fc = Executable(spack_fc)

            # Compile modules that generate .mod files first
            fc('-o', 'condcomp', 'condcomp.f')

            make(parallel=False)
 

    def install(self, spec, prefix):
        # Continue with the rest of the installation steps
        mkdirp(prefix.bin)
        install('flosic/nrlmol_exe', prefix.bin)

        # Handle other installation steps as needed
 
