with import <nixpkgs> {}; {
  pcmEnv = stdenv.mkDerivation {
    name = "PCMSolver";
    buildInputs = [
      boost
      bundler
      ccache
      clang
      clang-analyzer
      clang-tools
      cmake
      doxygen
      gcc
      gdb
      gfortran
      lldb
      openmpi
      python35Packages.breathe
      python35Packages.jupyter
      python35Packages.matplotlib
      python35Packages.numpy
      python35Packages.pyyaml
      python35Packages.recommonmark
      python35Packages.scipy
      python35Packages.sphinx
      python35Packages.sphinx_rtd_theme
      python35Packages.sympy
      valgrind
      zlib
    ];
  };
}
