let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "ac355040656de04f59406ba2380a96f4124ebdad";
    sha256 = "0frhc7mnx88sird6ipp6578k5badibsl0jfa22ab9w6qrb88j825";
  });
in
  with import nixpkgs {};
  stdenv.mkDerivation {
    name = "PCMSolver-dev";
    buildInputs = [
      boost
      bundler
      ccache
      clang
      clang-analyzer
      clang-tools
      cmake
      doxygen
      exa
      #gcc
      gdb
      gfortran
      graphviz
      lldb
      python35Packages.matplotlib
      python35Packages.numpy
      python35Packages.virtualenvwrapper
      valgrind
      zlib
    ];
  }
