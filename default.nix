let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "nixos-unstable";
    sha256 = "1d4q92jw42d51s9bn380jayy2rs1v0h1y8kvrbjg3i43f72ck5q5";
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
      gfortran
      graphviz
      lldb
      python3Packages.docopt
      python3Packages.matplotlib
      python3Packages.numpy
      python3Packages.pyyaml
      python3Packages.virtualenvwrapper
      valgrind
      zlib
    ];
  }
