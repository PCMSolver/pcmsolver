let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "nixos-unstable";
    sha256 = "1ig27lcx533jjdxg9f7wj7n9bwkd12yj5zpzkn4qx60kpbvla27c";
  });
in
  with import nixpkgs {
    overlays = [(self: super:
      {
        python3 = super.python3.override {
          packageOverrides = py-self: py-super: {
            matplotlib = py-super.matplotlib.override {
              enableTk = true;
              enableQt = true;
            };
          };
        };
      }
    )];
  };
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
      ffmpeg
      gfortran
      graphviz
      lcov
      pipenv
      python3Full
      python3Packages.docopt
      python3Packages.jupyter
      python3Packages.matplotlib
      python3Packages.numpy
      python3Packages.pyyaml
      valgrind
      zlib
    ];
    src = null;
    shellHook = ''
    export NINJA_STATUS="[Built edge %f of %t in %e sec]"
    SOURCE_DATE_EPOCH=$(date +%s)
    '';
  }
