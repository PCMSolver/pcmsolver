let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "nixos-unstable";
    sha256 = "0q47hjl01wylp6cdb52fyjbbx45djnz7bwfyj2mfv2lh04iph3s8";
  });
in
  with import nixpkgs {
    overlays = [(self: super:
      # Within the overlay we use a recursive set, though I think we can use `self` as well.
      rec {
        # nix-shell -p python.pkgs.my_stuff
        python3 = super.python3.override {
          # Careful, we're using a different self and super here!
          packageOverrides = py-self: py-super: {
            matplotlib = py-super.matplotlib.override {
              enableTk = true;
              enableQt = true;
            };
          };
        };
        python3Packages = python3.pkgs;
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
      lldb
      pipenv
      python3Full
      python3Packages.docopt
      python3Packages.jupyter
      python3Packages.matplotlib
      python3Packages.numpy
      python3Packages.pyyaml
      python3Packages.virtualenvwrapper
      valgrind
      zlib
    ];
    src = null;
    shellHook = ''
      NINJA_STATUS="[Built edge %f of %t in %e sec]"
      SOURCE_DATE_EPOCH=$(date +%s)
    '';
  }
