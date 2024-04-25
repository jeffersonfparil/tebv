{ pkgs ? import <nixpkgs> { } }:
with pkgs;
let
  simquantgen = pkgs.rPackages.buildRPackage {
    name = "simquantgen";
    src = pkgs.fetchFromGitHub {
      owner = "jeffersonfparil";
      repo = "simquantgen";
      rev = "bd5f07b13cb2ae0aecfec4d84c152bd664843f4d";
      sha256 = "hp3uXuc7LUTIMIVz8h2NF/C4PLTpSfVFjnDyNt5eDF8=";
    };
    propagatedBuildInputs = with rPackages; [
        MASS
        doParallel
        foreach
        txtplot
    ];
  }; ### Use lib.fakeSha256 on "sha256" to get the correct code from the error messages
  my-r-pkgs = rWrapper.override {
    packages = with rPackages; [
      remotes
      devtools
      testthat
      txtplot
      sommer
      MASS
      doParallel
      foreach
      simquantgen
    ];
  };
in mkShell {
  buildInputs = with pkgs; [git glibcLocales openssl which openssh curl wget ];
  inputsFrom = [ my-r-pkgs ];
  shellHook = ''
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
  '';
  GIT_SSL_CAINFO = "${cacert}/etc/ssl/certs/ca-bundle.crt";
#  LOCALE_ARCHIVE = stdenv.lib.optionalString stdenv.isLinux
#    "${glibcLocales}/lib/locale/locale-archive";
}

### nix-shell --run bash --pure