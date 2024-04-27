{ pkgs ? import <nixpkgs> { } }:
with pkgs;
mkShell {
  buildInputs = [
    gcc
    gnumake
    valgrind
    gdb
    linuxPackages_latest.perf
    firefox
  ];

  shellHook = ''
    # ...
  '';
}
