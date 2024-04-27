{ pkgs ? import <nixpkgs> { } }:
with pkgs;
mkShell {
  buildInputs = [
    stdenv
    valgrind
    gdb
    linuxPackages_latest.perf
  ];

  shellHook = ''
    # ...
  '';
}
