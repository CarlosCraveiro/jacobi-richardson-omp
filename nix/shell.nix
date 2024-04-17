{ pkgs ? import <nixpkgs> { } }:
with pkgs;
mkShell {
  buildInputs = [
    gcc
    gnumake
  ];

  shellHook = ''
    # ...
  '';
}
