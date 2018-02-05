# DeepSV
[![Build Status](https://travis-ci.org/gamazeps/deepSV.svg?branch=master)](https://travis-ci.org/gamazeps/deepSV)

Structural Variant caller based on deep learning.

This will hold the code for the various steps of deepSV, from training to varaint calling.

This will probably be constituted of multiple languages and binaries.

## Dependencies

### One-liner

A working installation of rust is needed for running the data generation.

```bash
curl https://sh.rustup.rs -sSf | sh
source $HOME/.cargo/env
```

Note thqt you also need to have `java`, `picardtools` and `samtools` installed for the many
calls to bioinformatics tool (before we reimplement them in rust of course).

## Build

All instructions assume that they are executed from the root of the repository.

### Preproc

```bash
cd preproc
cargo build
```

## Guarantees

This is a *very* work in progress, do not touch it for now !

This obviously comes with no guarantees at all, it may burn your computer or give you cancer for all
I know.
