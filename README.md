# DeepSV

Structural Variant caller based on deep learning.

This will hold the code for the various steps of deepSV, from training to varaint calling.

This will probably be constituted of multiple languages and binaries.

## Dependencies

### One-liner

```bash
curl https://sh.rustup.rs -sSf | sh
source $HOME/.cargo/env
```

### Explanations

- In order to build the `preproc` binary you will need to have `cargo` installed you can get it
  with.

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
