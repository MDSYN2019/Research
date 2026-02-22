# Rust Jarzynski port

This crate ports the serial Jarzynski processing workflow to Rust so it can be used as:

- a reusable Rust library (`jarzynski_rs`), and
- a command-line tool (`jarzynski-cli`) with the same pull-file pattern used by the C++ binary.

## Build and test

```bash
cargo test
```

## Run

```bash
cargo run --bin jarzynski-cli -- pull ./../data 4
```

Arguments:

1. file prefix (`pull`)
2. directory containing pull files
3. number of files
