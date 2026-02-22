use std::path::PathBuf;

use jarzynski_rs::{compute_bins, read_pull_files};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: jarzynski-cli <filename-prefix> <pull-directory> <number-of-files>");
        std::process::exit(1);
    }

    let prefix = &args[1];
    let pull_dir = PathBuf::from(&args[2]);
    let number_of_files: usize = args[3].parse().expect("number-of-files must be an integer");

    let paths: Vec<PathBuf> = (1..=number_of_files)
        .map(|idx| pull_dir.join(format!("{prefix}.{idx}")))
        .collect();

    let samples = match read_pull_files(&paths) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Error reading pull files: {e}");
            std::process::exit(1);
        }
    };

    let bins = match compute_bins(&samples, 303.0) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("Error computing free-energy bins: {e}");
            std::process::exit(1);
        }
    };

    println!("center lower upper raw raw_err taylor taylor_err alpha alpha_err");
    for bin in bins {
        println!(
            "{} {:.3} {:.3} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
            bin.center,
            bin.lower,
            bin.upper,
            bin.raw.value,
            bin.raw.stdev,
            bin.taylor.value,
            bin.taylor.stdev,
            bin.alpha.value,
            bin.alpha.stdev,
        );
    }
}
