use clap::Parser;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
struct Args {
    #[clap(short = 'i', long = "input", required = true, num_args = 1..)]
    input_files: Vec<PathBuf>,
}

fn generate_output_path(input_path: &Path) -> PathBuf {
    let stem = input_path.file_stem().unwrap_or_default().to_string_lossy();
    let ext = input_path
        .extension()
        .map(|e| e.to_string_lossy())
        .unwrap_or_default();

    let new_filename = if ext.is_empty() {
        format!("{}_cleaned", stem)
    } else {
        format!("{}_cleaned.{}", stem, ext)
    };

    input_path
        .parent()
        .unwrap_or(Path::new(""))
        .join(new_filename)
}

fn generate_csv_output_path(input_path: &Path) -> PathBuf {
    let stem = input_path.file_stem().unwrap_or_default().to_string_lossy();
    let new_filename = format!("{}_cleaned.csv", stem);
    input_path
        .parent()
        .unwrap_or(Path::new(""))
        .join(new_filename)
}

fn process_file(
    input_path: &PathBuf,
    output_path: &PathBuf,
    csv_output_path: &PathBuf,
) -> io::Result<()> {
    let input_file = File::open(input_path)?;
    let mut reader = std::io::BufReader::new(input_file);

    let output_file = File::create(output_path)?;
    let mut writer = std::io::BufWriter::new(output_file);

    let csv_file = File::create(csv_output_path)?;
    let mut csv_writer = std::io::BufWriter::new(csv_file);

    let mut line = String::new();
    loop {
        let length = reader.read_line(&mut line)?;
        if length == 0 {
            break;
        }

        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() == 9 {
            // Check if it's the header line
            if parts[0].trim() == "Freq(Hz)" {
                // Process header
                let new_header = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    parts[0].trim(),
                    parts[4].trim(),
                    parts[5].trim(),
                    parts[6].trim(),
                    parts[7].trim(),
                    parts[8].trim()
                );
                writer.write_all(new_header.as_bytes())?;
                writer.write_all(b"\n")?;

                // // Write CSV header
                // // For DRTtool compatibility, we comment this out, and give only data rows.
                // let csv_header = format!(
                //     "{},{},{},{},{},{}",
                //     parts[0].trim(),
                //     parts[4].trim(),
                //     parts[5].trim(),
                //     parts[6].trim(),
                //     parts[7].trim(),
                //     parts[8].trim()
                // );
                // csv_writer.write_all(csv_header.as_bytes())?;
                // csv_writer.write_all(b"\n")?;
            } else {
                let izr: Result<f64, _> = parts[5].parse();
                if izr.unwrap_or(1.0) <= 0.0 {
                    let new_data = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        parts[0].trim(),
                        parts[4].trim(),
                        parts[5].trim(),
                        parts[6].trim(),
                        parts[7].trim(),
                        parts[8].trim()
                    );
                    writer.write_all(new_data.as_bytes())?;
                    writer.write_all(b"\n")?;

                    // Write CSV data row
                    let csv_data = format!(
                        "{},{},{},{},{},{}",
                        parts[0].trim(),
                        parts[4].trim(),
                        parts[5].trim(),
                        parts[6].trim(),
                        parts[7].trim(),
                        parts[8].trim()
                    );
                    csv_writer.write_all(csv_data.as_bytes())?;
                    csv_writer.write_all(b"\n")?;
                }
                // else {
                //     println!("line jumped: {:?}", line);
                // }
            }
        } else {
            // Write other lines as-is (only to the original output file)
            writer.write_all(line.as_bytes())?;
        }

        line.clear();
    }

    Ok(())
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    for input_path in &args.input_files {
        let output_path = generate_output_path(input_path);
        let csv_output_path = generate_csv_output_path(input_path);
        process_file(input_path, &output_path, &csv_output_path)?;
        println!(
            "✅ 处理完成: {} → {} and {}",
            input_path.display(),
            output_path.display(),
            csv_output_path.display()
        );
    }

    Ok(())
}
