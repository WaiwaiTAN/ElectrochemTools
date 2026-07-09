use clap::Parser;
use regex::Regex;
use std::io::{self, Write};
use std::{fs, fs::File, path::Path, path::PathBuf};
use regex::Regex;

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
        format!("{}_trimmed", stem)
    } else {
        format!("{}_trimmed.{}", stem, ext)
    };

    input_path
        .parent()
        .unwrap_or(Path::new(""))
        .join(new_filename)
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    for input_path in &args.input_files {
        let output_path = generate_output_path(input_path);
        process_file(input_path, &output_path)?;
        println!(
            "✅ 处理完成: {} → {}",
            input_path.display(),
            output_path.display()
        );
    }

    Ok(())
}

fn process_file(input_path: &PathBuf, output_path: &PathBuf) -> io::Result<()> {
    // 读取文件内容
    let content = fs::read_to_string(input_path)?;
    let lines: Vec<&str> = content.lines().collect();

    // 解析实验参数
    let params = parse_parameters(&lines);
    println!("{:?}", params);

    let n_keep = calculate_points_to_keep(&params);

    // 定位数据起始位置
    let data_start_index = lines
        .iter()
        .position(|&line| line.starts_with("End Comments"))
        .map(|idx| idx + 1)
        .unwrap_or(0);

    // 提取并验证数据
    let (trimmed_data, warnings) =
        extract_and_validate_data(&lines[data_start_index..], n_keep, &params);

    // 输出警告信息
    if !warnings.is_empty() {
        eprintln!("⚠️ 数据完整性警告:");
        warnings.iter().for_each(|msg| eprintln!("  - {msg}"));
    }

    // 构建输出文件
    let mut output = File::create(output_path)?;
    write!(output, "{}", lines[..data_start_index].join("\n"))?;
    write!(output, "\n{}", trimmed_data.join("\n"))?;

    Ok(())
}

fn is_close(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance
}

fn extract_and_validate_data(
    data_lines: &[&str],
    n_keep: usize,
    params: &ExperimentParams,
) -> (Vec<String>, Vec<String>) {
    let total_points = data_lines.len();
    let start_index = if total_points > n_keep {
        total_points - n_keep
    } else {
        0
    };
    let mut warnings = Vec::new();

    // 提取首尾电压值
    let first_voltage = parse_voltage(data_lines[start_index]);
    let last_voltage = parse_voltage(data_lines[data_lines.len() - 1]);

    // 检查中间转折点（周期中点应接近极值）
    if n_keep > 10 {
        let mid_point = start_index + n_keep / 2;
        if let Some(mid_line) = data_lines.get(mid_point) {
            let mid_voltage = parse_voltage(mid_line);
            // case1: 首位电压接近 voltage_step1，中点应接近 voltage_step2
            if is_close(first_voltage, params.voltage_step1, 0.05)
                && is_close(last_voltage, params.voltage_step1, 0.05)
            {
                if !is_close(mid_voltage, params.voltage_step2, 0.05) {
                    warnings.push(format!(
                        "中点电压 {:.2} 不接近预期的 {:.2}（基于首位电压为 voltage_step1）",
                        mid_voltage, params.voltage_step2
                    ));
                }
            }
            // case2: 首位电压接近 voltage_step2，中点应接近 voltage_step1
            else if is_close(first_voltage, params.voltage_step2, 0.05)
                && is_close(last_voltage, params.voltage_step2, 0.05)
            {
                if !is_close(mid_voltage, params.voltage_step1, 0.05) {
                    warnings.push(format!(
                        "中点电压 {:.2} 不接近预期的 {:.2}（基于首位电压为 voltage_step2）",
                        mid_voltage, params.voltage_step1
                    ));
                }
            }
            // 无法判断趋势
            else {
                warnings.push("首位电压不一致，无法判断周期趋势".to_string());
            }
        }
    }

    (
        data_lines[start_index..]
            .iter()
            .map(|&s| s.to_string())
            .collect(),
        warnings,
    )
}

// 辅助函数：从数据行解析电压值
fn parse_voltage(line: &str) -> f64 {
    line.split('\t')
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(f64::NAN)
}

#[derive(Debug)]
struct ExperimentParams {
    scan_rate: f64,     // mV/s
    sample_rate: f64,   // Hz
    voltage_start: f64, // V
    voltage_step1: f64, // V
    voltage_step2: f64, // V
}

fn parse_parameters(lines: &[&str]) -> ExperimentParams {
    let mut params = ExperimentParams {
        scan_rate: 0.0,
        sample_rate: 0.0,
        voltage_start: 0.0,
        voltage_step1: 0.0,
        voltage_step2: 0.0,
    };

    for line in lines.iter() {
        if line.starts_with("ExpParmas:") {
            // 查找 ExpInfo 部分在行中的位置
            if let Some(exp_info_start) = line.find("ExpInfo:") {
                // 只处理 ExpInfo 部分（跳过前面的 "ExpInfo:"）
                let exp_info = &line[exp_info_start + 8..];

                // 使用正则表达式按逗号分割字符串（考虑可能的空格）
                let re = Regex::new(r",\s*").unwrap(); // 匹配逗号及后续空格
                let parts: Vec<&str> = re.split(exp_info).collect();

                for part in parts {
                    let part = part.trim(); // 去除两端空白
                    if part.starts_with("Scan Rate(mV/s):") {
                        // 提取数字部分（忽略单位描述）
                        if let Some(value) = part
                            .split(':')
                            .nth(1)
                            .and_then(|s| s.split_whitespace().next())
                        {
                            params.scan_rate = value.trim().parse().unwrap_or_default();
                        }
                    } else if part.starts_with("Freq(Hz):") {
                        if let Some(value) = part.split(':').nth(1) {
                            params.sample_rate = value.trim().parse().unwrap_or_default();
                        }
                    } else if part.starts_with("Init E(V):") {
                        if let Some(value) = part
                            .split(':')
                            .nth(1)
                            .and_then(|s| s.split_whitespace().next())
                        {
                            params.voltage_start = value.trim().parse().unwrap_or_default();
                        }
                    } else if part.starts_with("Step1 E(V):") {
                        if let Some(value) = part
                            .split(':')
                            .nth(1)
                            .and_then(|s| s.split_whitespace().next())
                        {
                            params.voltage_step1 = value.trim().parse().unwrap_or_default();
                        }
                    } else if part.starts_with("Step2 E(V):") {
                        if let Some(value) = part
                            .split(':')
                            .nth(1)
                            .and_then(|s| s.split_whitespace().next())
                        {
                            params.voltage_step2 = value.trim().parse().unwrap_or_default();
                        }
                    }
                }
            }
        }
    }

    params
}

fn calculate_points_to_keep(params: &ExperimentParams) -> usize {
    // 计算单程扫描的电压跨度
    let total_span = (params.voltage_step2 - params.voltage_step1).abs() * 2.0;

    // 计算一个周期所需时间（秒）
    let scan_rate_v_per_s = params.scan_rate / 1000.0; // 转换为V/s
    let cycle_time = total_span / scan_rate_v_per_s;

    // 计算一个周期的数据点数
    (cycle_time * params.sample_rate + 1.0) as usize
}
