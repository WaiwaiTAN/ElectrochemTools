use anyhow::{Context, Result, bail};
use chrono::{DateTime, FixedOffset};
use clap::Parser;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

/// 提取文件的元数据和测试时间
fn extract_metadata(file_path: &PathBuf) -> Result<(Vec<String>, Option<DateTime<FixedOffset>>)> {
    let file =
        File::open(file_path).with_context(|| format!("文件未找到: {}", file_path.display()))?;
    let reader = BufReader::new(file);
    let mut metadata = Vec::new();
    let mut test_datetime = None;
    let mut data_start = false;
    let mut end_comments_found = false;

    for line in reader.lines() {
        let line = line.with_context(|| format!("读取文件失败: {}", file_path.display()))?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if end_comments_found {
            break;
        }

        if trimmed.contains("End Comments") {
            end_comments_found = true;
            metadata.push(line);
            continue;
        }

        metadata.push(line.clone());

        if data_start {
            continue;
        }

        if trimmed.starts_with("Data:") && trimmed.contains("Time:") {
            let parts: Vec<&str> = trimmed.split("Time:").collect();
            if parts.len() != 2 {
                println!("警告: {} 中的时间格式无效", file_path.display());
                continue;
            }

            let cleaned = parts[0].replace("Data:", "");
            let date_str = cleaned.trim(); 

            let time_str = parts[1].trim(); 

            match parse_datetime(&format!("{} {}", date_str, time_str)) {
                Ok(dt) => {
                    test_datetime = Some(dt);
                    data_start = true;
                }
                Err(_) => {
                    println!(
                        "警告: {} 中的时间格式无效: {}",
                        file_path.display(),
                        trimmed
                    );
                }
            }
        }
    }

    Ok((metadata, test_datetime))
}

/// 解析日期时间字符串
fn parse_datetime(datetime_str: &str) -> Result<DateTime<FixedOffset>> {
    DateTime::parse_from_str(&format!("{} +08:00", datetime_str), "%Y-%m-%d %H:%M:%S %z")
        .with_context(|| format!("无效的时间格式: {}", datetime_str))
}

/// 提取文件的数据部分
fn extract_data(file_path: &PathBuf) -> Result<Vec<String>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut data = Vec::new();
    let mut data_start = false;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if trimmed.contains("End Comments") {
            data_start = true;
            continue;
        }

        if data_start && !trimmed.starts_with('#') {
            data.push(line);
        }
    }

    Ok(data)
}

/// 调整时间列偏移
fn adjust_time(data_line: &str, time_offset: i64) -> Option<String> {
    let parts: Vec<&str> = data_line.split('\t').collect();
    if parts.len() != 3 {
        println!("警告: 无效的数据行格式 - {}", data_line);
        return None;
    }

    match parts[2].parse::<f64>() {
        Ok(t) => {
            let adjusted_time = t + (time_offset as f64) / 1000.0;
            Some(format!("{}\t{}\t{:.5E}", parts[0], parts[1], adjusted_time))
        }
        Err(_) => {
            println!("警告: 无效的时间值 - {}", parts[2]);
            None
        }
    }
}

/// 命令行参数结构
#[derive(Parser, Debug)]
#[clap(author, version, about = "合并COR文件并修正时间")]
struct Args {
    /// 输入文件列表
    #[clap(short = 'i', long = "input", required = true, value_delimiter = ' ', num_args = 1..)]
    input_files: Vec<PathBuf>,

    /// 输出文件
    #[clap(short = 'o', long = "output", required = true)]
    output_file: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // 验证输入文件是否存在
    for file in &args.input_files {
        if !file.exists() {
            bail!("错误: 输入文件不存在 - {}", file.display());
        }
    }

    // 收集文件信息 (路径, 测试时间)
    let mut files_info = Vec::new();
    for file in &args.input_files {
        match extract_metadata(file) {
            Ok((_, Some(dt))) => {
                files_info.push((file.clone(), dt));
            }
            Ok((_, None)) => {
                bail!("错误: {} 中未找到有效时间戳", file.display());
            }
            Err(e) => {
                bail!("错误: 处理 {} 失败 - {}", file.display(), e);
            }
        }
    }

    // 按测试时间排序
    if files_info.is_empty() {
        bail!("错误: 没有有效的输入文件");
    }

    files_info.sort_by(|a, b| a.1.cmp(&b.1));
    let (base_file, base_time) = &files_info[0];
    println!("基础文件: {}, 测试时间: {}", base_file.display(), base_time);

    // 创建输出文件
    let out_file = File::create(&args.output_file)
        .with_context(|| format!("无法创建输出文件: {}", args.output_file.display()))?;
    let mut writer = BufWriter::new(out_file);

    // 写入基础元数据
    let (base_metadata, _) = extract_metadata(base_file)?;
    for line in base_metadata {
        writeln!(writer, "{}", line)?;
    }

    // 处理并写入数据
    for (file, test_time) in &files_info {
        let time_offset = (*test_time - *base_time)
            .to_std()
            .map(|d| d.as_millis() as i64)
            .unwrap_or(0);

        println!(
            "处理文件: {}, 时间偏移: {:.2} 秒",
            file.display(),
            time_offset as f64 / 1000.0
        );

        let data_lines = extract_data(file)?;
        for line in data_lines {
            if let Some(adjusted) = adjust_time(&line, time_offset) {
                writeln!(writer, "{}", adjusted)?;
            }
        }
    }

    println!("合并完成! 结果已保存到: {}", args.output_file.display());

    Ok(())
}
