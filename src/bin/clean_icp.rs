use anyhow::{Context, Result, bail};
use calamine::{DataType, Range, Reader, Xlsx, open_workbook};
use clap::Parser;
use csv::Writer;
use regex::Regex;
use std::path::PathBuf;

/// 命令行参数结构
#[derive(Parser, Debug)]
#[clap(author, version, about = "处理Excel数据并过滤特定格式的样本")]
struct Args {
    /// 输入Excel文件路径 (.xls 或 .xlsx)
    #[clap(short = 'i', long)]
    input: PathBuf,

    /// 输出CSV文件路径 (.csv)
    #[clap(short = 'o', long)]
    output: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // 验证输入文件是否存在
    if !args.input.exists() {
        bail!("错误: 输入文件不存在 - {}", args.input.display());
    }

    // 处理Excel数据
    let (headers, filtered_rows) = process_excel(&args.input)?;

    // 写入结果到CSV文件
    write_output_csv(&args.output, headers, filtered_rows)?;

    println!("处理成功! 结果已保存至: {}", args.output.display());
    Ok(())
}

/// 处理Excel文件的核心逻辑
fn process_excel(input_path: &PathBuf) -> Result<(Vec<String>, Vec<Vec<calamine::Data>>)> {
    // 打开工作簿
    let mut workbook: Xlsx<_> = open_workbook(input_path)
        .with_context(|| format!("无法打开文件: {}", input_path.display()))?;
    
    
    let range = workbook
        .worksheet_range("Sheet1")
        .map_err(|e| anyhow::anyhow!("读取 Sheet1 失败: {}", e))?;

    // 提取表头 (第一行)
    let headers = get_headers(&range)?;

    // 编译正则表达式匹配 T[1-9]-[0-9] 格式
    let re = Regex::new(r"^T[1-9]-[0-9]$")?;

    // 过滤有效数据行
    let mut filtered_rows = Vec::new();
    for row in range.rows().skip(1) {
        // 跳过空行
        if row.iter().all(|cell| cell.is_empty()) {
            continue;
        }

        // 检查Label列是否符合格式 (第4列，索引3)
        if let Some(label_cell) = row.get(3) {
            if let calamine::Data::String(label) = label_cell {
                if re.is_match(&label) {
                    filtered_rows.push(row.to_vec());
                }
            }
        }
    }

    Ok((headers, filtered_rows))
}

/// 提取表头行
fn get_headers(range: &Range<calamine::Data>) -> Result<Vec<String>> {
    if range.is_empty() {
        bail!("工作表为空");
    }

    range
        .rows()
        .next()
        .unwrap()
        .iter()
        .map(|cell| {
            Ok(match cell {
                calamine::Data::String(s) => s.clone(),
                _ => cell.to_string(),
            })
        })
        .collect()
}

/// 写入结果到CSV文件
fn write_output_csv(
    output_path: &PathBuf,
    headers: Vec<String>,
    rows: Vec<Vec<calamine::Data>>,
) -> Result<()> {
    // 创建CSV写入器
    let mut wtr = Writer::from_path(output_path)?;

    // 写入表头
    wtr.write_record(&headers)?;

    // 写入数据行
    for row in rows {
        let record: Vec<String> = row.iter().map(|cell| cell.to_string()).collect();
        wtr.write_record(&record)?;
    }

    wtr.flush()?;
    Ok(())
}
