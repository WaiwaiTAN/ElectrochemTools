use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::fs;
use std::io::Read;
use std::path::{Component, Path, PathBuf};

pub const SCHEMA_VERSION: u32 = 1;
const RESUME_HINT: &str = "existing output is not resumable; inspect/remove it, use --overwrite, or choose another --out-root";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    Running,
    Success,
    Failed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ManifestInput {
    pub display_path: String,
    pub size_bytes: u64,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    pub schema_version: u32,
    pub program_version: String,
    pub command: String,
    pub status: RunStatus,
    pub input: ManifestInput,
    pub configuration: Value,
    pub configuration_sha256: String,
    pub outputs: Vec<PathBuf>,
    pub error: Option<String>,
}

impl RunManifest {
    pub fn new(command: &str, input: &Path, configuration: Value) -> Result<Self> {
        Ok(Self {
            schema_version: SCHEMA_VERSION,
            program_version: env!("CARGO_PKG_VERSION").to_string(),
            command: command.to_string(),
            status: RunStatus::Running,
            input: inspect_input(input)?,
            configuration_sha256: value_sha256(&configuration)?,
            configuration,
            outputs: Vec::new(),
            error: None,
        })
    }
}

pub fn inspect_input(path: &Path) -> Result<ManifestInput> {
    let metadata =
        fs::metadata(path).with_context(|| format!("failed to inspect {}", path.display()))?;
    Ok(ManifestInput {
        display_path: path.display().to_string(),
        size_bytes: metadata.len(),
        sha256: file_sha256(path)?,
    })
}

pub fn file_sha256(path: &Path) -> Result<String> {
    let mut file = fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

pub fn value_sha256(value: &Value) -> Result<String> {
    let bytes = serde_json::to_vec(value)?;
    Ok(format!("{:x}", Sha256::digest(bytes)))
}

pub fn write_atomic(output_dir: &Path, manifest: &RunManifest) -> Result<()> {
    let path = output_dir.join("run.json");
    let temporary = output_dir.join("run.json.tmp");
    fs::write(&temporary, serde_json::to_vec_pretty(manifest)?)?;
    if path.exists() {
        fs::remove_file(&path)?;
    }
    fs::rename(&temporary, &path)?;
    Ok(())
}

pub fn verify_resume(output_dir: &Path, expected: &RunManifest) -> Result<bool> {
    let reject = |reason: &str| anyhow::anyhow!("{reason}; {RESUME_HINT}");
    let bytes = fs::read(output_dir.join("run.json"))
        .map_err(|e| reject(&format!("missing run.json: {e}")))?;
    let existing: RunManifest =
        serde_json::from_slice(&bytes).map_err(|e| reject(&format!("invalid run.json: {e}")))?;
    if existing.schema_version != SCHEMA_VERSION {
        return Err(reject("unsupported run.json schema"));
    }
    if existing.status != RunStatus::Success {
        return Err(reject("run status is not success"));
    }
    if existing.command != expected.command {
        return Err(reject("command differs"));
    }
    if existing.input.sha256 != expected.input.sha256 {
        return Err(reject("input SHA-256 differs"));
    }
    if existing.configuration_sha256 != expected.configuration_sha256 {
        return Err(reject("configuration differs"));
    }
    if existing.outputs.is_empty() {
        return Err(reject("manifest declares no outputs"));
    }
    let canonical_output_dir = fs::canonicalize(output_dir)
        .map_err(|e| reject(&format!("cannot resolve sample directory: {e}")))?;
    for relative in &existing.outputs {
        if relative.is_absolute()
            || relative.components().any(|part| {
                matches!(
                    part,
                    Component::ParentDir | Component::RootDir | Component::Prefix(_)
                )
            })
        {
            return Err(reject("manifest output path escapes the sample directory"));
        }
        let declared = output_dir.join(relative);
        if !declared.is_file() {
            return Err(reject(&format!(
                "declared output is missing: {}",
                relative.display()
            )));
        }
        let canonical_declared = fs::canonicalize(&declared)
            .map_err(|e| reject(&format!("cannot resolve declared output: {e}")))?;
        if !canonical_declared.starts_with(&canonical_output_dir) {
            return Err(reject(
                "manifest output resolves outside the sample directory",
            ));
        }
    }
    Ok(false)
}

pub fn execute_manifested<F>(
    output_dir: &Path,
    mut manifest: RunManifest,
    operation: F,
) -> Result<()>
where
    F: FnOnce() -> Result<Vec<PathBuf>>,
{
    write_atomic(output_dir, &manifest)?;
    match operation() {
        Ok(outputs) => {
            manifest.status = RunStatus::Success;
            manifest.outputs = outputs;
            write_atomic(output_dir, &manifest)?;
            Ok(())
        }
        Err(error) => {
            manifest.status = RunStatus::Failed;
            manifest.error = Some(format!("{error:#}"));
            let _ = write_atomic(output_dir, &manifest);
            Err(error)
        }
    }
}

pub fn ensure_relative_outputs(output_dir: &Path, names: &[&str]) -> Result<Vec<PathBuf>> {
    let outputs = names.iter().map(PathBuf::from).collect::<Vec<_>>();
    for output in &outputs {
        if !output_dir.join(output).is_file() {
            bail!("expected output was not created: {}", output.display());
        }
    }
    Ok(outputs)
}

pub fn collect_output_files(output_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut outputs = fs::read_dir(output_dir)?
        .filter_map(|entry| entry.ok())
        .filter(|entry| {
            entry
                .file_type()
                .map(|kind| kind.is_file())
                .unwrap_or(false)
        })
        .map(|entry| PathBuf::from(entry.file_name()))
        .filter(|path| path != Path::new("run.json") && path != Path::new("run.json.tmp"))
        .collect::<Vec<_>>();
    outputs.sort();
    if outputs.is_empty() {
        bail!("calculation produced no output files");
    }
    Ok(outputs)
}
