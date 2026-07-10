use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum EisFormat {
    CorrTestZ60,
    CorrTestText,
    Csv,
    Tsv,
    HeaderlessThreeColumn,
    LegacyCleanCsv,
}
