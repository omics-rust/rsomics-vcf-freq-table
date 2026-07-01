use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, ToolMeta};
use rsomics_vcf_freq_table::{Mode, freq_table};

use crate::table;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-vcf-freq-table",
    version,
    about = "Allele frequency and count tables from a VCF (vcftools --freq/--freq2/--counts/--counts2)"
)]
pub struct Cli {
    /// Input VCF or VCF.gz file; omit or pass `-` to read stdin.
    #[arg(value_name = "VCF")]
    pub input: Option<PathBuf>,

    /// Output mode: freq | freq2 | counts | counts2.
    #[arg(long, default_value = "freq", value_name = "MODE")]
    pub mode: Mode,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let path = self.input.as_deref();
        let mode = self.mode;
        let rows = freq_table(path, mode)?;
        let tsv = table::to_tsv(&rows, mode);

        if self.common.json {
            let env = serde_json::json!({
                "schema_version": rsomics_common::SCHEMA_VERSION,
                "tool": META.name,
                "tool_version": META.version,
                "status": "ok",
                "result": {
                    "mode": mode,
                    "rows": rows,
                },
            });
            println!("{}", serde_json::to_string(&env).unwrap_or_default());
        } else {
            print!("{tsv}");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
