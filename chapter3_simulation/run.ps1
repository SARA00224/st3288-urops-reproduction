# Run using: powershell -File .\papers\01_plasma_proteome_MR\pqtl_beta_reproduction\run.ps1
$ErrorActionPreference = 'Stop'
& Rscript --vanilla (Join-Path $PSScriptRoot 'run_all.R')
if ($LASTEXITCODE -ne 0) { throw "R analysis failed: $LASTEXITCODE" }
