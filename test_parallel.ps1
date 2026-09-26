# ============================================================================
# Run masterScript.R in batches of 2
# ============================================================================

$projectDir = "C:\GitHub\breedingBirdIndicator_OccupancyModelling"

Set-Location $projectDir

# ---------------------------------------------------------------------------
# Define jobs
# ---------------------------------------------------------------------------

$jobs = @(
    @{ version = 1; option = "TRUE"; name = "v1_TRUE" }
    @{ version = 1; option = "FALSE"; name = "v1_FALSE" }
)

# ---------------------------------------------------------------------------
# Remove any existing PowerShell jobs
# ---------------------------------------------------------------------------

Get-Job | Remove-Job -Force -ErrorAction SilentlyContinue

# ---------------------------------------------------------------------------
# Process jobs in batches of 2
# ---------------------------------------------------------------------------

$batchSize = 2

for ($i = 0; $i -lt $jobs.Count; $i += $batchSize) {

    $end = [Math]::Min($i + $batchSize - 1, $jobs.Count - 1)

    $batch = $jobs[$i..$end]

    Write-Host ""
    Write-Host "============================================" -ForegroundColor Cyan
    Write-Host "Starting batch $([int]($i / $batchSize) + 1)" -ForegroundColor Cyan
    Write-Host "============================================" -ForegroundColor Cyan

    # -----------------------------------------------------------------------
    # Start jobs in this batch
    # -----------------------------------------------------------------------

    foreach ($job in $batch) {

        Write-Host "Starting $($job.name) ..." -ForegroundColor Yellow

        Start-Job -Name $job.name -ScriptBlock {

            param(
                $projectDir,
                $version,
                $option
            )

            Set-Location $projectDir

            & Rscript --vanilla `
                "scripts/02_masterScriptRstudio.R" `
                $version `
                $option

            # Return R exit code
            $LASTEXITCODE

        } -ArgumentList `
            $projectDir,
            $job.version,
            $job.option
    }

    # -----------------------------------------------------------------------
    # Wait for ALL jobs in this batch
    # -----------------------------------------------------------------------

    Write-Host ""
    Write-Host "Waiting for batch to finish..." -ForegroundColor Yellow

    $batchNames = $batch.name

    Get-Job |
        Where-Object { $batchNames -contains $_.Name } |
        Wait-Job | Out-Null

    # -----------------------------------------------------------------------
    # Check results
    # -----------------------------------------------------------------------

    Write-Host ""
    Write-Host "Batch completed:" -ForegroundColor Green

    foreach ($job in $batch) {

        $psJob = Get-Job -Name $job.name

        $result = Receive-Job -Name $job.name -Keep

        Write-Host "$($job.name): " -NoNewline

        if ($psJob.State -eq "Completed" -and $result -contains 0) {
            Write-Host "SUCCESS" -ForegroundColor Green
        }
        else {
            Write-Host "FAILED" -ForegroundColor Red
            Write-Host "State: $($psJob.State)" -ForegroundColor Red
            Write-Host "Output:" -ForegroundColor Red
            $result
        }
    }

    # -----------------------------------------------------------------------
    # Remove completed jobs before starting next batch
    # -----------------------------------------------------------------------

    Get-Job |
        Where-Object { $batchNames -contains $_.Name } |
        Remove-Job -Force

    Write-Host ""
    Write-Host "Batch $([int]($i / $batchSize) + 1) finished." -ForegroundColor Green
}

# ---------------------------------------------------------------------------
# Finished
# ---------------------------------------------------------------------------

Write-Host ""
Write-Host "============================================" -ForegroundColor Green
Write-Host "ALL FOUR R JOBS FINISHED" -ForegroundColor Green
Write-Host "============================================" -ForegroundColor Green