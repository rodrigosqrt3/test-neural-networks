# Baixa a EUPPBench v1.0, dados de ESTAÇÕES (Demaeyer et al. 2023, ESSD; Zenodo 7708362):
# 18,3 GB zipados, Zarr; ~35 GB descompactados. Critério 2 da fase 1 (PLANO.md).
# Retoma de onde parou se cair (curl -C -), confere o md5 e descompacta.
#
# Uso (na pasta do projeto; pode rodar junto com o exp6, só usa rede e disco):
#   powershell -ExecutionPolicy Bypass -File analises\19_baixar_euppbench.ps1

$raiz = Split-Path $PSScriptRoot -Parent
$pasta = Join-Path $raiz "dados\brutos\euppbench"
New-Item -ItemType Directory -Force $pasta | Out-Null
$zip = Join-Path $pasta "EUPPBench-stations.zip"
$url = "https://zenodo.org/records/7708362/files/EUPPBench-stations.zip?download=1"
$md5_certo = "e409457279b3494d18f2dfb41f3f449b"

Write-Host "$(Get-Date -Format HH:mm) baixando para $zip"
$tentativa = 0
do {
  $tentativa++
  & curl.exe -L -C - --retry 20 --retry-delay 30 --retry-all-errors -o $zip $url
  $ok = $LASTEXITCODE -eq 0
  if (-not $ok) { Write-Host "curl saiu com codigo $LASTEXITCODE; tentando de novo ($tentativa)"; Start-Sleep -Seconds 60 }
} while (-not $ok -and $tentativa -lt 10)
if (-not $ok) { Write-Host "ERRO: download nao terminou. Rode de novo que ele continua de onde parou."; exit 1 }

Write-Host "$(Get-Date -Format HH:mm) conferindo md5 (demora alguns minutos)"
$md5 = (Get-FileHash $zip -Algorithm MD5).Hash.ToLower()
if ($md5 -ne $md5_certo) { Write-Host "ERRO: md5 $md5 diferente de $md5_certo. Apague o zip e rode de novo."; exit 1 }

Write-Host "$(Get-Date -Format HH:mm) descompactando"
& tar.exe -xf $zip -C $pasta
if ($LASTEXITCODE -ne 0) { Write-Host "ERRO ao descompactar"; exit 1 }
Get-ChildItem $pasta | Select-Object Name, @{n='GB';e={[math]::Round(((Get-ChildItem $_.FullName -Recurse -File -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum)/1GB,2)}}
Write-Host "$(Get-Date -Format HH:mm) pronto. O zip pode ser apagado depois de conferir os dados."
