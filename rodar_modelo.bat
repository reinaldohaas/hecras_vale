@echo off
setlocal enabledelayedexpansion
rem ===========================================================================
rem  Roda o modelo hidrodinamico da Bacia do Itajai (HEC-RAS 7.0.1).
rem
rem    rodar_modelo.bat                 cheia sintetica
rem    rodar_modelo.bat 1983            evento historico com chuva real
rem    rodar_modelo.bat --todos         1983, 2008, 2011 e 2023
rem    rodar_modelo.bat 1983 --sem-barragens
rem
rem  Nao depende de nenhum assistente: e so Python + HEC-RAS.
rem ===========================================================================
cd /d "%~dp0"

rem --- 1. Localiza o Python -------------------------------------------------
set "PY="
for %%P in (
  "%USERPROFILE%\miniforge3\python.exe"
  "%USERPROFILE%\miniconda3\python.exe"
  "%USERPROFILE%\anaconda3\python.exe"
  "%LOCALAPPDATA%\miniforge3\python.exe"
  "C:\ProgramData\miniforge3\python.exe"
) do (
  if not defined PY if exist %%P set "PY=%%~P"
)
if not defined PY (
  for /f "delims=" %%P in ('where python 2^>nul') do (
    if not defined PY set "PY=%%P"
  )
)
if not defined PY (
  echo [ERRO] Python nao encontrado. Instale o Miniforge ou ajuste este .bat.
  exit /b 1
)
echo   Python: %PY%

rem --- 2. Blinda o PATH ------------------------------------------------------
rem O MSYS2/Git Bash coloca C:\msys64\...\bin no PATH, e as DLLs de libpng e
rem freetype de la se sobrepoem as do conda: o matplotlib quebra na hora de
rem salvar figura e o rasterio derruba o processo. Colocando o Python e suas
rem pastas de biblioteca na FRENTE, as DLLs corretas vencem.
for %%I in ("%PY%") do set "PYDIR=%%~dpI"
set "PATH=%PYDIR%;%PYDIR%Library\bin;%PYDIR%Library\mingw-w64\bin;%PYDIR%Scripts;%PATH%"

rem --- 3. Confere o HEC-RAS --------------------------------------------------
if not exist "C:\Program Files (x86)\HEC\HEC-RAS\7.0.1\Ras.exe" (
  echo   [AVISO] HEC-RAS 7.0.1 nao encontrado no caminho padrao.
  echo           A simulacao vai falhar; a geracao dos arquivos ainda funciona.
)

rem --- 4. Dependencias ------------------------------------------------------
"%PY%" -c "import numpy, geopandas, rasterio, shapely, pyproj, h5py, scipy" 2>nul
if errorlevel 1 (
  echo   Instalando dependencias que faltam...
  "%PY%" -m pip install --quiet numpy geopandas rasterio shapely pyproj h5py scipy
)
"%PY%" -c "import win32com.client" 2>nul
if errorlevel 1 (
  echo   Instalando pywin32 ^(interface COM do HEC-RAS^)...
  "%PY%" -m pip install --quiet pywin32
)

rem --- 5. Roda --------------------------------------------------------------
"%PY%" rodar_modelo.py %*
set "RC=%ERRORLEVEL%"

if not "%RC%"=="0" (
  echo.
  echo   [ERRO] O pipeline terminou com erro ^(codigo %RC%^).
)
rem A dica vale mesmo quando o HEC-RAS nao fecha: o motor proprio roda sempre e
rem as paginas sao geradas de qualquer jeito. Antes ela so aparecia em caso de
rem sucesso, e mandava servir app/ enquanto as paginas saiam na raiz.
echo.
echo   Para ver no navegador:
echo     "%PY%" -m http.server 8050 --directory app
echo.
echo   e abra:
echo     http://localhost:8050/                              interface completa
echo     http://localhost:8050/^<PROJETO^>_cheia_motor.html    perfil animado ^(motor^)
echo     http://localhost:8050/^<PROJETO^>_cheia_hecras.html   idem, pelo HEC-RAS
exit /b %RC%
