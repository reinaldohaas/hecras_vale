@echo off
REM ===================================================================
REM  Modelo do Vale do Itajai -- execucao em serie, autocorretiva.
REM
REM  Uso:
REM     rodar_vale.bat                 Copernicus, tudo, com autocorrecao
REM     rodar_vale.bat 5-10            SO os passos 5 a 10 (reaproveita 1-4)
REM     rodar_vale.bat 8               SO o passo 8
REM     rodar_vale.bat passos          lista os passos e o que ja foi feito
REM     rodar_vale.bat 1-10 selecao=Iraputa projeto=so_iraputa
REM                                    UM RIO SO -- roda em minutos, e o que
REM                                    quebrar quebra sozinho. Qualquer opcao
REM                                    (chave=valor) vai adiante; veja-as em
REM                                    "python -m vale opcoes".
REM     rodar_vale.bat ambiente        so confere o ambiente e o que falta
REM     rodar_vale.bat qaqc            qualifica o SIG-SC contra o Copernicus
REM     rodar_vale.bat sigsc [passos]  roda com o MDT do SIG-SC a 10 m
REM     rodar_vale.bat 1983  [passos]  evento de 1983 (Sul e Norte operando)
REM
REM  REFAZER SO O QUE MUDOU. Cada passo grava o resultado em modelo\estado.pkl
REM  e o seguinte le de la, entao dar um intervalo reaproveita tudo que veio
REM  antes. Mexeu no perfil (passo 5) ou na calha (passo 6)?  "5-10" pula
REM  catalogo, eixos, terreno e o corte das secoes. Mexeu so na visualizacao?
REM  "10" sozinho. Rodar "tudo" a cada ajuste era o modo errado de usar isto,
REM  e a culpa era deste arquivo, que so sabia dizer "tudo".
REM
REM  O interpretador NAO fica cravado aqui. Ele e procurado, e o que falta
REM  nele e listado com o comando de instalacao.
REM
REM  Tudo que a execucao produz vai para modelo\ .
REM ===================================================================
setlocal enabledelayedexpansion

cd /d "%~dp0"
if not exist modelo mkdir modelo

REM ---- procura um interpretador que consiga ao menos importar o pacote ----
set PY=
for %%P in (
    "C:\Users\haas\miniforge3\envs\hecras-qc\python.exe"
    "C:\Users\haas\miniforge3\envs\vale\python.exe"
    "C:\Users\haas\miniforge3\python.exe"
) do (
    if exist %%~P if not defined PY set PY=%%~P
)
if not defined PY (
    for /f "delims=" %%A in ('where python 2^>nul') do (
        if not defined PY set PY=%%A
    )
)
if not defined PY (
    echo.
    echo Nao encontrei nenhum Python. Instale o Miniforge ou ponha o
    echo python.exe no PATH e rode este arquivo de novo.
    exit /b 2
)
echo Interpretador: %PY%

REM ---- confere as dependencias ANTES de tentar rodar ----
"%PY%" -m vale.ambiente
if errorlevel 1 (
    echo.
    echo Instale o que falta acima e rode de novo.
    echo Se preferir criar um ambiente novo e limpo:
    echo    mamba create -n vale -c conda-forge -y python=3.12 numpy shapely geopandas pyproj rasterio pandas scipy h5py
    echo    mamba activate vale
    echo    pip install ras-commander
    exit /b 1
)

if "%1"=="ambiente" goto FIM
if "%1"=="passos" goto LISTA
if "%1"=="qaqc" goto QAQC

REM ---- separa o modo do intervalo de passos ----
REM  Sem argumento nenhum, ou com um intervalo no primeiro argumento, o modo e
REM  Copernicus. Com "sigsc" ou "1983", o intervalo vem no segundo.
REM  set com aspas: "set PASSOS=%2 )" deixa o espaco de antes do parentese
REM  dentro da variavel, e "sigsc" sozinho virava um passo chamado " ".
set "MODO=copernicus"
set "PASSOS=%~1"
if "%~1"=="sigsc" ( set "MODO=sigsc" & set "PASSOS=%~2" )
if "%~1"=="1983"  ( set "MODO=1983"  & set "PASSOS=%~2" )
if "%PASSOS%"=="" set "PASSOS=tudo"

REM  Opcoes extras seguem adiante: "rodar_vale.bat 1-10 selecao=Iraputa
REM  projeto=so_iraputa" roda um rio so, em projeto proprio, sem tocar no
REM  modelo cheio. Rio isolado quebra sozinho e em minutos.
REM  Reconstroi a linha em vez de iterar %1: o batch trata "=" como separador
REM  de argumentos, entao "selecao=Iraputa" chegaria partido em dois e a opcao
REM  se perderia calada. Consome o que ja foi lido (o modo, quando ha, e o
REM  intervalo) e passa o resto adiante.
set "TODOS=%*"
set "RESTO="
set "EXTRA="
if defined TODOS for /f "tokens=1,* delims= " %%a in ("!TODOS!") do set "RESTO=%%b"
if "%MODO%"=="copernicus" (
    set "EXTRA=!RESTO!"
) else (
    if defined RESTO for /f "tokens=1,* delims= " %%a in ("!RESTO!") do set "EXTRA=%%b"
)
if defined EXTRA if not "!EXTRA!"=="" echo  opcoes: !EXTRA!

echo.
echo  passos: %PASSOS%
if not "%PASSOS%"=="tudo" (
    echo  os passos anteriores sao lidos de modelo\estado.pkl, nao refeitos
)

if "%MODO%"=="sigsc" goto SIGSC
if "%MODO%"=="1983"  goto EVENTO
goto COPERNICUS

:COPERNICUS
echo.
echo ============================================================
echo  COPERNICUS -- 30 m, cobre tudo, sem vazio
echo  Modelo de SUPERFICIE: contem copa de mata e a lamina d'agua.
echo  Por isso a escavacao da calha fica DESLIGADA (o programa avisa).
echo ============================================================
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '%PY%' -u -m vale %PASSOS% --auto fonte=copernicus !EXTRA! 2>&1 | Tee-Object -FilePath modelo\execucao.log"
goto FIM

:SIGSC
echo.
echo ============================================================
echo  SIG-SC -- MDT a 10 m, solo exposto
echo  Qualifique ANTES:  rodar_vale.bat qaqc
echo  O passo do terreno le 55 GB e leva cerca de duas horas.
echo ============================================================
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '%PY%' -u -m vale %PASSOS% --auto fonte=sigsc res_sigsc=10 terreno_hdf=false !EXTRA! 2>&1 | Tee-Object -FilePath modelo\execucao.log"
goto FIM

:EVENTO
echo.
echo ============================================================
echo  EVENTO DE 1983 -- so Sul e Norte operavam (Oeste em construcao)
echo ============================================================
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '%PY%' -u -m vale %PASSOS% --auto fonte=copernicus evento=1983 projeto=vale_1983 !EXTRA! 2>&1 | Tee-Object -FilePath modelo\execucao.log"
goto FIM

:LISTA
echo.
"%PY%" -m vale
goto FIM

:QAQC
echo.
echo ============================================================
echo  QA/QC do SIG-SC contra o Copernicus
echo ============================================================
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '%PY%' -u -m vale.qaqc 2>&1 | Tee-Object -FilePath modelo\qaqc.log"
goto FIM

:FIM
echo.
echo ------------------------------------------------------------
echo  Log em modelo\execucao.log
echo  Passos e estado:   rodar_vale.bat passos
echo  So um intervalo:   rodar_vale.bat 5-10
echo  Opcoes:            "%PY%" -m vale opcoes
echo ------------------------------------------------------------
endlocal
