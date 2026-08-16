#!/usr/bin/env bash
# ===========================================================================
#  Roda o modelo hidrodinamico da Bacia do Itajai (HEC-RAS 7.0.1).
#
#    ./rodar_modelo.sh                cheia sintetica
#    ./rodar_modelo.sh 1983           evento historico com chuva real
#    ./rodar_modelo.sh --todos        1983, 2008, 2011 e 2023
#    ./rodar_modelo.sh 1983 --sem-barragens
#
#  Nao depende de nenhum assistente: e so Python + HEC-RAS.
#
#  Nota: a simulacao usa a interface COM do HEC-RAS, que so existe no
#  Windows. Rodando por Git Bash / WSL com o HEC-RAS instalado, funciona.
#  Em Linux puro, as etapas de geracao e de mancha rodam; a simulacao nao.
# ===========================================================================
set -u
cd "$(dirname "$0")"

# --- 1. Blinda o PATH (ANTES de escolher o Python) ---------------------------
# O MSYS2/Git Bash coloca /c/msys64/*/bin no PATH, e as DLLs de libpng e
# freetype de la se sobrepoem as do conda: o matplotlib quebra ao salvar figura
# e o rasterio derruba o processo.
#
# A limpeza tem de vir ANTES da deteccao, nao depois. Com o msys64 na frente, o
# proprio teste de import do miniforge falha (rasterio nao carrega as DLLs
# certas), o script conclui "esse Python nao serve" e cai no /ucrt64/bin/python3
# -- que nao tem numpy nenhum. Limpando primeiro, o teste mede o Python, nao o
# PATH; e de quebra o python do msys some da lista de candidatos.
PATH="$(echo "$PATH" | tr ':' '\n' | grep -viE '^/c/msys64|^/mingw|^/ucrt64' | paste -sd: -)"
export PATH

# --- 2. Localiza o Python ----------------------------------------------------
# Nao basta o interpretador EXISTIR: cada candidato e testado importando o que
# o modelo usa, e o primeiro que passar vence.
#
# A raiz do perfil precisa de mais de uma via: num shell MSYS2, $HOME e
# /home/haas (nao o perfil do Windows) e $USERPROFILE nem existe -- entao
# procurar so por $HOME nao acha o miniforge e o script cai no python do
# Store, que tambem nao tem numpy. Tenta todas e usa a primeira que servir.
RAIZES="${HOME:-}"
[ -n "${USERPROFILE:-}" ] && \
  RAIZES="$RAIZES $(cygpath -u "$USERPROFILE" 2>/dev/null || echo "${USERPROFILE:-}")"
QUEM="$(whoami 2>/dev/null | tr -d '\\r' | sed 's|.*\\\\||')"
[ -n "$QUEM" ] && RAIZES="$RAIZES /c/Users/$QUEM $(cygpath -H 2>/dev/null)/$QUEM"

CANDIDATOS=""
for r in $RAIZES; do
  [ -n "$r" ] || continue
  for dist in miniforge3 miniconda3 anaconda3; do
    CANDIDATOS="$CANDIDATOS $r/$dist/python.exe $r/$dist/bin/python"
  done
done
CANDIDATOS="$CANDIDATOS $(command -v python3 || true) $(command -v python || true)"

PY=""
PY_FALLBACK=""
for c in $CANDIDATOS; do
  [ -n "$c" ] && [ -x "$c" ] || continue
  if "$c" -c "import numpy, geopandas, rasterio, shapely, pyproj, h5py" 2>/dev/null; then
    PY="$c"; break
  fi
  [ -z "$PY_FALLBACK" ] && PY_FALLBACK="$c"
done
if [ -z "$PY" ]; then
  PY="${PY_FALLBACK:-}"
  [ -n "$PY" ] && echo "  [aviso] nenhum Python com as dependencias; tentando $PY"
fi
if [ -z "$PY" ]; then
  echo "[ERRO] Python nao encontrado. Instale o Miniforge."
  exit 1
fi
echo "  Python: $PY"

# --- 3. Garante %TEMP% -------------------------------------------------------
# O win32com grava o cache do gen_py em %TEMP%. Num shell MSYS2 essa variavel
# nao existe, ele cai no default C:\WINDOWS\gen_py e a simulacao morre com
# PermissionError antes de abrir o HEC-RAS.
if [ -z "${TEMP:-}" ] || [ -z "${TMP:-}" ]; then
  T="/c/Users/$QUEM/AppData/Local/Temp"
  [ -d "$T" ] || T="/tmp"
  TEMP="$(cygpath -w "$T" 2>/dev/null || echo "$T")"
  export TEMP TMP="$TEMP"
  echo "  TEMP: $TEMP"
fi

# --- 4. Poe as bibliotecas do Python escolhido na frente ---------------------
PYDIR="$(dirname "$PY")"
export PATH="$PYDIR:$PYDIR/Library/bin:$PYDIR/Library/mingw-w64/bin:$PYDIR/Scripts:$PATH"

# --- 3. Dependencias --------------------------------------------------------
if ! "$PY" -c "import numpy, geopandas, rasterio, shapely, pyproj, h5py, scipy" 2>/dev/null; then
  echo "  Instalando dependencias que faltam..."
  "$PY" -m pip install --quiet numpy geopandas rasterio shapely pyproj h5py scipy
fi
if ! "$PY" -c "import win32com.client" 2>/dev/null; then
  echo "  (pywin32 ausente -- necessario apenas para simular no HEC-RAS)"
  "$PY" -m pip install --quiet pywin32 2>/dev/null || true
fi

# --- 4. Roda ----------------------------------------------------------------
"$PY" rodar_modelo.py "$@"
RC=$?

[ "$RC" -eq 0 ] || echo "
  [ERRO] O pipeline terminou com erro (codigo $RC)."

# A dica vale mesmo quando o HEC-RAS nao fecha: o motor proprio roda sempre e as
# paginas sao geradas de qualquer jeito. Ela so aparecia em caso de sucesso, e
# apontava para app/ enquanto as paginas saiam na raiz -- ou seja, mandava
# servir uma pasta onde nada do que fora gerado estava.
echo "
  Para ver no navegador:
    $PY -m http.server 8050 --directory app

  e abra:
    http://localhost:8050/                              interface completa
    http://localhost:8050/<PROJETO>_cheia_motor.html    perfil animado (motor)
    http://localhost:8050/<PROJETO>_cheia_hecras.html   idem, pelo HEC-RAS
"
exit $RC
