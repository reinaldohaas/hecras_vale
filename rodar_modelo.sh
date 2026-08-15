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

# --- 1. Localiza o Python ---------------------------------------------------
PY=""
for c in "$HOME/miniforge3/python.exe" "$HOME/miniconda3/python.exe" \
         "$HOME/anaconda3/python.exe" "$HOME/miniforge3/bin/python" \
         "$HOME/miniconda3/bin/python"; do
  [ -x "$c" ] && PY="$c" && break
done
[ -z "$PY" ] && PY="$(command -v python3 || command -v python || true)"
if [ -z "$PY" ]; then
  echo "[ERRO] Python nao encontrado. Instale o Miniforge."
  exit 1
fi
echo "  Python: $PY"

# --- 2. Blinda o PATH -------------------------------------------------------
# O MSYS2/Git Bash coloca /c/msys64/*/bin no PATH, e as DLLs de libpng e
# freetype de la se sobrepoem as do conda: o matplotlib quebra ao salvar
# figura e o rasterio derruba o processo. Tirando o msys64 do PATH e pondo o
# Python na frente, as bibliotecas corretas vencem.
PYDIR="$(dirname "$PY")"
PATH="$(echo "$PATH" | tr ':' '\n' | grep -viE '^/c/msys64|^/mingw|^/ucrt64' | paste -sd: -)"
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

if [ "$RC" -eq 0 ]; then
  echo
  echo "  Concluido. Para ver no navegador:"
  echo "    $PY -m http.server 8050 --directory app"
else
  echo
  echo "  [ERRO] O pipeline terminou com erro (codigo $RC)."
fi
exit $RC
