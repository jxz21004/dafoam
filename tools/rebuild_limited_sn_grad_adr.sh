#!/usr/bin/env bash
# Apply an idempotent, layout-preserving OpenFOAM-AD fix and rebuild ADR.
set -eo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
foam=${OPENFOAM_AD_ROOT:-/home/dafoamuser/dafoam/OpenFOAM/OpenFOAM-AD}
patchfile="$repo/tools/patches/openfoam-ad-limitedSnGrad-construction.patch"
backup=${SA_PATCH_BACKUP:-/mdouser/test_aerostructure/reports/openfoam_ad_limited_fix_backup}
mkdir -p "$backup"
header=src/finiteVolume/finiteVolume/snGradSchemes/limitedSnGrad/limitedSnGrad.H
lib="$foam/platforms/linux64GccDPInt32OptADR/lib/libfiniteVolumeADR.so"
if patch -d "$foam" -p1 --dry-run --forward < "$patchfile" >/dev/null 2>&1; then
 cp -n "$foam/$header" "$backup/limitedSnGrad.H.original"
 cp -n "$lib" "$backup/libfiniteVolumeADR.so.original"
 patch -d "$foam" -p1 < "$patchfile"
elif patch -d "$foam" -p1 --dry-run --reverse < "$patchfile" >/dev/null 2>&1; then
 echo 'limitedSnGrad fix already applied'
else
 echo 'Source does not match the patch; stopping' >&2; exit 1
fi
set +e
source /home/dafoamuser/dafoam/loadDAFoam.sh
source "$foam/etc/bashrc" WM_AD_MODE=ADR
set -e
export WM_AD_MODE=ADR WM_OPTIONS=linux64GccDPInt32OptADR
export FOAM_LIBBIN="$foam/platforms/$WM_OPTIONS/lib"
export LD_LIBRARY_PATH="$FOAM_LIBBIN:$FOAM_LIBBIN/sys-openmpi:$FOAM_LIBBIN/dummy:$LD_LIBRARY_PATH"
wmake -j "${BUILD_JOBS:-4}" objects "$foam/src/finiteVolume"
# Installed OpenFOAM-AD libraries are suffixed ADR, whereas upstream Make
# files still use unsuffixed names. Override only the final link target/libs.
(
 cd "$foam/src/finiteVolume"
 make -f "$WM_DIR/makefiles/general" MAKE_DIR=Make \
   "OBJECTS_DIR=$foam/build/$WM_OPTIONS/src/finiteVolume" \
   "LIB=$FOAM_LIBBIN/libfiniteVolumeADR" \
   "LIB_LIBS=-lOpenFOAMADR -lfileFormatsADR -lsurfMeshADR -lmeshToolsADR" libso
)
sha256sum "$foam/$header" "$lib" > "$backup/installed.sha256"
echo 'FINITE_VOLUME_ADR_REBUILD_COMPLETE'
