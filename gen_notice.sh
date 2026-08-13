#!/usr/bin/env bash
#
# Build THIRD_PARTY_LICENSES.md from an installed pixi environment.
#
# License text is only ever COPIED FROM DISK -- never downloaded, never
# generated.
#
# Layout:
#   licenses/manual/<pkg>/...      hand-written text.
#   licenses/generated/<pkg>/...   extracted from the env and regenerated each time.
#   licenses/PREAMBLE.md          project's own license and any notes about copyleft components.
#   licenses/INVENTORY.txt         package/version/license table.
#   licenses/MISSING.txt           packages with no license text anywhere.
#   THIRD_PARTY_LICENSES.md        the final notice. 

set -euo pipefail

ENVIRONMENT=default
PREFIX=""
CHECK=0
OUT=THIRD_PARTY_LICENSES.md
LICDIR=licenses
GEN="$LICDIR/generated"
MANUAL="$LICDIR/manual"

while [ $# -gt 0 ]; do
  case "$1" in
    --environment) ENVIRONMENT="$2"; shift 2 ;;
    --prefix)      PREFIX="$2"; shift 2 ;;
    --check)       CHECK=1; shift ;;
    -o|--output)   OUT="$2"; shift 2 ;;
    -h|--help)     sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

TAB=$(printf '\t')

# hashing
if command -v shasum >/dev/null 2>&1; then
  hash_file() { shasum -a 256 "$1" | cut -d' ' -f1; }
elif command -v sha256sum >/dev/null 2>&1; then
  hash_file() { sha256sum "$1" | cut -d' ' -f1; }
else
  hash_file() { wc -c < "$1" | tr -d ' '; } 
fi

#  locate the prefix
if [ -z "$PREFIX" ]; then
  if [ -d ".pixi/envs/$ENVIRONMENT" ]; then
    PREFIX=".pixi/envs/$ENVIRONMENT"
  elif command -v pixi >/dev/null 2>&1; then
    # Detached-environments config can move envs outside the project.
    PREFIX=$(pixi run --environment "$ENVIRONMENT" bash -c 'echo "$CONDA_PREFIX"' 2>/dev/null | tail -1)
  fi
fi

[ -n "$PREFIX" ] || { echo "ERROR: could not locate environment '$ENVIRONMENT'. Run 'pixi install'." >&2; exit 1; }
[ -d "$PREFIX/conda-meta" ] || { echo "ERROR: $PREFIX is not a conda/pixi env (no conda-meta/)." >&2; exit 1; }
[ -x "$PREFIX/bin/python" ] || { echo "ERROR: no python in $PREFIX/bin." >&2; exit 1; }

SP=$("$PREFIX/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')
[ -d "$SP" ] || { echo "ERROR: site-packages not found: $SP" >&2; exit 1; }

# The resolved site-packages must live inside the prefix. If it does not, we are
# reading some other interpreter (an activated conda base, a system python, a
# symlink) and the notice would describe the wrong environment.
ABS_PREFIX=$(cd "$PREFIX" && pwd -P)
ABS_SP=$(cd "$SP" && pwd -P)
case "$ABS_SP" in
  "$ABS_PREFIX"/*) : ;;
  *) echo "ERROR: site-packages is outside the environment." >&2
     echo "  prefix:        $ABS_PREFIX" >&2
     echo "  site-packages: $ABS_SP" >&2
     echo "Deactivate any active conda env and re-run." >&2
     exit 1 ;;
esac

echo "environment:   $ENVIRONMENT"
echo "prefix:        $PREFIX"
echo "site-packages: $SP"

# reset dirs
rm -rf "$GEN"
mkdir -p "$GEN" "$MANUAL" "$LICDIR"
: > "$LICDIR/MISSING.txt"
: > "$LICDIR/INVENTORY.txt"

# Canonical package key: lowercase, underscores to hyphens, no version.
canon() { echo "$1" | tr '[:upper:]_' '[:lower:]-'; }

# Copy a file into the package dir unless identical content is already there.
copy_dedup() {
  src="$1"; destdir="$2"; relname="$3"
  mkdir -p "$destdir"
  h=$(hash_file "$src")
  for existing in "$destdir"/*; do
    [ -f "$existing" ] || continue
    if [ "$(hash_file "$existing")" = "$h" ]; then return 0; fi
  done
  target="$destdir/$relname"
  mkdir -p "$(dirname "$target")"
  if [ -e "$target" ]; then target="$target.2"; fi
  cp "$src" "$target"
}

#conda inventory
"$PREFIX/bin/python" - "$PREFIX" <<'PY' > "$LICDIR/conda-inventory.tsv"
import json, glob, sys
files = sorted(glob.glob(sys.argv[1] + '/conda-meta/*.json'))
if not files:
    sys.exit('ERROR: no conda-meta/*.json found')
for f in files:
    try:
        r = json.load(open(f))
    except Exception as exc:
        print('WARN: unreadable %s: %s' % (f, exc), file=sys.stderr)
        continue
    print('\t'.join([r.get('name', ''), r.get('version', ''),
                     r.get('license') or 'UNKNOWN',
                     r.get('extracted_package_dir') or '']))
PY

conda_total=$(wc -l < "$LICDIR/conda-inventory.tsv" | tr -d ' ')
echo "conda records: $conda_total"

#conda license text
conda_ok=0
while IFS="$TAB" read -r name ver lic dir; do
  [ -n "$name" ] || continue
  key=$(canon "$name")
  printf 'conda\t%s\t%s\t%s\n' "$name" "$ver" "$lic" >> "$LICDIR/INVENTORY.txt"
  if [ -n "$dir" ] && [ -d "$dir/info/licenses" ]; then
    found=$(find "$dir/info/licenses" -type f 2>/dev/null || true)
    if [ -n "$found" ]; then
      echo "$found" | while IFS= read -r f; do
        copy_dedup "$f" "$GEN/$key" "$(basename "$f")"
      done
      conda_ok=$((conda_ok + 1))
      continue
    fi
  fi
  printf '%s\t%s\t%s\tNO TEXT (conda)\n' "$name" "$ver" "$lic" >> "$LICDIR/MISSING.txt"
done < "$LICDIR/conda-inventory.tsv"

# pypi license text
for d in "$SP"/*.dist-info "$SP"/*.egg-info; do
  [ -d "$d" ] || continue
  base=$(basename "$d")
  stem=${base%.dist-info}; stem=${stem%.egg-info}
  # egg-info dirs are name-version-pyX.Y; dist-info are name-version.
  stem=${stem%-py3.[0-9]}; stem=${stem%-py3.[0-9][0-9]}
  name=${stem%%-*}; ver=${stem#*-}
  [ "$ver" = "$name" ] && ver=""
  key=$(canon "$name")
  printf 'pypi\t%s\t%s\t-\n' "$name" "$ver" >> "$LICDIR/INVENTORY.txt"
  found=$(find "$d" -type f \( -iname 'LICEN[CS]E*' -o -iname 'COPYING*' \
                            -o -iname 'NOTICE*' -o -iname 'AUTHORS*' \) 2>/dev/null || true)
  if [ -n "$found" ]; then
    echo "$found" | while IFS= read -r f; do
      copy_dedup "$f" "$GEN/$key" "$(basename "$f")"
    done
  else
    # Only a gap if the conda half did not already supply text for this package.
    if [ ! -d "$GEN/$key" ] && [ ! -d "$MANUAL/$key" ]; then
      printf '%s\t%s\t-\tNO TEXT (pypi)\n' "$name" "$ver" >> "$LICDIR/MISSING.txt"
    fi
  fi
done

sort -u "$LICDIR/INVENTORY.txt" -o "$LICDIR/INVENTORY.txt"
sort -u "$LICDIR/MISSING.txt"   -o "$LICDIR/MISSING.txt"

if [ -s "$LICDIR/MISSING.txt" ]; then
  tmp=$(mktemp)
  while IFS="$TAB" read -r name rest; do
    [ -n "$name" ] || continue
    key=$(canon "$name")
    if [ -d "$MANUAL/$key" ] && [ -n "$(ls -A "$MANUAL/$key" 2>/dev/null)" ]; then continue; fi
    printf '%s\t%s\n' "$name" "$rest" >> "$tmp"
  done < "$LICDIR/MISSING.txt"
  mv "$tmp" "$LICDIR/MISSING.txt"
fi

NEW=$(mktemp)

if [ -f "$LICDIR/PREAMBLE.md" ]; then
  cat "$LICDIR/PREAMBLE.md" >> "$NEW"
  printf '\n' >> "$NEW"
else
  cat >> "$NEW" <<'HDR'
# Third Party Notice

This file lists the licenses of third-party software distributed with this
project. License text is copied verbatim from the installed packages and is
never generated or rewritten.

Create `licenses/PREAMBLE.md` to replace this header with your own -- name
your project's own license there, and call out the copyleft components you
bundle. Markdown is allowed in that file.

HDR
fi

# write the packages' license text verbatim, each in a fenced code block.
emit_section() {
  srcroot="$1"; label="$2"
  [ -d "$srcroot" ] || return 0
  any=0
  for pkgdir in "$srcroot"/*/; do [ -d "$pkgdir" ] && any=1 && break; done
  [ "$any" = 1 ] || return 0
  printf '## %s\n\n' "$label" >> "$NEW"
  for pkgdir in "$srcroot"/*/; do
    [ -d "$pkgdir" ] || continue
    pkg=$(basename "$pkgdir")
    printf '### %s\n\n' "$pkg" >> "$NEW"
    meta=$(grep -i "$TAB$pkg$TAB" "$LICDIR/INVENTORY.txt" 2>/dev/null | head -1 || true)
    if [ -n "$meta" ]; then
      kind=$(printf '%s' "$meta" | cut -f1)
      ver=$(printf '%s' "$meta" | cut -f3)
      lic=$(printf '%s' "$meta" | cut -f4)
      [ -n "$ver" ] && printf -- '- **Version:** %s\n' "$ver" >> "$NEW"
      [ "$lic" != "-" ] && [ -n "$lic" ] && printf -- '- **License:** %s\n' "$lic" >> "$NEW"
      printf -- '- **Source:** %s\n' "$kind" >> "$NEW"
    fi
    printf '\n' >> "$NEW"
    find "$pkgdir" -type f | sort | while IFS= read -r f; do
      rel="${f#"$pkgdir"}"
      printf '#### `%s`\n\n' "$rel" >> "$NEW"
      # widen the fence past any backtick run present in the file
      fence='```'
      while grep -qF "$fence" "$f" 2>/dev/null; do fence="$fence"'`'; done
      printf '%s text\n' "$fence" >> "$NEW"
      cat "$f" >> "$NEW"
      printf '\n%s\n\n' "$fence" >> "$NEW"
    done
  done
}

emit_section "$GEN"    "Packages from the installed environment"
emit_section "$MANUAL" "Packages documented manually"

gen_count=$(find "$GEN" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
man_count=$(find "$MANUAL" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
gap_count=$(wc -l < "$LICDIR/MISSING.txt" | tr -d ' ')

echo "packages with text (extracted): $gen_count"
echo "packages with text (manual):    $man_count"
echo "gaps remaining:                 $gap_count"

if [ "$CHECK" = 1 ]; then
  if [ ! -f "$OUT" ] || ! cmp -s "$NEW" "$OUT"; then
    rm -f "$NEW"
    echo "ERROR: $OUT is out of date. Regenerate it." >&2
    exit 1
  fi
  rm -f "$NEW"
  echo "$OUT is up to date."
else
  mv "$NEW" "$OUT"
  echo "wrote $OUT"
fi

if [ "$gap_count" -gt 0 ]; then
  echo
  echo "Packages shipping no license text -- add files under $MANUAL/<name>/:"
  sed 's/^/  /' "$LICDIR/MISSING.txt"
  exit 3
fi