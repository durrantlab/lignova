#!/usr/bin/env bash
#
# Build THIRD_PARTY_LICENSES.md from installed pixi environments.
# License text is only ever COPIED FROM DISK -- never downloaded, never
# generated.  
#
# Layout:
#   licenses/manual/<pkg>/...      hand-written text.
#   licenses/generated/<pkg>/...   extracted from the env and regenerated each time.
#   docs/licensing-preamble.md         project's own license and any notes about copyleft components.
#   licenses/INVENTORY.txt         package/version/license table.
#   licenses/MISSING.txt           packages with no license text anywhere.
#   THIRD_PARTY_LICENSES.md        the final notice.
set -euo pipefail

ENVS=()
PREFIXES=()
CHECK=0
VERIFY=1
SELF=lignova
OUT=THIRD_PARTY_LICENSES.md
LICDIR=licenses
GEN="$LICDIR/generated"
MANUAL="$LICDIR/manual"

while [ $# -gt 0 ]; do
  case "$1" in
    --environment) ENVS+=("$2"); shift 2 ;;
    --prefix)      PREFIXES+=("$2"); shift 2 ;;
    --check)       CHECK=1; shift ;;
    --no-verify)   VERIFY=0; shift ;;
    --self)        SELF="$2"; shift 2 ;;
    -o|--output)   OUT="$2"; shift 2 ;;
    -h|--help)     sed -n '2,48p' "$0"; exit 0 ;;
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

if [ "${#ENVS[@]}" -eq 0 ]; then
  ENVS=(default)
fi
if [ "${#PREFIXES[@]}" -gt 0 ] && [ "${#PREFIXES[@]}" -ne "${#ENVS[@]}" ]; then
  echo "ERROR: pass one --environment per --prefix (got ${#ENVS[@]} and ${#PREFIXES[@]})." >&2
  exit 2
fi

# resolve_prefix <environment> <explicit-prefix-or-empty> -> prints the prefix
resolve_prefix() {
  env="$1"; explicit="${2:-}"
  prefix="$explicit"
  if [ -z "$prefix" ]; then
    if [ -d ".pixi/envs/$env" ]; then
      prefix=".pixi/envs/$env"
    elif command -v pixi >/dev/null 2>&1; then
      # Detached-environments config can move envs outside the project.
      prefix=$(pixi run --environment "$env" bash -c 'echo "$CONDA_PREFIX"' 2>/dev/null | tail -1)
    fi
  fi
  [ -n "$prefix" ] || { echo "ERROR: could not locate environment '$env'. Run 'pixi install --environment $env'." >&2; exit 1; }
  [ -d "$prefix/conda-meta" ] || { echo "ERROR: $prefix is not a conda/pixi env (no conda-meta/)." >&2; exit 1; }
  echo "$prefix"
}

# site_packages_of <prefix> -> prints purelib, after checking it is inside the env
site_packages_of() {
  prefix="$1"
  [ -x "$prefix/bin/python" ] || { echo "" ; return 0; }
  sp=$("$prefix/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])' 2>/dev/null || true)
  [ -n "$sp" ] || { echo "" ; return 0; }
  [ -d "$sp" ] || { echo "ERROR: site-packages not found: $sp" >&2; exit 1; }
  # The resolved site-packages must live inside the prefix. 
  abs_prefix=$(cd "$prefix" && pwd -P)
  abs_sp=$(cd "$sp" && pwd -P)
  case "$abs_sp" in
    "$abs_prefix"/*) : ;;
    *) echo "ERROR: site-packages is outside the environment." >&2
       echo "  prefix:        $abs_prefix" >&2
       echo "  site-packages: $abs_sp" >&2
       echo "Deactivate any active conda env and re-run." >&2
       exit 1 ;;
  esac
  echo "$sp"
}

pick_tool_python() {
  for cand in "$@"; do
    [ -n "$cand" ] || continue
    if "$cand" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 6) else 1)' \
        >/dev/null 2>&1; then
      echo "$cand"
      return 0
    fi
  done
  return 1
}

echo "environments:  ${ENVS[*]}"

# reset dirs 
rm -rf "$GEN"
mkdir -p "$GEN" "$MANUAL" "$LICDIR"
: > "$LICDIR/MISSING.txt"
: > "$LICDIR/index.tsv"

if [ ! -f "$LICDIR/OVERRIDES.tsv" ]; then
  cat > "$LICDIR/OVERRIDES.tsv" <<'OVR'
# Hand-held license decisions. Tab-separated:
#
#   <package-key>\t<SPDX expression>\t<reason>
#
# The key is the lowercased package name with underscores turned into hyphens.
# An entry here beats both the channel metadata and the license text detected on
# disk, and its reason is printed in THIRD_PARTY_LICENSES.md. Use it for stale
# metadata and for elections under a dual license.
OVR
fi

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

# Inventory every environment first, tagging each record with the environment it came from.
: > "$LICDIR/conda-inventory.tsv"
: > "$LICDIR/pypi-inventory.tsv"

# Resolve every prefix first, so the helper interpreter can be chosen from all of
# them before any of them is read.
PRE=()
env_index=0
for env in "${ENVS[@]}"; do
  explicit=""
  if [ "${#PREFIXES[@]}" -gt 0 ]; then explicit="${PREFIXES[$env_index]}"; fi
  env_index=$((env_index + 1))
  prefix=$(resolve_prefix "$env" "$explicit")
  PRE+=("$prefix")
  echo "  $env -> $prefix"
done

CANDIDATES=()
for prefix in "${PRE[@]}"; do CANDIDATES+=("$prefix/bin/python"); done
CANDIDATES+=(python3 python)
TOOL_PYTHON=$(pick_tool_python "${CANDIDATES[@]}") || {
  echo "ERROR: no Python 3.6+ interpreter available to run this script's helpers." >&2
  echo "Tried: ${CANDIDATES[*]}" >&2
  exit 1
}
echo "helper python: $TOOL_PYTHON"

env_index=0
for env in "${ENVS[@]}"; do
  prefix="${PRE[$env_index]}"
  env_index=$((env_index + 1))
  sp=$(site_packages_of "$prefix")

  #conda inventory
  "$TOOL_PYTHON" - "$prefix" <<'PY' > "$LICDIR/.conda.tsv"
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
  awk -F"$TAB" -v e="$env" 'BEGIN { OFS = "\t" } NF { print e, $0 }' \
    "$LICDIR/.conda.tsv" >> "$LICDIR/conda-inventory.tsv"

  # pypi inventory. Unlike conda-meta, a dist-info carries the license in its
  # METADATA, so read it instead of leaving the column blank. 
  if [ -z "$sp" ]; then
    echo "  note: $env has no usable python; skipping its site-packages" >&2
    : > "$LICDIR/.pypi.tsv"
  else
  "$TOOL_PYTHON" - "$sp" "$SELF" <<'PY' > "$LICDIR/.pypi.tsv"
import os, re, sys

sp, selfname = sys.argv[1], sys.argv[2]

def canon(s):
    return s.strip().lower().replace('_', '-')

# Only the classifiers that name an unambiguous license. "BSD License" and
# "GNU General Public License (GPL)" name a family, not a license, so they are
# left as-is and the audit will not try to check them.
CLASSIFIERS = {
    'mit license': 'MIT',
    'mit no attribution license (mit-0)': 'MIT-0',
    'apache software license': 'Apache-2.0',
    'bsd license': 'BSD',
    'isc license (iscl)': 'ISC',
    'mozilla public license 2.0 (mpl 2.0)': 'MPL-2.0',
    'mozilla public license 1.1 (mpl 1.1)': 'MPL-1.1',
    'gnu general public license v2 (gplv2)': 'GPL-2.0-only',
    'gnu general public license v2 or later (gplv2+)': 'GPL-2.0-or-later',
    'gnu general public license v3 (gplv3)': 'GPL-3.0-only',
    'gnu general public license v3 or later (gplv3+)': 'GPL-3.0-or-later',
    'gnu lesser general public license v2 (lgplv2)': 'LGPL-2.0-only',
    'gnu lesser general public license v2 or later (lgplv2+)': 'LGPL-2.0-or-later',
    'gnu lesser general public license v3 (lgplv3)': 'LGPL-3.0-only',
    'gnu lesser general public license v3 or later (lgplv3+)': 'LGPL-3.0-or-later',
    'python software foundation license': 'PSF-2.0',
    'the unlicense (unlicense)': 'Unlicense',
    'zlib/libpng license': 'Zlib',
    'boost software license 1.0 (bsl-1.0)': 'BSL-1.0',
}
# free-text License: values worth normalizing; anything else is passed through.
LOOSE = {
    'mit': 'MIT', 'mit license': 'MIT', 'expat': 'MIT',
    'apache 2.0': 'Apache-2.0', 'apache-2.0': 'Apache-2.0',
    'apache license 2.0': 'Apache-2.0', 'apache software license': 'Apache-2.0',
    'bsd-3-clause': 'BSD-3-Clause', 'bsd 3-clause': 'BSD-3-Clause',
    'bsd-2-clause': 'BSD-2-Clause', 'bsd 2-clause': 'BSD-2-Clause',
    'lgplv2.1': 'LGPL-2.1-only', 'lgpl-2.1': 'LGPL-2.1-only',
    'lgplv3': 'LGPL-3.0-only', 'lgpl-3.0': 'LGPL-3.0-only',
    'gplv2': 'GPL-2.0-only', 'gplv3': 'GPL-3.0-only',
    'mpl-2.0': 'MPL-2.0', 'mpl 2.0': 'MPL-2.0',
    'psf': 'PSF-2.0', 'python software foundation license': 'PSF-2.0',
    'isc': 'ISC', 'unlicense': 'Unlicense', 'zlib': 'Zlib',
}

def headers(path):
    """Parse the RFC822 header block of METADATA/PKG-INFO."""
    out = []
    try:
        fh = open(path, encoding='utf-8', errors='replace')
    except OSError:
        return out
    with fh:
        key = None
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                break
            if line[:1] in (' ', '\t') and key:
                out[-1] = (key, out[-1][1] + ' ' + line.strip())
                continue
            if ':' not in line:
                continue
            key, val = line.split(':', 1)
            key = key.strip().lower()
            out.append((key, val.strip()))
    return out

def license_of(hdrs):
    got = dict()
    ids = []
    for k, v in hdrs:
        got.setdefault(k, v)
        if k == 'classifier' and v.lower().startswith('license ::'):
            leaf = v.split('::')[-1].strip().lower()
            if leaf in CLASSIFIERS:
                ids.append(CLASSIFIERS[leaf])
    # PEP 639: a real SPDX expression.
    expr = got.get('license-expression', '').strip()
    if expr:
        return expr
    if ids:
        seen = []
        for i in ids:
            if i not in seen:
                seen.append(i)
        return ' AND '.join(seen)
    # Legacy License: field. Some projects paste the whole license in here;
    # those are useless as identifiers, so keep only short values.
    loose = got.get('license', '').strip()
    if loose and len(loose) <= 40 and '\n' not in loose:
        return LOOSE.get(loose.lower(), loose)
    return ''

for entry in sorted(os.listdir(sp)):
    if not (entry.endswith('.dist-info') or entry.endswith('.egg-info')):
        continue
    d = os.path.join(sp, entry)
    if not os.path.isdir(d):
        continue
    stem = re.sub(r'\.(dist|egg)-info$', '', entry)
    stem = re.sub(r'-py3\.\d+$', '', stem)
    name, _, ver = stem.partition('-')
    hdrs = []
    for meta in ('METADATA', 'PKG-INFO'):
        if os.path.isfile(os.path.join(d, meta)):
            hdrs = headers(os.path.join(d, meta))
            break
    got = {}
    for k, v in hdrs:
        got.setdefault(k, v)
    name = got.get('name', name)
    ver = got.get('version', ver)
    if canon(name) == canon(selfname):
        continue
    # '-' rather than '': a tab-split empty field collapses in the shell loop.
    print('\t'.join([name, ver, license_of(hdrs) or '-', d]))
PY
  fi
  awk -F"$TAB" -v e="$env" 'BEGIN { OFS = "\t" } NF { print e, $0 }' \
    "$LICDIR/.pypi.tsv" >> "$LICDIR/pypi-inventory.tsv"
done
rm -f "$LICDIR/.conda.tsv" "$LICDIR/.pypi.tsv"

conda_total=$(wc -l < "$LICDIR/conda-inventory.tsv" | tr -d ' ')
echo "conda records: $conda_total"

#conda license text
while IFS="$TAB" read -r env name ver lic dir; do
  [ -n "$name" ] || continue
  key=$(canon "$name")
  [ "$key" = "$(canon "$SELF")" ] && continue
  printf 'conda\t%s\t%s\t%s\t%s\t%s\n' "$key" "$name" "$ver" "$lic" "$env" >> "$LICDIR/index.tsv"
  if [ -n "$dir" ] && [ -d "$dir/info/licenses" ]; then
    found=$(find "$dir/info/licenses" -type f 2>/dev/null || true)
    if [ -n "$found" ]; then
      echo "$found" | while IFS= read -r f; do
        copy_dedup "$f" "$GEN/$key" "$(basename "$f")"
      done
      continue
    fi
  fi
  printf '%s\t%s\t%s\tNO TEXT (conda, %s)\n' "$name" "$ver" "$lic" "$env" >> "$LICDIR/MISSING.txt"
done < "$LICDIR/conda-inventory.tsv"


# pypi license text
while IFS="$TAB" read -r env name ver lic d; do
  [ -n "$name" ] || continue
  key=$(canon "$name")
  printf 'pypi\t%s\t%s\t%s\t%s\t%s\n' "$key" "$name" "$ver" "$lic" "$env" >> "$LICDIR/index.tsv"
  found=$(find "$d" -type f \( -iname 'LICEN[CS]E*' -o -iname 'COPYING*' \
                            -o -iname 'NOTICE*' -o -iname 'AUTHORS*' \) 2>/dev/null || true)
  if [ -n "$found" ]; then
    echo "$found" | while IFS= read -r f; do
      copy_dedup "$f" "$GEN/$key" "$(basename "$f")"
    done
  else
    # Only a gap if the conda half did not already supply text for this package.
    if [ ! -d "$GEN/$key" ] && [ ! -d "$MANUAL/$key" ]; then
      printf '%s\t%s\t%s\tNO TEXT (pypi, %s)\n' "$name" "$ver" "$lic" "$env" >> "$LICDIR/MISSING.txt"
    fi
  fi
done < "$LICDIR/pypi-inventory.tsv"

sort -u "$LICDIR/MISSING.txt" -o "$LICDIR/MISSING.txt"

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

# Reads licenses/index.tsv, licenses/manual/*/META.txt and licenses/OVERRIDES.tsv,
# fingerprints the copied text, and writes:
#   licenses/INVENTORY.txt       the table, with effective licenses
#   licenses/LICENSE-AUDIT.txt   every decision that was not a straight pass-through
#   licenses/meta.tsv            per-package metadata for the emitter below
#   licenses/CONFLICTS.txt       unresolved disagreements (exit 4)
#   licenses/PLACEHOLDERS.txt    unfilled manual scaffolding (exit 5)
"$TOOL_PYTHON" - "$LICDIR" "$GEN" "$MANUAL" "$VERIFY" <<'PY'
import os, sys

LICDIR, GEN, MANUAL, VERIFY = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == '1'

PLACEHOLDER = 'REPLACE-WITH-VERBATIM-LICENSE-TEXT'
SIDECARS = {'META.txt', 'NOTES.md'}


def canon(s):
    return s.strip().lower().replace('_', '-')


def read(path):
    try:
        with open(path, 'rb') as fh:
            return fh.read().decode('utf-8', 'replace')
    except OSError:
        return ''


def flat(text):
    return ' '.join(text.split()).lower()


def detect(text):
    """Identify well-known license texts by a distinctive phrase from each."""
    t = flat(text)
    ids = set()
    if 'mozilla public license version 2.0' in t:
        ids.add('MPL-2.0')
    if 'mozilla public license version 1.1' in t:
        ids.add('MPL-1.1')
    if 'apache license' in t and 'version 2.0, january 2004' in t:
        ids.add('Apache-2.0')
    if 'gnu lesser general public license version 3, 29 june 2007' in t:
        ids.add('LGPL-3.0')
    if 'gnu lesser general public license version 2.1, february 1999' in t:
        ids.add('LGPL-2.1')
    if 'gnu library general public license version 2, june 1991' in t:
        ids.add('LGPL-2.0')
    if 'gnu general public license version 3, 29 june 2007' in t:
        ids.add('GPL-3.0')
    if 'gnu general public license version 2, june 1991' in t:
        ids.add('GPL-2.0')
    if 'gnu affero general public license version 3, 19 november 2007' in t:
        ids.add('AGPL-3.0')
    if 'permission is hereby granted, free of charge' in t:
        ids.add('MIT')
    if 'permission to use, copy, modify, and/or distribute this software for any purpose with or without fee' in t:
        ids.add('ISC' if 'appear in all copies' in t else '0BSD')
    if ('permission to use, copy, modify, and distribute this software and its '
            'documentation for any purpose and without fee is hereby granted') in t:
        ids.add('HPND')
    if 'redistribution and use in source and binary forms' in t:
        if 'all advertising materials mentioning features or use' in t:
            ids.add('BSD-4-Clause')
        elif 'neither the name' in t or 'nor the names of its contributors' in t:
            ids.add('BSD-3-Clause')
        else:
            ids.add('BSD-2-Clause')
    if 'altered source versions must be plainly marked as such' in t:
        ids.add('Zlib')
    if 'boost software license' in t:
        ids.add('BSL-1.0')
    if 'sil open font license' in t:
        ids.add('OFL-1.1')
    if 'python software foundation license' in t:
        ids.add('PSF-2.0')
    if 'free and unencumbered software released into the public domain' in t:
        ids.add('CC0-1.0' if 'creative commons' in t else 'Unlicense')
    return ids


# Everything detect() can emit; a declared identifier outside this set cannot be
# contradicted by the fingerprints, so it is passed through untouched.
VERIFIABLE = {
    'mpl-2.0', 'mpl-1.1', 'apache-2.0', 'lgpl-3.0', 'lgpl-2.1', 'lgpl-2.0',
    'gpl-3.0', 'gpl-2.0', 'agpl-3.0', 'mit', 'isc', '0bsd', 'hpnd',
    'bsd-4-clause', 'bsd-3-clause', 'bsd-2-clause', 'zlib', 'bsl-1.0',
    'ofl-1.1', 'psf-2.0', 'cc0-1.0', 'unlicense',
}
ALIAS = {'python-2.0': 'psf-2.0', 'zlib-acknowledgement': 'zlib',
         'mit-0': 'mit', 'x11': 'hpnd'}
# Family names, not licenses: a PyPI "BSD License" classifier says nothing about
# which BSD variant applies, so such a claim carries no information here.
AMBIGUOUS = {'bsd', 'gpl', 'lgpl', 'agpl', 'apache', 'mpl', 'cc', 'python',
             'unknown', 'other/proprietary', 'dfsg-approved', 'osi-approved'}
# "X-or-later" also permits the later licenses of the same family.
LATER = {
    'lgpl-2.0': ['lgpl-2.0', 'lgpl-2.1', 'lgpl-3.0'],
    'lgpl-2.1': ['lgpl-2.1', 'lgpl-3.0'],
    'lgpl-3.0': ['lgpl-3.0'],
    'gpl-2.0': ['gpl-2.0', 'gpl-3.0'],
    'gpl-3.0': ['gpl-3.0'],
    'agpl-3.0': ['agpl-3.0'],
}
# License files that describe someone else's code, not the package's own terms.
NOT_EVIDENCE = ('bundled', 'third-party', 'third_party', 'thirdparty',
                '3rdparty', 'vendor', 'authors', 'notice')


def norm_id(tok):
    tok = tok.strip().strip('()').lower()
    for suffix in ('-only', '-or-later', '+'):
        if tok.endswith(suffix):
            tok = tok[: -len(suffix)]
    return ALIAS.get(tok, tok)


def declared_ids(expr):
    """SPDX expression -> (identifiers, identifiers the grant also permits)."""
    ids, allow = set(), set()
    if not expr:
        return ids, allow
    tokens = expr.replace('(', ' ').replace(')', ' ').split()
    skip_next = False
    for tok in tokens:
        low = tok.lower()
        if skip_next:            # the operand of WITH is an exception, not a license
            skip_next = False
            continue
        if low == 'with':
            skip_next = True
            continue
        if low in ('and', 'or'):
            continue
        base = norm_id(tok)
        if not base or base in ('-', 'none'):
            continue
        ids.add(base)
        if low.endswith('-or-later') or low.endswith('+'):
            allow.update(LATER.get(base, [base]))
        else:
            allow.add(base)
    return ids, allow


class Claim(object):
    """One party's statement of a package's license."""

    def __init__(self, source, expr):
        self.source = source
        self.expr = expr
        self.ids, self.allow = declared_ids(expr)
        self.ambiguous = (not self.ids) or bool(self.ids & AMBIGUOUS)
        # Verifiable = every identifier is one the fingerprints can recognize,
        # so silence from the text is meaningful. "curl" or "LicenseRef-Qhull"
        # are specific and correct but unverifiable; leave those alone.
        self.verifiable = bool(self.ids) and self.ids <= VERIFIABLE

    def matches(self, detected):
        return bool(self.allow & detected)


def pkg_files(root, key):
    d = os.path.join(root, key)
    out = []
    if not os.path.isdir(d):
        return out
    for dirpath, _dirs, files in os.walk(d):
        for f in sorted(files):
            full = os.path.join(dirpath, f)
            out.append((os.path.relpath(full, d), full))
    return sorted(out)


# ---- inputs ---------------------------------------------------------------
records = {}   # key -> dict


def rec(key):
    return records.setdefault(key, {
        'kinds': [], 'names': [], 'versions': [], 'declared': {},
        'meta': {}, 'note': '', 'status': 'OK', 'rows': [], 'envs': [],
    })


for line in read(os.path.join(LICDIR, 'index.tsv')).splitlines():
    if not line.strip():
        continue
    parts = line.split('\t')
    while len(parts) < 6:
        parts.append('')
    kind, key, name, ver, lic, env = parts[:6]
    r = rec(key)
    if kind not in r['kinds']:
        r['kinds'].append(kind)
    if name and name not in r['names']:
        r['names'].append(name)
    if ver and ver not in r['versions']:
        r['versions'].append(ver)
    if lic and lic not in ('-', 'UNKNOWN'):
        r['declared'].setdefault(kind, lic)
    if env and env not in r['envs']:
        r['envs'].append(env)
    r['rows'].append([kind, name, ver, env])

placeholders = []
if os.path.isdir(MANUAL):
    for key in sorted(os.listdir(MANUAL)):
        d = os.path.join(MANUAL, key)
        if not os.path.isdir(d) or not os.listdir(d):
            continue
        r = rec(key)
        if 'manual' not in r['kinds']:
            r['kinds'].append('manual')
        meta = {}
        for line in read(os.path.join(d, 'META.txt')).splitlines():
            line = line.strip()
            if not line or line.startswith('#') or ':' not in line:
                continue
            k, v = line.split(':', 1)
            meta[k.strip().lower()] = v.strip()
        r['meta'] = meta
        if meta.get('name') and meta['name'] not in r['names']:
            r['names'].insert(0, meta['name'])
        if meta.get('version') and meta['version'] not in r['versions']:
            r['versions'].append(meta['version'])
        if meta.get('license'):
            r['declared']['manual'] = meta['license']
        for rel, full in pkg_files(MANUAL, key):
            if os.path.basename(rel) in SIDECARS:
                continue
            if PLACEHOLDER in read(full):
                placeholders.append('%s\t%s' % (key, rel))

overrides = {}
for line in read(os.path.join(LICDIR, 'OVERRIDES.tsv')).splitlines():
    if not line.strip() or line.lstrip().startswith('#'):
        continue
    parts = line.split('\t')
    key = canon(parts[0])
    spdx = parts[1].strip() if len(parts) > 1 else ''
    why = parts[2].strip() if len(parts) > 2 else ''
    if key and spdx:
        overrides[key] = (spdx, why)
        rec(key)

# ---- decide the effective license -----------------------------------------
audit = []
conflicts = []

SOURCE_ORDER = ['manual', 'conda', 'pypi']

for key, r in records.items():
    claims = [Claim(s, r['declared'][s]) for s in SOURCE_ORDER if s in r['declared']]
    declared = claims[0].expr if claims else ''
    found = set()
    for root in (GEN, MANUAL):
        for rel, full in pkg_files(root, key):
            base = os.path.basename(rel).lower()
            if base in {s.lower() for s in SIDECARS}:
                continue
            if any(tag in base for tag in NOT_EVIDENCE):
                continue
            found |= detect(read(full))
    r['detected'] = sorted(found)
    effective = declared

    if key in overrides:
        effective, why = overrides[key]
        r['status'] = 'OVERRIDE'
        r['note'] = why or 'set by hand in licenses/OVERRIDES.tsv'
        if (declared and canon(declared) != canon(effective)
                and canon(declared) not in r['note'].lower()):
            r['note'] += ' (package metadata says %s)' % declared
    elif not VERIFY or not found:
        r['status'] = 'OK' if declared else 'NO-EVIDENCE'
    else:
        det = {norm_id(i) for i in found}
        agreeing = [c for c in claims if c.matches(det)]
        informative = [c for c in claims if not c.ambiguous]
        checkable = [c for c in informative if c.verifiable]
        if agreeing:
            effective = agreeing[0].expr
            r['status'] = 'OK'
            if agreeing[0] is not claims[0]:
                # e.g. conda says BSD-3-Clause, PyPI says MIT, the text is MIT.
                r['note'] = ('license as declared by %s and confirmed by the '
                             'license text on file; %s claims %s'
                             % (agreeing[0].source, claims[0].source, claims[0].expr))
            elif agreeing[0].ids != det and agreeing[0].allow != agreeing[0].ids:
                # An "or later" grant, with the text pinning one version of it.
                r['status'] = 'RANGE'
        elif not informative:
            # Nothing declared, or only a family name such as "BSD License".
            if len(found) == 1:
                effective = r['detected'][0]
                r['status'] = 'FILLED'
                r['note'] = ('license taken from the license text on file; '
                             'the package metadata declares %s'
                             % (('only %s' % declared) if declared else 'none'))
            else:
                r['status'] = 'NO-DECLARATION'
        elif not checkable:
            # e.g. "IJG AND BSD-3-Clause AND Zlib", "curl", "LicenseRef-Qhull":
            # specific, but not something the fingerprints can speak to.
            r['status'] = 'UNVERIFIED'
        elif len(found) == 1:
            effective = r['detected'][0]
            r['status'] = 'CORRECTED'
            r['note'] = ('license taken from the license text on file; '
                         'the %s metadata claims %s'
                         % (checkable[0].source, checkable[0].expr))
        else:
            r['status'] = 'CONFLICT'
            conflicts.append('%s\tdeclared: %s\tdetected: %s'
                             % (key, declared or '(none)', ', '.join(r['detected'])))

    r['declared_str'] = declared
    r['effective'] = effective
    if r['status'] in ('OVERRIDE', 'FILLED', 'CORRECTED', 'CONFLICT', 'RANGE'):
        audit.append('%-24s %-10s declared=%-28s detected=%-28s -> %s'
                     % (key, r['status'], declared or '(none)',
                        ', '.join(r['detected']) or '(none)', effective or '(none)'))

# ---- outputs --------------------------------------------------------------
inv = {}
for key, r in records.items():
    lic = r['effective'] or '-'
    rows = [tuple(x) for x in r['rows']]
    if not rows:                           # manual-only entry
        rows = [('manual',
                 r['names'][0] if r['names'] else key,
                 r['versions'][0] if r['versions'] else '', '-')]
    for kind, name, ver, env in rows:
        envs = inv.setdefault((kind, name, ver, lic), [])
        if env and env not in envs:
            envs.append(env)
with open(os.path.join(LICDIR, 'INVENTORY.txt'), 'w') as fh:
    for (kind, name, ver, lic), envs in sorted(inv.items()):
        fh.write('\t'.join([kind, name, ver, lic,
                            ' '.join(envs) if envs else '-']) + '\n')

# A hand-written entry carrying only notes, with no license text from any
# environment, would otherwise pass silently -- that is how a package drops out
# of the notice when the environment holding it is not scanned.
missing_path = os.path.join(LICDIR, 'MISSING.txt')
gaps = [ln for ln in read(missing_path).splitlines() if ln.strip()]
for key in sorted(records):
    if not os.path.isdir(os.path.join(MANUAL, key)):
        continue
    texts = [rel for rel, _f in pkg_files(GEN, key)]
    texts += [rel for rel, _f in pkg_files(MANUAL, key)
              if os.path.basename(rel) not in SIDECARS]
    if texts:
        continue
    r = records[key]
    gaps.append('\t'.join([
        r['names'][0] if r['names'] else key,
        r['versions'][0] if r['versions'] else '',
        r['effective'] or '-',
        'NO TEXT (only notes under licenses/manual; add the text, or scan the '
        'environment that ships this package)']))
with open(missing_path, 'w') as fh:
    for line in sorted(set(gaps)):
        fh.write(line + '\n')

with open(os.path.join(LICDIR, 'LICENSE-AUDIT.txt'), 'w') as fh:
    if audit:
        fh.write('# Packages whose recorded license is not a straight copy of the\n'
                 '# channel metadata. Regenerated by gen_notice.sh.\n\n')
        for line in sorted(audit):
            fh.write(line + '\n')

with open(os.path.join(LICDIR, 'CONFLICTS.txt'), 'w') as fh:
    for line in sorted(conflicts):
        fh.write(line + '\n')

with open(os.path.join(LICDIR, 'PLACEHOLDERS.txt'), 'w') as fh:
    for line in sorted(placeholders):
        fh.write(line + '\n')

# key, kinds, display name, version, license, note, homepage
with open(os.path.join(LICDIR, 'meta.tsv'), 'w') as fh:
    for key in sorted(records):
        r = records[key]
        meta = r['meta']
        display = meta.get('name') or (r['names'][0] if r['names'] else key)
        source = meta.get('source') or ', '.join(r['kinds'])
        version = meta.get('version') or (r['versions'][0] if r['versions'] else '')
        note = r['note'] or meta.get('note', '')
        fh.write('\t'.join([key, source, display, version,
                            r['effective'] or '', note,
                            meta.get('homepage', ''),
                            ' '.join(r['envs'])]) + '\n')

counts = {}
for r in records.values():
    counts[r['status']] = counts.get(r['status'], 0) + 1
print('license audit:  ' + ', '.join('%s=%d' % kv for kv in sorted(counts.items())))
PY

if [ -s "$LICDIR/PLACEHOLDERS.txt" ]; then
  echo "ERROR: unfilled placeholders under $MANUAL -- paste the real license text:" >&2
  sed 's/^/  /' "$LICDIR/PLACEHOLDERS.txt" >&2
  exit 5
fi

NEW=$(mktemp)

if [ -f "docs/licensing-preamble.md" ]; then
  cat "docs/licensing-preamble.md" >> "$NEW"
  printf '\n' >> "$NEW"
else
  cat >> "$NEW" <<'HDR'
# Third Party Notice

This file lists the licenses of third-party software used by this project.
License text is copied verbatim from the installed packages and is never
generated or rewritten.

Create `docs/licensing-preamble.md` to replace this header with your own -- name
your project's own license there, say whether these packages are distributed
with the project or installed from upstream channels, and call out the
copyleft ones. Markdown is allowed in that file.

HDR
fi

# metadata for one package, from the audit's meta.tsv
meta_field() {
  awk -F"$TAB" -v k="$1" -v n="$2" '$1==k {print $n; exit}' "$LICDIR/meta.tsv"
}

# write one package's license text verbatim, each file in a fenced code block.
emit_files() {
  pkgdir="$1"; origin="$2"
  [ -d "$pkgdir" ] || return 0
  find "$pkgdir" -type f | sort | while IFS= read -r f; do
    rel="${f#"$pkgdir"/}"
    case "$(basename "$rel")" in META.txt|NOTES.md) continue ;; esac
    if [ -n "$origin" ]; then
      printf '#### `%s` %s\n\n' "$rel" "$origin" >> "$NEW"
    else
      printf '#### `%s`\n\n' "$rel" >> "$NEW"
    fi
    # widen the fence past any backtick run present in the file
    fence='```'
    while grep -qF "$fence" "$f" 2>/dev/null; do fence="$fence"'`'; done
    printf '%s text\n' "$fence" >> "$NEW"
    cat "$f" >> "$NEW"
    printf '\n%s\n\n' "$fence" >> "$NEW"
  done
}

# One section per package, whether its text came out of the environment, was
# written by hand, or both.
keys=$( { for d in "$GEN"/*/ "$MANUAL"/*/; do [ -d "$d" ] || continue; basename "$d"; done; } | LC_ALL=C sort -u )

if [ -n "$keys" ]; then
  printf '## Third-party packages\n\n' >> "$NEW"
  printf '%s\n' "$keys" | while IFS= read -r pkg; do
    [ -n "$pkg" ] || continue
    display=$(meta_field "$pkg" 3); [ -n "$display" ] || display="$pkg"
    source=$(meta_field "$pkg" 2)
    version=$(meta_field "$pkg" 4)
    lic=$(meta_field "$pkg" 5)
    note=$(meta_field "$pkg" 6)
    home=$(meta_field "$pkg" 7)
    envs=$(meta_field "$pkg" 8)
    printf '### %s\n\n' "$display" >> "$NEW"
    [ -n "$version" ] && printf -- '- **Version:** %s\n' "$version" >> "$NEW"
    [ -n "$lic" ] && printf -- '- **License:** %s\n' "$lic" >> "$NEW"
    [ -n "$source" ] && printf -- '- **Source:** %s\n' "$source" >> "$NEW"
    [ -n "$envs" ] && printf -- '- **Environment:** %s\n' "$(echo "$envs" | tr ' ' ',' | sed 's/,/, /g')" >> "$NEW"
    [ -n "$home" ] && printf -- '- **Homepage:** %s\n' "$home" >> "$NEW"
    [ -n "$note" ] && printf -- '- **License note:** %s\n' "$note" >> "$NEW"
    printf '\n' >> "$NEW"
    if [ -f "$MANUAL/$pkg/NOTES.md" ]; then
      cat "$MANUAL/$pkg/NOTES.md" >> "$NEW"
      printf '\n' >> "$NEW"
    fi
    emit_files "$GEN/$pkg" ''
    emit_files "$MANUAL/$pkg" '(added by hand)'
  done
fi

gen_count=$(find "$GEN" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
man_count=$(find "$MANUAL" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
gap_count=$(wc -l < "$LICDIR/MISSING.txt" | tr -d ' ')
conflict_count=$(wc -l < "$LICDIR/CONFLICTS.txt" | tr -d ' ')

echo "packages with text (extracted): $gen_count"
echo "packages with text (manual):    $man_count"
echo "gaps remaining:                 $gap_count"
echo "license conflicts:              $conflict_count"

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

if [ "$conflict_count" -gt 0 ]; then
  echo
  echo "License text disagrees with the metadata and more than one license was"
  echo "found, so the winner is not obvious. Decide each in $LICDIR/OVERRIDES.tsv:"
  sed 's/^/  /' "$LICDIR/CONFLICTS.txt"
  exit 4
fi

if [ "$gap_count" -gt 0 ]; then
  echo
  echo "Packages shipping no license text -- add files under $MANUAL/<name>/:"
  sed 's/^/  /' "$LICDIR/MISSING.txt"
  exit 3
fi