#!/bin/bash

usage() {
  if [ -n "$1" ]; then echo "ERROR: $1" >&2; echo >&2; fi
  cat >&2 <<'EOF'
Run every test_*.R script beside this one, serially, against an h0testr R/ source directory.
One Rscript process per test, since the DEqMS, msqrob2 and prolfqua engines cannot share a
session, and concurrent Rscript processes segfault on this platform.

Usage: run_all_tests.sh <log_dir> [r_dir] [rscript]

Required positional arguments:
  <log_dir>  Directory for output; created if absent. Existing *.log and summary.txt in it
               are deleted first.

Optional positional arguments:
  [r_dir]    Path to the h0testr package R/ source directory. Default: the R directory two
               levels above this script, i.e. the package this script ships with.
  [rscript]  Rscript executable. Default: /c/Programs/R/bin/Rscript

Output: <log_dir>/<test_name>.log per test, holding that test's stdout and stderr;
  <log_dir>/summary.txt, one line per test with its exit code, PASS and FAIL counts and
  elapsed seconds, ending in DONE. The same summary lines go to stdout as tests finish.

Exit codes: 0 every test exited 0; 1 one or more tests failed; 2 usage error; 3 r_dir not a
  directory; 4 log_dir not writable; 5 no test_*.R found.

Examples:
  run_all_tests.sh ./testlogs
  run_all_tests.sh /c/Users/me/data/h0testr_tests/1/testlogs
  run_all_tests.sh ./testlogs /c/other/h0testr/R /c/Programs/R-4.4.1/bin/Rscript
EOF
  exit 2
}

if [ $# -lt 1 ] || [ $# -gt 3 ]; then usage "wrong number of arguments"; fi
case "$1" in -h|--help) usage;; esac

TDIR=$(cd "$(dirname "$0")" && pwd)
LDIR=$1
RDIR=${2:-$(cd "$TDIR/../.." && pwd)/R}
RS=${3:-/c/Programs/R/bin/Rscript}

[ -d "$RDIR" ] || { echo "ERROR: r_dir not a directory: $RDIR" >&2; exit 3; }
mkdir -p "$LDIR" 2>/dev/null || { echo "ERROR: cannot create log_dir: $LDIR" >&2; exit 4; }
LDIR=$(cd "$LDIR" && pwd) || exit 4
rm -f "$LDIR"/*.log "$LDIR"/summary.txt

cd "$TDIR" || exit 4
ls test_*.R > /dev/null 2>&1 || { echo "ERROR: no test_*.R in $TDIR" >&2; exit 5; }

echo "r_dir:   $RDIR"
echo "log_dir: $LDIR"
n_bad=0

for f in test_*.R; do
  b=${f%.R}
  echo "running $b ..."
  s0=$(date +%s)
  "$RS" "$f" "$RDIR" > "$LDIR/$b.log" 2>&1
  rc=$?
  s1=$(date +%s)
  [ "$rc" -eq 0 ] || n_bad=$((n_bad + 1))
  np=$(grep -c '^PASS' "$LDIR/$b.log")
  nf=$(grep -c '^FAIL' "$LDIR/$b.log")
  printf '%-28s exit %-3s pass %-4s fail %-4s %ss\n' "$b" "$rc" "$np" "$nf" "$((s1 - s0))" \
    | tee -a "$LDIR/summary.txt"
done

echo "DONE" >> "$LDIR/summary.txt"
exit $([ "$n_bad" -eq 0 ] && echo 0 || echo 1)
