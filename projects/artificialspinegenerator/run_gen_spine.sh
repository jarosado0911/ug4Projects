#!/usr/bin/env bash
set -euo pipefail

# Path to ugshell and the Lua script (override via env if you like)
UGSHELL="${UGSHELL:-ugshell}"
SCRIPT="${SCRIPT:-gen_spine.lua}"

# Prefer local project paths so helperf.lua / ug_util.lua resolve as expected
export ROOT_PATH="${ROOT_PATH:-$PWD}"
export SCRIPT_PATH="${SCRIPT_PATH:-$PWD/scripts:$PWD}"
export APPS_PATH="${APPS_PATH:-$PWD/apps}"

# Help for this wrapper
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: run_gen_spine.sh [params...]

Params are forwarded to gen_spine.lua, e.g.:
  -outName output_dir
  -neckRad 0.08  -neckLen 0.7
  -headRad 0.29  -headLen 0.58
  -appRad 0.025  -appLen 1.175
  -appHeadLen 0.3 -appHeadRad 0.15

Env overrides:
  UGSHELL=/path/to/ugshell
  SCRIPT=gen_spine.lua (or relative/absolute path)
  ROOT_PATH, SCRIPT_PATH, APPS_PATH
EOF
  exit 0
fi

# Build and run: put params after '--' so UGShell doesn't eat them
cmd=( "$UGSHELL" -ex "$SCRIPT"  "$@" )
printf '>> ' && printf '%q ' "${cmd[@]}" && printf '\n'
exec "${cmd[@]}"

