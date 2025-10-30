#!/usr/bin/env bash
# Library code sourced into bash scripts.

# Logging function
function _log() {
	local program
	program=${0##*/}
	local now
	now=$(date '+%Y-%m-%d %H:%M:%S.%3N')

	cat > /dev/null # Absorb stdin. Prevents issues with back-to-back `tee`s.

	echo -e "[${now}] (${program})" "$@"

} >&2