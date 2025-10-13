#!/usr/bin/env bash
# Library code sourced into bash scripts.

function _log() {
	local program
	local now
	program=${0##*/}
	now=$(date '+%Y-%m-%d %H:%M:%S.%3N')

	< /dev/null : # Ignore stdin
	echo -e "[${now}] (${program})" "$@"
} >&2
