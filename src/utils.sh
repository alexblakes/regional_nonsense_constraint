# Library code sourced into bash scripts.

MyLog() {
	# Timestamp a message. Send it to $LOGFILE and stderr.
	echo -e $(date +"%y-%m-%dT%T") "${@}" "\n" | tee -a $FILE_LOG 1>&2
}