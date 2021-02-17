#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore-get-tools.sh	
#	AUTHOR:		David Rich
#	DESC: 		Finds associated tools (HMMER, MMSEQS, MMORE) for pipeline.
#              Should be source'd into main script.
###########################################################################

# only echoes if has higher verbosity level
function echo_v
{
	if (( $VERBOSE >= $1 ))
	then
		echo $2
	fi
}

echo_v 1 "In 'mmore-get-tools.sh'..."

# default whether to use system or local tools
USE_LOCAL_TOOLS="${USE_LOCAL_TOOLS:-0}"

# if system tools allowed, then try to using system installs
{
	echo_v 3 "SEARCHING FOR SYSTEM TOOLS..."
	# system installs
	MMSEQS_SYSTEM=$(which mmseqs)
	HMMER_SYSTEM=$(which hmmbuild)
	MMORE_SYSTEM=$(which mmore)
}

# check if program is installed on system
{
	# which mmseqs
	if [ -z "$MMSEQS_SYSTEM" ]
	then
		echo "MMSEQS is not installed on system."
	fi
	# which hmmbuild
	if [ -z "$HMMER_SYSTEM" ]
	then
		echo "HMMER is not installed on system."
	fi
	# which mmore
	if [ -z "$MMORE_SYSTEM" ]
	then
		echo "MMORE is not installed on system."
	fi
}

# if not installed on system, add local install to path (paths passed as environmental vars in main)
{
	echo_v 3 "SEARCHING FOR LOCAL TOOL..."
	echo_v 3 "MMORE_DIR: 	${MMORE_DIR}"
	echo_v 3 "MMSEQS_DIR: 	${MMSEQS_DIR}"
	echo_v 3 "HMMER_DIR: 	${HMMER_DIR}"

	# add local tools to path
	export PATH="${MMSEQS_DIR}:$PATH"
	export PATH="${HMMER_DIR}:$PATH"
	export PATH="${MMORE_DIR}:$PATH"
	# local installs
	MMSEQS_LOCAL=$(which mmseqs)
	HMMER_LOCAL=$(which hmmbuild)
	MMORE_LOCAL=$(which mmore)
}

# look for other builds of mmore
{
	MMORE_PROGS=$(ls ${MMORE_DIR}/mmore-*)
	echo PROGS: $(ls ${MMORE_DIR}/mmore-*)
	for PROG in $MMORE_PROGS
	do 
		echo $PROG
		ALT_MMORE=$(which $PROG)
		echo LOCAL: $MMORE_LOCAL
		if [ -z "$MMORE_LOCAL" ]
		then
			echo_v 1 "USING NON-STANDARD BUILD OF MMORE: $(GET_FILE_NAME $ALT_MMORE)"
			MMORE_LOCAL="$ALT_MMORE"
		fi
	done
}

# check if program is installed locally
{
	# which mmseqs
	if [ -z "$MMSEQS_LOCAL" ]
	then
		echo "MMSEQS is not installed locally."
	fi
	# which hmmbuild
	if [ -z "$HMMER_LOCAL" ]
	then
		echo "HMMER is not installed locally."
	fi
	# which mmore
	if [ -z "$MMORE_LOCAL" ]
	then
		echo "MMORE is not installed locally."
	fi
}

# set default tool (which tool to prioritize)
{
	if (( $USE_LOCAL_TOOLS != 1 )); then 
		# use local tools 
		MMSEQS="${MMSEQS:-$MMSEQS_SYSTEM}"
		HMMER="${HMMER:-$HMMER_SYSTEM}"
		MMORE="${MMORE:-$MMORE_SYSTEM}"
		# use system tools 
		MMSEQS="${MMSEQS:-$MMSEQS_LOCAL}"
		HMMER="${HMMER:-$HMMER_LOCAL}"
		MMORE="${MMORE:-$MMORE_LOCAL}"
	else 
		# use system tools 
		MMSEQS="${MMSEQS:-$MMSEQS_LOCAL}"
		HMMER="${HMMER:-$HMMER_LOCAL}"
		MMORE="${MMORE:-$MMORE_LOCAL}"
		# use local tools 
		MMSEQS="${MMSEQS:-$MMSEQS_SYSTEM}"
		HMMER="${HMMER:-$HMMER_SYSTEM}"
		MMORE="${MMORE:-$MMORE_SYSTEM}"
	fi 
}

# if program still can't be found, throw error.
{
	# which mmseqs
	if [ -z "$MMSEQS" ]
	then
		echo "ERROR: MMSEQS is not installed on system or locally."
		ERROR="MMSEQS_NOT_INSTALLED"
	fi
	# which hmmbuild
	if [ -z "$HMMER" ]
	then
		echo "ERROR: HMMER is not installed on system or locally."
		ERROR="HMMER_NOT_INSTALLED"
	fi
	# which mmore
	if [ -z "$MMORE" ]
	then
		echo "ERROR: MMORE is not installed on system or locally."
		ERROR="MMORE_NOT_INSTALLED"
	fi

	if [ ! -z "$MMSEQS" ] && [ ! -z "$HMMER" ] && [ ! -z "$MMORE" ]
	then 
		echo_v 3 "All tools successfully found."
	fi
}

echo_v 3 "'mmore_get-tools.sh' Loaded."
