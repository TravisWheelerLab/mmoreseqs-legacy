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

echo_v 3 "In 'mmore-get-tools.sh'..."

# default whether to use system or local tools
USE_LOCAL_TOOLS="${USE_LOCAL_TOOLS:-0}"

# if system tools allowed, then try to using system installs
{
	echo_v 3 "SEARCHING FOR SYSTEM TOOLS..."
	# system installs
	MMSEQS_SYSTEM=$(which mmseqs)
	HMMER_SYSTEM=$(which hmmbuild)
	HMMBUILD_SYSTEM=$(which hmmbuild)
	HMMEMIT_SYSTEM=$(which hmmemit)
	HMMSEARCH_SYSTEM=$(which hmmsearch)
	MMORE_SYSTEM=$(which mmore)
}

# check if program is installed on system
{
	# which mmseqs
	if [ -z "$MMSEQS_SYSTEM" ]
	then
		echo_v 2 "Warning: MMSEQS is not installed on system."
	fi
	# which hmmbuild
	if [ -z "$HMMER_SYSTEM" ]
	then
		echo_v 2 "Warning: HMMER is not installed on system."
	fi
	# which mmore
	if [ -z "$MMORE_SYSTEM" ]
	then
		echo_v 2 "Warning: MMORE is not installed on system."
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
	HMMBUILD_LOCAL=$(which hmmbuild)
	HMMEMIT_LOCAL=$(which hmmemit)
	HMMSEARCH_LOCAL=$(which hmmsearch)
	MMORE_LOCAL=$(which mmore)
}

# look for other builds of mmore
{
	MMORE_PROGS=$(ls ${MMORE_DIR}/mmore-*)
	echo_v 3 "PROGS: $(ls ${MMORE_DIR}/mmore-*)"
	for PROG in $MMORE_PROGS
	do 
		echo_v 3 "$PROG"
		ALT_MMORE=$(which $PROG)
		echo_v 3 "LOCAL: $MMORE_LOCAL"
		if [ -z "$MMORE_LOCAL" ]
		then
			echo_v 2 "Warning: Using non-standard version of MMORE => $(GET_FILE_NAME $ALT_MMORE)"
			MMORE_LOCAL="$ALT_MMORE"
		fi
	done
}

# check if program is installed locally
{
	# which mmseqs
	if [ -z "$MMSEQS_LOCAL" ]
	then
		echo_v 3 "MMSEQS is not installed locally."
	fi
	# which hmmbuild
	if [ -z "$HMMER_LOCAL" ]
	then
		echo_v 3 "HMMER is not installed locally."
	fi
	# which mmore
	if [ -z "$MMORE_LOCAL" ]
	then
		echo_v 3 "MMORE is not installed locally."
	fi
}

# set default tool (which tool to prioritize)
{
	if (( $USE_LOCAL_TOOLS != 1 )); then 
		# use local tools 
		MMSEQS="${MMSEQS_LOCAL:-$MMSEQS_SYSTEM}"
		HMMER="${HMMER_LOCAL:-$HMMER_SYSTEM}"
		HMMBUILD="${HMMBUILD_LOCAL:-$HMMBUILD_SYSTEM}"
		HMMEMIT="${HMMEMIT_LOCAL:-$HMMEMIT_SYSTEM}"
		HMMSEARCH="${HMMSEARCH_LOCAL:-$HMMSEARCH_SYSTEM}"
		MMORE="${MMORE_LOCAL:-$MMORE_SYSTEM}"
	else 
		# use local tools 
		MMSEQS="${MMSEQS_SYSTEM:-$MMSEQS_LOCAL}"
		HMMER="${HMMER_SYSTEM:-$HMMER_LOCAL}"
		HMMBUILD="${HMMBUILD_SYSTEM:-$HMMBUILD_LOCAL}"
		HMMEMIT="${HMMEMIT_SYSTEM:-$HMMEMIT_LOCAL}"
		HMMSEARCH="${HMMSEARCH_SYSTEM:-$HMMSEARCH_LOCAL}"
		MMORE="${MMORE_SYSTEM:-$MMORE_LOCAL}"
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
echo_v 3 "# HMMBUILD: $HMMBUILD"
echo_v 3 "# HMMEMIT: $HMMEMIT"
echo_v 3 "# MMSEQS: $MMSEQS"
echo_v 3 "# MMORE: $MMORE"