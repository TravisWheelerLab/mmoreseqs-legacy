#!/usr/bin/bash
###########################################################################
#	NAME: 		mmore-get-tools.sh	
#	AUTHOR:		David Rich
#	DESC: 		Finds associated tools (HMMER, MMSEQS, MMORE) for pipeline.
#              Should be 'source'd into main script.
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

# if system tools allowed, then try to using system installs
USE_LOCAL_TOOLS="${USE_LOCAL_TOOLS:-0}"
if (( USE_LOCAL_TOOLS == 0 ))
then
	echo_v 3 "SEARCHING FOR SYSTEM TOOLS..."
	# which mmseqs
	MMSEQS_DEFAULT=$(which mmseqs)
	# which hmmbuild
	HMMER_DEFAULT=$(which hmmbuild)
	# which mmore
	MMORE_DEFAULT=$(which mmore)
fi

# if not install on system, add local install to path (paths passed as environmental vars in main)
if [ -z "$MMSEQS_DEFAULT" ]
then
	echo_v 3 "USING LOCAL MMSEQS: ${MMSEQS_DIR}"
	export PATH="${MMSEQS_DIR}:$PATH"
fi
if [ -z "$HMMER_DEFAULT" ]
then
	echo_v 3 "USING LOCAL HMMER: ${HMMER_DIR}"
	export PATH="${HMMER_DIR}:$PATH"
fi
if [ -z "$MMORE_DEFAULT" ]
then
	echo_v 3 "USING LOCAL MMORE: ${MMORE_DIR}"
	export PATH="${MMORE_DIR}:$PATH"
fi

# if program still can't be found, throw error.
# which mmseqs
MMSEQS=$(which mmseqs)
if [ -z "$MMSEQS" ]
then
	echo "ERROR: hmmer is not installed on system or locally."
	exit 1
fi
# which hmmbuild
HMMBUILD=$(which hmmbuild)
if [ -z "$HMMBUILD" ]
then
	echo "ERROR: hmmer is not installed on system or locally."
	exit 1
fi
# which mmore
MMORE=$(which mmore)
if [ -z "$MMORE" ]
then
	echo "ERROR: mmore is not installed on system or locally."
	exit 1
fi
