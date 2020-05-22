#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

INPUT="$INPUT"

echo "CREATEDB_PAR: $CREATEDB_QUERY_PAR"
echo "CREATELININDEX_PAR: $CREATELININDEX_PAR"
echo "CONVERT_PAR: $CONVERT_PAR"

echo "COMMAND 1: Create Query Database"
echo "$MMSEQS" createdb "$@" "${TMP_PATH}/query" ${CREATEDB_QUERY_PAR}

if notExists "${TMP_PATH}/query.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "$@" "${TMP_PATH}/query" ${CREATEDB_QUERY_PAR} \
        || fail "query createdb died"
fi

echo "COMMAND 2: Create Target Database"
echo "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR}

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"
fi

echo "COMMAND 3: Create Linear Index"
echo  "$MMSEQS" createlinindex "${TARGET}" "${TMP_PATH}/index_tmp" ${CREATELININDEX_PAR}

if [ -n "${LINSEARCH}" ] && notExists "${TARGET}.linidx"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createlinindex "${TARGET}" "${TMP_PATH}/index_tmp" ${CREATELININDEX_PAR} \
        || fail "createlinindex died"
fi

echo "COMMAND 4: Search"
echo  "$MMSEQS" "${SEARCH_MODULE}" "${TMP_PATH}/query" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR}

INTERMEDIATE="${TMP_PATH}/result"
if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${SEARCH_MODULE}" "${TMP_PATH}/query" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

echo "COMMAND 5: Summarize Result"
echo $RUNNER "$MMSEQS" summarizeresult "${TMP_PATH}/result" "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} 

if [ -n "${GREEDY_BEST_HITS}" ]; then
    if notExists "${TMP_PATH}/result_best.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" summarizeresult "${TMP_PATH}/result" "${TMP_PATH}/result_best" ${SUMMARIZE_PAR} \
            || fail "Search died"
    fi
    INTERMEDIATE="${TMP_PATH}/result_best"
fi

echo "COMMAND 6: Convert Alignments "
echo "$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}${INDEXEXT}" "${INTERMEDIATE}" "${RESULTS}" ${CONVERT_PAR}
echo "PARAMS: ${CONVERT_PAR}"

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${TMP_PATH}/query" "${TARGET}${INDEXEXT}" "${INTERMEDIATE}" "${RESULTS}" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

CLOUD_SEARCH=fb-pruner
echo "COMMAND 6.5: Cloud Search"
echo ${CLOUD_SEARCH} ${QUERY} ${TARGET} -p 4 -i ${RESULTS} 

# ${CLOUD_SEARCH} ${QUERY} ${TARGET} -p 4 -i ${RESULTS} 

echo "COMMAND 7: Remove Directories"
echo "$MMSEQS" rmdb "${TMP_PATH}/result_best" ${VERBOSITY}

# # remove temporary folders
# if [ -n "${REMOVE_TMP}" ]; then
#     if [ -n "${GREEDY_BEST_HITS}" ]; then
#         # shellcheck disable=SC2086
#         "$MMSEQS" rmdb "${TMP_PATH}/result_best" ${VERBOSITY}
#     fi
#     # shellcheck disable=SC2086
#     "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
#     if [ -z "${LEAVE_INPUT}" ]; then
#         if [ -f "${TMP_PATH}/target" ]; then
#             # shellcheck disable=SC2086
#             "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
#             # shellcheck disable=SC2086
#             "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
#         fi
#         # shellcheck disable=SC2086
#         "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
#         # shellcheck disable=SC2086
#         "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
#     fi
#     rm -rf "${TMP_PATH}/search_tmp"
#     rm -f "${TMP_PATH}/easysearch.sh"
# fi
