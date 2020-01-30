#!/bin/sh -e
# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";


INPUT="$1"
TARGET="$2"
TMP_PATH="$4"

echo "INPUT: $INPUT, TARGET: $TARGET, TMP: $TMP_PATH, RUNNER: $RUNNER"
echo "THIS IS A TEST FROM data/blastp.sh !"

STEP=0
STEPS=${STEPS:-1}
ALN_RES_MERGE="$TMP_PATH/aln_0"
while [ "$STEP" -lt "$STEPS" ]; do
    SENS_PARAM=SENSE_${STEP}
    eval SENS="\$$SENS_PARAM"

    # call prefilter module
    if notExists "$TMP_PATH/pref_$STEP.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "$INPUT" "$TARGET" "$TMP_PATH/pref_$STEP" $PREFILTER_PAR -s "$SENS" \
            || fail "Prefilter died"
    fi

    # call alignment module
    if [ "$STEPS" -eq 1 ]; then
        if notExists "$3.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$TARGET${ALIGNMENT_DB_EXT}" "$TMP_PATH/pref_$STEP" "$3" $ALIGNMENT_PAR  \
                || fail "Alignment died"
        fi
    else
        if notExists "$TMP_PATH/aln_$STEP.dbtype"; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "$INPUT" "$TARGET${ALIGNMENT_DB_EXT}" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" $ALIGNMENT_PAR  \
                || fail "Alignment died"
        fi
    fi

    # INSERT CLOUD SEARCH HERE!!
    echo "THIS IS WHERE CLOUD SEARCH GOES..."
    CLOUD_SEARCH_RUNNER="/Users/Devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/build/cloud_fwdbck.exe"

    # only merge results after first step
    if [ "$STEP" -gt 0 ]; then
        if notExists "$TMP_PATH/aln_${SENS}.hasmerged"; then
            if [ "$STEP" -lt $((STEPS-1)) ]; then
                "$MMSEQS" mergedbs "$1" "$TMP_PATH/aln_merge" "$ALN_RES_MERGE" "$TMP_PATH/aln_$STEP" \
                || fail "Mergedbs died"
                ALN_RES_MERGE="$TMP_PATH/aln_merge"
            else
                "$MMSEQS" mergedbs "$1" "$3" "$ALN_RES_MERGE" "$TMP_PATH/aln_$STEP" \
                || fail "Mergedbs died"
            fi
            touch "$TMP_PATH/aln_${STEP}.hasmerged"
        fi
    fi

    NEXTINPUT="$TMP_PATH/input_step$STEP"
    #do not create subdb at last step
    if [ "$STEP" -lt "$((STEPS-1))" ]; then
        if notExists "$TMP_PATH/order_step$STEP.dbtype"; then
            awk '$3 < 2 { print $1 }' "$TMP_PATH/aln_$STEP.index" > "$TMP_PATH/order_step$STEP" \
                || fail "Awk step $STEP died"
        fi

        if [ ! -s "$TMP_PATH/order_step$STEP" ]; then break; fi

        if notExists "$NEXTINPUT.dbtype"; then
            "$MMSEQS" createsubdb "$TMP_PATH/order_step$STEP" "$INPUT" "$NEXTINPUT" \
                || fail "Order step $STEP died"
        fi
    fi
    INPUT="$NEXTINPUT"
    STEP="$((STEP+1))"
done


if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    STEP=0
    while [ "$STEP" -lt "$STEPS" ]; do
        "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP"
        "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP"
        "$MMSEQS" rmdb "${TMP_PATH}/input_step$STEP"
        #NEXTINPUT="$TMP_PATH/input_step$STEP" # this line is unused
        STEP="$((STEP+1))"
    done
    rm -f "$TMP_PATH/blastp.sh"
fi
