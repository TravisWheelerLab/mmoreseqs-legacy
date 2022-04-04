/*******************************************************************************
 *  - FILE:      work_loader.c
 *  - DESC:    Pipelines Workflow Subroutines.
 *             WORK interfaces between pipeline WORKER object and various
 *functions. Loads query sequences and target hmm profiles.
 *******************************************************************************/

#ifndef _WORK_LOADER
#define _WORK_LOADER

/*! FUNCTION:  	WORK_load_mmseqs_input()
 *  SYNOPSIS:  	Loads mmseqs input m8 file into <results_in>, located at
 * <mmseqs_m8_filein>. Verifies that search index range does not exceed bounds
 * of list.
 */
void WORK_load_mmseqs_file(WORKER* worker);

/*! FUNCTION:  	WORK_load_mmseqs_by_id()
 *  SYNOPSIS:  	Load <i>th mmseqs input from .m8 <mmseqs_data> list into
 * <worker>.
 */
void WORK_load_mmseqs_by_id(WORKER* worker, int i);

/*! FUNCTION:  	WORK_load_mmseqs_alignment()
 *  SYNOPSIS:  	Load mmseq's viterbi alignment <trace_vit> into <worker>.
 *                Query and Target should already be loaded.
 */
void WORK_load_mmseqs_alignment(WORKER* worker);

/*! FUNCTION:  	WORK_load_target()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <mmseqs_data>'s target name field.
 */
void WORK_load_target(WORKER* worker);

/*! FUNCTION:  	WORK_load_query()
 *  SYNOPSIS:  	Loads <query> HMM_PROFILE by <mmseqs_data>'s query name field.
 */
void WORK_load_query(WORKER* worker);

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_target_by_id(WORKER* worker, int id);

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_query_by_id(WORKER* worker, int id);

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_target_by_name(WORKER* worker, char* name);

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <name> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_query_by_name(WORKER* worker, char* name);

/*! FUNCTION:  	WORK_load_target_by_id()
 *  SYNOPSIS:  	Loads <target> HMM_PROFILE by <t_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_target_by_findex_id(WORKER* worker, int index_id);

/*! FUNCTION:  	WORK_load_query_by_id()
 *  SYNOPSIS:  	Loads <query> SEQUENCE by <q_index> F_INDEX <id> field.
 *                Depends on <task> settings in <worker>.
 */
void WORK_load_query_by_findex_id(WORKER* worker, int index_id);

#endif /* _WORK_LOADER */
