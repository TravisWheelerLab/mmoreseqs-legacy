/*******************************************************************************
 *  @file parser.h
 *  @brief Parses .hmm and .fasta files
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _HMM_PARSER_H
#define _HMM_PARSER_H

void hmmprofile_Create(HMM_PROFILE *prof, char *_filename_);
void hmmprofile_Config(HMM_PROFILE *prof, int mode);
void hmmprofile_CalcOccupancy(HMM_PROFILE *prof, float *occ);
void hmmprofile_ReconfigLength(HMM_PROFILE *prof, int L);
void hmmprofile_Display(HMM_PROFILE *prof);
void hmmprofile_Save(HMM_PROFILE *prof, char *_filename_);

void seq_Create(SEQ *seq, char *_filename_);
void seq_Display(SEQ *seq);

void submat_Create(SUBMAT *submat, char *_filename_);
int submat_Keymap(char a, char b);
float submat_Get(SUBMAT *submat, char q, char t);
void submat_Display(SUBMAT *submat);

void results_Display(RESULTS *results);

float negln2real(float negln_prob);
float real2negln(float real_prob);

#endif /* _HMM_PARSER_H */
