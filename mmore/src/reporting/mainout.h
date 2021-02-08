/*******************************************************************************
 *  FILE:      mainout.h
 *  PURPOSE:   Reporting Output for Main Pipeline output.
 *             Styled after HMMER output.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _MAINOUT_H
#define _MAINOUT_H

/* === MAINOUT OUTPUT === */
/* EXAMPLE:
 *
   [header]
   # hmmsearch :: search profile(s) against a sequence database
   # HMMER 3.3 (Nov 2019); http://hmmer.org/
   # Copyright (C) 2019 Howard Hughes Medical Institute.
   # Freely distributed under the BSD open source license.
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # query HMM file:                  test-input/3-PAP.hmm
   # target sequence database:        test-input/3-PAP.fa
   # per-seq hits tabular output:     tblout.csv
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   [entry]
   Query:       3-PAP  [M=133]
   Accession:   PF12578.1
   Description: Myotubularin-associated protein
   Scores for complete sequences (score includes all domains):
      --- full sequence ---   --- best 1 domain ---    -#dom-
       E-value  score  bias    E-value  score  bias    exp  N  Sequence                 Description
       ------- ------ -----    ------- ------ -----   ---- --  --------                 -----------
       3.2e-12   34.5   0.0    1.5e-08   22.6   0.0    2.5  2  3-PAP/16/510-647/718-827  domains: MTMRA_DANRE/548-685 C3Z9W9
       2.7e-11   31.5   0.0    3.9e-10   27.7   0.0    2.4  2  3-PAP/13/1-136/374-513    domains: MTMRB_MOUSE/553-688 B7Q8P0
       1.1e-07   19.9   0.0    2.1e-07   18.9   0.0    1.5  1  3-PAP/14/86-218/365-501   domains: MTMRC_PONAB/559-691 A4HUS9

   Domain annotation for each sequence (and alignments):
   >> 3-PAP/16/510-647/718-827  domains: MTMRA_DANRE/548-685 C3Z9W9_BRAFL/506-615
      #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
    ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
      1 !    9.1   0.0   0.00022   0.00022      68     106 ..     582     618 ..     557     637 .. 0.78
      2 !   22.6   0.0   1.5e-08   1.5e-08      75     110 ..     769     803 ..     719     824 .. 0.73

     Alignments for each domain:
     == domain 1  score: 9.1 bits;  conditional E-value: 0.00022
                        3-PAP  68 LePqcrilDlevWdqCYfRWlPvLeikgGGqpqvDLfnR 106
                                  L Pq     l+vW+  + RW+P  +i +GG   v  f++
     3-PAP/16/510-647/718-827 582 LLPQLLPSHLSVWKLYFLRWVPEAQIPHGGP--VTAFHK 618
                                  667777789********************95..344444 PP

     == domain 2  score: 22.6 bits;  conditional E-value: 1.5e-08
                        3-PAP  75 lDlevWdqCYfRWlPvLeikgGGqpqvDLfnRllLs 110
                                    l++W qCY RW+P     gGG p  + f+  lL 
     3-PAP/16/510-647/718-827 769 AGLKLWTQCYMRWIPWAHTVGGGPPS-EYFHQCLLV 803
                                  5799************9998888665.777777664 PP

   ...
   [footer]
   Internal pipeline statistics summary:
   -------------------------------------
   Query model(s):                            1  (133 nodes)
   Target sequences:                          3  (2475 residues searched)
   Passed MSV filter:                         3  (1); expected 0.1 (0.02)
   Passed bias filter:                        3  (1); expected 0.1 (0.02)
   Passed Vit filter:                         3  (1); expected 0.0 (0.001)
   Passed Fwd filter:                         3  (1); expected 0.0 (1e-05)
   Initial search space (Z):                  3  [actual number of targets]
   Domain search space  (domZ):               3  [number of targets reported over threshold]
   # CPU time: 0.01u 0.00s 00:00:00.01 Elapsed: 00:00:00.01
   # Mc/sec: 24.80
   //
   [ok]
 *
 */

/*! FUNCTION:   REPORT_stdout_header()
 *  SYNOPSIS:   Print Header to Output <fp>.
 *              (modeled after HMMER output, see example)
 */
void 
REPORT_stdout_header(   WORKER*     worker,
                        FILE*       fp );

/*!   FUNCTION:   REPORT_stdout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void 
REPORT_stdout_entry(    WORKER*  worker,
                        RESULT*  result,
                        FILE*    fp );

/*!   FUNCTION:   REPORT_stdout_footer()
 *    SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void 
REPORT_stdout_footer(   WORKER*  worker,
                        FILE*    fp );

/*!   FUNCTION:   REPORT_stdout_footer_search_summary()
 *    SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void 
REPORT_stdout_footer_search_summary(   WORKER*  worker,
                                       FILE*    fp );

/*!   FUNCTION:   REPORT_stdout_footer_time_summary()
 *    SYNOPSIS:   Print Runtime Summary.
 */
void 
REPORT_stdout_footer_time_summary(     WORKER*  worker,
                                       FILE*    fp );



#endif /* _REPORT_H */