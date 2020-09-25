/*******************************************************************************
 *  FILE:      report.h
 *  PURPOSE:   Reporting Subroutines for generating output.
 *
 *  AUTHOR:    Dave Rich
 *******************************************************************************/

#ifndef _REPORT_H
#define _REPORT_H


/* === STDOUT OUTPUT === */
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

/* FUNCTION:   REPORT_stdout_header()
 * SYNOPSIS:   Print Header to Output <fp>.
 *             (modeled after HMMER, see example)
 */
void REPORT_stdout_header(    WORKER*     worker,
                              FILE*       fp );

/*    FUNCTION:   REPORT_stdout_footer()
 *    SYNOPSIS:   Print Summary Statistics after all searches completed.
 *                (modeled after HMMER, see example)
 */
void REPORT_stdout_footer(    WORKER*  worker,
                              FILE*    fp );

/*    FUNCTION:   REPORT_stdout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_stdout_entry(  WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp );

/* === TBLOUT OUTPUT === */
/* EXAMPLE:
 *
      [header]
      #                                                                   --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
      # target name            accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
      #    ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
      [entry]
      3-PAP/16/510-647/718-827 -          3-PAP                PF12578.1    3.2e-12   34.5   0.0   1.5e-08   22.6   0.0   2.5   2   0   0   2   2   2   2 domains: MTMRA_DANRE/548-685 C3Z9W9_BRAFL/506-615
		3-PAP/13/1-136/374-513   -          3-PAP                PF12578.1    2.7e-11   31.5   0.0   3.9e-10   27.7   0.0   2.4   2   0   0   2   2   2   1 domains: MTMRB_MOUSE/553-688 B7Q8P0_IXOSC/504-643
		3-PAP/14/86-218/365-501  -          3-PAP                PF12578.1    1.1e-07   19.9   0.0   2.1e-07   18.9   0.0   1.5   1   0   0   1   1   1   1 domains: MTMRC_PONAB/559-691 A4HUS9_LEIIN/8-144
      [footer]
      #
      # Program:         hmmsearch
      # Version:         3.3 (Nov 2019)
      # Pipeline mode:   SEARCH
      # Query file:      test-input/3-PAP.hmm
      # Target file:     test-input/3-PAP.fa
      # Option settings: hmmsearch --tblout tblout.csv test-input/3-PAP.hmm test-input/3-PAP.fa 
      # Current dir:     /home/devreckas/Google-Drive/Wheeler-Labs/Personal_Work/fb-pruner/fb-pruner
      # Date:            Fri Aug 28 00:42:33 2020
      # [ok]
 *
 */


/*    FUNCTION:   REPORT_tblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_tblout_header( WORKER*  worker,
                           FILE*    fp );

/*    FUNCTION:   REPORT_alignment()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_tblout_entry(  WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp );

/*    FUNCTION:   REPORT_tblout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)              
 */
void REPORT_tblout_footer(    WORKER*  worker,
                              FILE*    fp );

/* === M8 OUTPUT === */
/* Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1,2) identifiers query and target sequences/profiles, 
   - (3) sequence identity, 
   - (4) alignment leng(5) number of mismatches, 
   - (6) number of gap openings, 
   - (7,8) domain (start,end)-position in query and in target, 
   - (9,10) domain
   - (11) E-value, and 
   - (12) bit score.
 */

/* EXAMPLE:
 *
   [header/]
   #  qid         tid                                  %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ----------  -------------------------------      ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*    FUNCTION:   REPORT_m8out_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_m8out_header(  WORKER*  worker,
                           FILE*    fp );

/*    FUNCTION:   REPORT_m8out_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_m8out_entry(   WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp );

/*    FUNCTION:   REPORT_m8out_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_m8out_footer(    WORKER*  worker,
                              FILE*    fp );

/* === MYOUT FUNCTIONS === */

/* === M OUTPUT === */
/* BLAST-style output 
   Description: The file is formatted as a tab-separated list with 12 columns: 
   - (1,2) identifiers query and target sequences/profiles, 
   - (3) sequence identity, 
   - (4) alignment leng(5) number of mismatches, 
   - (6) number of gap openings, 
   - (7,8) domain (start,end)-position in query and in target, 
   - (9,10) domain
   - (11) E-value, and 
   - (12) bit score.
 */

/* EXAMPLE:
 *
   [header/]
   #  qid         tid                                  %id    aln   mis  gap  qb   qe   tb    te    eval       bit 
   #  ----------  -------------------------------      ---    ---   ---  ---  ---  ---  ---   ---   ----       ---
   [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21   171  51    191   2.631E-07  54
      Aldolase_II  decoy94676                          0.230  91    69   0    66   156  164   253   2.572E+01  29
      Aldolase_II  decoy95461                          0.553  24    10   0    16   39   102   124   4.583E+01  29
   [footer/]
   # [ok]
 *
 */

/*    FUNCTION:   REPORT_myout_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_myout_header(  WORKER*  worker,
                           FILE*    fp );

/*    FUNCTION:   REPORT_myout_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
void REPORT_myout_entry(   WORKER*  worker,
                           RESULT*  result,
                           FILE*    fp );

/*    FUNCTION:   REPORT_myout_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)              
 */
void REPORT_myout_footer(    WORKER*  worker,
                              FILE*    fp );

/* === UTILTITY FUNCTIONS === */

/*    FUNCTION:   REPORT_horizontal_rule()
 *    SYNOPSIS:   Print a horizontal rule. 
 */
void REPORT_horizontal_rule( FILE* fp );

/*    FUNCTION:   REPORT_horizontal_rule()
 *    SYNOPSIS:   Print a horizontal rule of specified length.
 */
void REPORT_hr_size( FILE* fp, 
                     int   length );


#endif /* _REPORT_H */