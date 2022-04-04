/*******************************************************************************
 *  - FILE:      m8out.h
 *  - DESC:    Methods for generating BLAST-style .m8 output.
 *******************************************************************************/

#ifndef _M8OUT_H
#define _M8OUT_H

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
   #  qid         tid                                  %id    aln   mis  gap  qb
 qe   tb    te    eval       bit #  ----------  -------------------------------
 ---    ---   ---  ---  ---  ---  ---   ---   ----       --- [entry/]
      Aldolase_II  Aldolase_II/544/29-203/211-382      0.252  151   105  0    21
 171  51    191   2.631E-07  54 Aldolase_II  decoy94676 0.230  91    69   0 66
 156  164   253   2.572E+01  29 Aldolase_II  decoy95461 0.553  24    10   0 16
 39   102   124   4.583E+01  29 [footer/] # [ok]
 *
 */

/*!   FUNCTION:   REPORT_m8out_header()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_m8out_header(WORKER* worker, FILE* fp);

/*!   FUNCTION:   REPORT_m8out_entry()
 *    SYNOPSIS:   Print all alignment data for current search
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_m8out_entry(WORKER* worker, RESULT* result, FILE* fp);

/*!   FUNCTION:   REPORT_m8out_footer()
 *    SYNOPSIS:   Print footer
 *                (modeled after HMMER, see example)
 */
STATUS_FLAG
REPORT_m8out_footer(WORKER* worker, FILE* fp);

#endif /* _M8OUT_H */
