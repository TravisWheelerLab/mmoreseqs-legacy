/*******************************************************************************
 *  - FILE:  domain_def.c
 *  - DESC:  DOMAIN_DEF object.
 *******************************************************************************/

#ifndef _DOMAIN_DEF_H
#define _DOMAIN_DEF_H

/** FUNCTION:  DOMAIN_DEF_Create()
 *  SYNOPSIS:  Create new DOMAIN_DEF object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
DOMAIN_DEF* DOMAIN_DEF_Create();

/** FUNCTION:  DOMAIN_DEF_Destroy()
 *  SYNOPSIS:  Frees DOMAIN_DEF object and returns pointer.
 */
DOMAIN_DEF* DOMAIN_DEF_Destroy(DOMAIN_DEF* dom_def);

/** FUNCTION:  DOMAIN_DEF_Reuse()
 *  SYNOPSIS:
 */
int DOMAIN_DEF_Reuse(DOMAIN_DEF* dom_def);

/** FUNCTION:  DOMAIN_DEF_Resize()
 *  SYNOPSIS:
 */
int DOMAIN_DEF_Resize(DOMAIN_DEF* dom_def, int size);

/** FUNCTION:  DOMAIN_DEF_Create()
 *  SYNOPSIS:  Create new DOMAIN_DEF object and returns pointer.
 *             Most data is left NULL to be supplied by WORK_init().
 */
int DOMAIN_DEF_GrowTo(DOMAIN_DEF* dom_def, int size);

#endif /* _DOMAIN_DEF_H */
