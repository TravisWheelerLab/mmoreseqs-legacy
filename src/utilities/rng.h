/*******************************************************************************
 *  - FILE:      rng.h
 *  - DESC:    Random Number Generation functions
 *******************************************************************************/

#ifndef _RNG_H
#define _RNG_H

/*! FUNCTION:  RNG_Init()
 *  SYNOPSIS:  Initializes random number generator.
 */
void RNG_Init();

/*! FUNCTION:  RNG_Generate()
 *  SYNOPSIS:  Generate random number.
 */
int RNG_Generate();

/*! FUNCTION:  RNG()
 *  SYNOPSIS:  Generate random int.
 */
int RNG();

/*! FUNCTION:  RNG_Range()
 *  SYNOPSIS:  Generate random int in range (beg, end]
 */
int RNG_Range(int beg, int end);

/*! FUNCTION:  RNG_INT()
 *  SYNOPSIS:  Generate random int.
 */
INT RNG_INT();

/*! FUNCTION:  RNG_INT_Range()CHAR
 *  SYNOPSIS:  Generate random int from range (beg, end]
 */
INT RNG_INT_Range(INT beg, INT end);

#endif /* _RNG_H */
