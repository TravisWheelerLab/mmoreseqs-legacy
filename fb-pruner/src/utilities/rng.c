/*******************************************************************************
 *  FILE:      rng.c
 *  PURPOSE:   RNG Object: generate random numbers for different data types.
 *
 *  AUTHOR:    Dave Rich
 *  BUG:     
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* local imports */
#include "../objects/structs.h"
#include "utilities.h"
#include "../objects/objects.h"

/* header */
#include "rng.h"

/*
 *  FUNCTION:  RNG_Init()
 *  SYNOPSIS:  Initializes random number generator.
 */
void RNG_Init()
{
   srand(time(NULL));
}

/*
 *  FUNCTION:  RNG_Generate()
 *  SYNOPSIS:  Generate random number.
 */
inline
int RNG_Generate()
{
   return rand();
}

/*
 *  FUNCTION:  RNG_INT()
 *  SYNOPSIS:  Generate random int.
 */
inline
int RNG()
{
   int r = RNG_Generate();
   return r;
}

/*
 *  FUNCTION:  RNG_Range()
 *  SYNOPSIS:  Generate random int in range (beg, end]
 */
inline
int RNG_Range(    int beg,
                  int end )
{
   INT range = end - beg;
   INT r = ( RNG_INT() % range ) + beg;
   return r;
}


/*
 *  FUNCTION:  RNG_INT()
 *  SYNOPSIS:  Generate random int.
 */
inline
INT RNG_INT()
{
   INT r = (INT) RNG_Generate();
   return r;
}

/*
 *  FUNCTION:  RNG_INT_Range()
 *  SYNOPSIS:  Generate random int in range (beg, end]
 */
inline
INT RNG_INT_Range(   INT beg,
                     INT end )
{
   INT range = end - beg;
   INT r = ( RNG_INT() % range ) + beg;
   return r;
}

/*
 *  FUNCTION:  RNG_CHAR()
 *  SYNOPSIS:  Generate random int.
 */
inline
CHAR RNG_CHAR()
{
   CHAR r = (CHAR) RNG_Generate();
   return r;
}

/*FLT
 *  FUNCTION:  RNG_CHAR_Range()
 *  SYNOPSIS:  Generate random int in range (beg, end]
 */
inline
CHAR RNG_CHAR_Range(    CHAR beg,
                        CHAR end )
{
   CHAR range = end - beg;
   CHAR r = ( RNG_CHAR() % range ) + beg;
   return r;
}

/*
 *  FUNCTION:  RNG_FLT()
 *  SYNOPSIS:  Generate random int.
 */
inline
FLT RNG_FLT()
{
   FLT r = (FLT) RNG_Generate();
   return r;
}

/*
 *  FUNCTION:  RNG_FLT_Range()
 *  SYNOPSIS:  Generate random int in range (beg, end]
 */
inline
FLT RNG_FLT_Range(   FLT beg,
                     FLT end )
{
   FLT range = end - beg;
   FLT r = ( RNG_FLT() * range ) / (FLT)RAND_MAX + beg;
   return r;
}

/*
 *  FUNCTION:  RNG_DBL()
 *  SYNOPSIS:  Generate random int.
 */
inline
DBL RNG_DBL()
{
   DBL r = (DBL) RNG_Generate();
   return r;
}

/*
 *  FUNCTION:  RNG_DBL_Range()
 *  SYNOPSIS:  Generate random int in range (beg, end]
 */
inline
DBL RNG_DBL_Range(   DBL beg,
                     DBL end )
{
   DBL range = end - beg;
   DBL r = ( RNG_DBL() * range ) / (FLT)RAND_MAX + beg;
   return r;
}
