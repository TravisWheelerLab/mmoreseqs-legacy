#include <stdio.h>
#include "funcs.h"

#define xstr(s) str(s)
#define str(s) #s

/* === MAIN ENTRY-POINT TO PROGRAM === */
int 
main ( int argc, char *argv[] )
{
   printf("Hello Makefile\n");

   #ifndef BUILD
   #define BUILD test
   #endif
   #define BUILD_STR xstr(BUILD)

   printf("BUILD DEFINED: %s\n", BUILD_STR );


   float ans;

   ans = MATH_Add( 2.0, 3.0 );

   printf("ANSWER: %f\n", ans);

   return 0;
}