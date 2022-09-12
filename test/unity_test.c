/*******************************************************************************
 *  - FILE:  unity_test.c
 *  - DESC:  Entry Point to Unity unit testing.
 *******************************************************************************/

/* import stdlib */
#include <stdio.h>
#include <stdlib.h>

/* import local libs */
#include "unity.h"
// #include "easel.h"

/* import local files */
#include "../src/macros/_macros.h"


void setUp(void)
{
  // set stuff up here
}

void tearDown(void)
{
  // clean stuff up here
}

void test_SimpleTestPasses(void)
{
  int a = 5;
  TEST_ASSERT_EQUAL_INT(a, 5);
}

void test_SimpleTestFails(void)
{
  int b = 10;
  TEST_ASSERT_EQUAL_INT(b, 10);
  TEST_ASSERT_EQUAL_INT(b, 5);
  TEST_ASSERT_EQUAL_INT(b, 2);
}

void test_SimpleTestSegfault(void)
{
  int c[5];
  int d[5] = {0, 1, 2, 3, 3};
  for (int i = 0; i < 5; i++) {
    c[i] = i;
  }
  for (int i = 0; i < 5; i++) {
    printf("[%d]: %d vs %d\n", i, c[i], d[i]);
  }
  TEST_ASSERT_INT_ARRAY_WITHIN(1, d, c, 5);
}

int main(void) 
{
  printf(BUILD_COPYRIGHT);
  UNITY_BEGIN();

  RUN_TEST(test_SimpleTestSegfault);
  RUN_TEST(test_SimpleTestPasses);
  RUN_TEST(test_SimpleTestFails);

  return UNITY_END();
}
