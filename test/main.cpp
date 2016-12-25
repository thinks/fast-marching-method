#include <gtest/gtest.h>

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  //::testing::GTEST_FLAG(list_tests) = true;
  //::testing::GTEST_FLAG(filter) = "SignedDistanceTest/3*DifferentUniformSpeed";
  ::testing::GTEST_FLAG(filter) = "*/0.BoxBoundary*";
  return RUN_ALL_TESTS();
}
