#include <gtest/gtest.h>

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  //::testing::GTEST_FLAG(list_tests) = true;
  ::testing::GTEST_FLAG(filter) =
    //"UnsignedDistanceTest/1*OverlappingCircles" ":"
    //"UnsignedDistanceTest/1*CircleInsideCircle" ":"
    //"UnsignedDistanceTest/1*Checkerboard" ":"
    "UnsignedDistanceTest/1*HighAccuracy" ":"
    "SignedDistanceTest/0*HighAccuracy" ":"
    //"SignedDistanceTest/0*Checkerboard" ":"
    //"SignedDistanceTest/0*OverlappingCircles" ":"
    //"SignedDistanceTest/0*CircleInsideCircleThrows" ":"
      ;
  return RUN_ALL_TESTS();
}
