#include <gtest/gtest.h>

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  //::testing::GTEST_FLAG(list_tests) = true;
  ::testing::GTEST_FLAG(filter) =
    //"UnsignedDistanceTest/1*OverlappingCircles" ":"
    //"UnsignedDistanceTest/1*CircleInsideCircle" ":"
    //"UnsignedDistanceTest/1*Checkerboard" ":"

#if 0
    "UnsignedDistanceTest/0*HighAccuracy" ":"
    "UnsignedDistanceTest/1*HighAccuracy" ":"
    "UnsignedDistanceTest/2*HighAccuracy" ":"
    "UnsignedDistanceTest/3*HighAccuracy" ":"
    "SignedDistanceTest/0*HighAccuracy" ":"
    "SignedDistanceTest/1*HighAccuracy" ":"
    "SignedDistanceTest/2*HighAccuracy" ":"
#endif

    //"SignedDistanceTest/0*Checkerboard" ":"
    //"SignedDistanceTest/0*OverlappingCircles" ":"

#if 1
    //"UnsignedDistanceTest/0*VaryingSpeed" ":"
    "UnsignedDistanceTest/1*VaryingSpeed" ":"
    //"UnsignedDistanceTest/2*VaryingSpeed" ":"
    //"UnsignedDistanceTest/3*VaryingSpeed" ":"
#endif

#if 0
    "UnsignedDistanceTest/1*CircleInsideCircle" ":"
    "SignedDistanceTest/0*CircleInsideCircleThrows" ":"
#endif
      ;
  return RUN_ALL_TESTS();
}
