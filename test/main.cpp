#include <gtest/gtest.h>

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  //::testing::GTEST_FLAG(list_tests) = true;
#if 1
  ::testing::GTEST_FLAG(filter) =
    //"UnsignedDistanceTest/1*OverlappingCircles" ":"
    //"UnsignedDistanceTest/1*CircleInsideCircle" ":"
    //"UnsignedDistanceTest/1*Checkerboard" ":"

#if 0
    "TimingTest*" ":"
#endif

#if 0
    "UnsignedDistanceTest/0*HighAccuracy" ":"
    "UnsignedDistanceTest/1*HighAccuracy" ":"
    "UnsignedDistanceTest/2*HighAccuracy" ":"
    "UnsignedDistanceTest/3*HighAccuracy" ":"
    //"SignedDistanceTest/0*HighAccuracy" ":"
    //"SignedDistanceTest/1*HighAccuracy" ":"
    //"SignedDistanceTest/2*HighAccuracy" ":"
#endif

#if 0
    "UnsignedDistanceTest*" ":"
#endif

#if 1
    "SignedDistanceTest*" ":"
#endif


#if 0
    "UniformSpeedEikonalSolverTest*" ":"
    "HighAccuracyUniformSpeedEikonalSolverTest*" ":"
    "VaryingSpeedEikonalSolverTest*" ":"
    "HighAccuracyVaryingSpeedEikonalSolverTest*" ":"
    "DistanceSolverTest*" ":"
#endif

#if 0
    "UnsignedDistanceTest/*EikonalSolverFailThrows" ":"
#endif

#if 0
    "UnsignedDistanceTest/0*NonUniformGridSpacing" ":"
    "UnsignedDistanceTest/1*NonUniformGridSpacing" ":"
    "UnsignedDistanceTest/2*NonUniformGridSpacing" ":"
    "UnsignedDistanceTest/3*NonUniformGridSpacing" ":"
#endif

#if 0
    "BridsonDistanceTest*" ":"
#endif

#if 1
    //"UnsignedDistanceTest/0*VaryingSpeed" ":"
    //"UnsignedDistanceTest/1*VaryingSpeed" ":"
    //"UnsignedDistanceTest/2*VaryingSpeed" ":"
    //"UnsignedDistanceTest/3*VaryingSpeed" ":"
#endif

#if 0
    "UnsignedDistanceTest/1*CircleInsideCircle" ":"
    "SignedDistanceTest/0*CircleInsideCircleThrows" ":"
#endif
      ;
#endif
  return RUN_ALL_TESTS();
}
