// Copyright 2017 Tommy Hinks
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

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

#if 1
    "UnsignedArrivalTimeTest/1*" ":"
#endif

#if 0
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
    "*OverlappingCircles" ":"
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

#if 0
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
