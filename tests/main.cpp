#include <gtest/gtest.h>
#include "MemoryLayoutPATest.h"
#include "MemoryLayoutAoATest.h"
#include "MemoryLayoutVectorizedPATest.h"
#include "MemoryLayoutVectorizedGroupedPATest.h"
#include "MemoryLayoutVectorizedAoATest.h"
#include "TriangulatorTest.h"
#include "ComprehensiveTriangulatorOrientationTest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
