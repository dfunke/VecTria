#include <gtest/gtest.h>
#include "MemoryLayoutPATest.h"
#include "MemoryLayoutAoATest.h"
#include "MemoryLayoutVectorizedPATest.h"
#include "MemoryLayoutVectorizedGroupedPATest.h"
#include "MemoryLayoutVectorizedAoATest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
