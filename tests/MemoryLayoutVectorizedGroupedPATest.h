#pragma once

#include <gtest/gtest.h>
#include "Datastructures.h"

TEST(MemoryLayoutVectorizedGroupedPA, ReadWriteLength1Dimension1) {
    MemoryLayoutVectorizedGroupedPA<int, 1> layout;
    layout.ensure(0);
    layout(0, 0) = 11;
    EXPECT_EQ(layout(0, 0), 11);
}

TEST(MemoryLayoutVectorizedGroupedPA, ReadWriteLength2Dimension1) {
    MemoryLayoutVectorizedGroupedPA<int, 1> layout;
    layout.ensure(1);
    layout(0, 0) = 11;
    layout(1, 0) = 21;
    EXPECT_EQ(layout(0, 0), 11);
    EXPECT_EQ(layout(1, 0), 21);
}