#pragma once

#include <gtest/gtest.h>
#include "Datastructures.h"

TEST(MemoryLayoutAoA, ReadWriteLength1Dimension1) {
    MemoryLayoutAoA<int, 1> layout;
    layout.ensure(0);
    layout(0, 0) = 11;
    EXPECT_EQ(layout(0, 0), 11);
}

TEST(MemoryLayoutAoA, ReadWriteLength2Dimension1) {
    MemoryLayoutAoA<int, 1> layout;
    layout.ensure(1);
    layout(0, 0) = 11;
    layout(1, 0) = 21;
    EXPECT_EQ(layout(0, 0), 11);
    EXPECT_EQ(layout(1, 0), 21);
}