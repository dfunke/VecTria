#pragma once

#include <gtest/gtest.h>
#include "Datastructures.h"

TEST(Geometry3D, inSphere) {
    MemoryLayoutPA<int, 1> layout;
    layout.ensure(0);
    layout(0, 0) = 11;
    EXPECT_EQ(layout(0, 0), 11);
}

TEST(MemoryLayoutPA, ReadWriteLength2Dimension1) {
    MemoryLayoutPA<int, 1> layout;
    layout.ensure(1);
    layout(0, 0) = 11;
    layout(1, 0) = 21;
    EXPECT_EQ(layout(0, 0), 11);
    EXPECT_EQ(layout(1, 0), 21);
}