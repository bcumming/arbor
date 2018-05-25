//#include <iostream>
#include <type_traits>

#include "../gtest.h"
#include "common.hpp"

#include <epoch.hpp>
#include <epoch_buffer.hpp>

using namespace arb;

TEST(epoch_buffer, index2) {
    epoch_buffer<int, 2> b;
    b[0] = 42;
    b[1] = 76;

    EXPECT_EQ(b[-2], 42);
    EXPECT_EQ(b[0], 42);
    EXPECT_EQ(b[2], 42);
    EXPECT_EQ(b[4], 42);
    EXPECT_EQ(b[6], 42);

    EXPECT_EQ(b[-1], 76);
    EXPECT_EQ(b[1], 76);
    EXPECT_EQ(b[3], 76);
    EXPECT_EQ(b[5], 76);
    EXPECT_EQ(b[7], 76);
}

TEST(epoch_buffer, index3) {
    epoch_buffer<int, 3> b;
    b[0] = 42;
    b[1] = 76;
    b[2] = 17;

    EXPECT_EQ(b[-3], 42);
    EXPECT_EQ(b[-2], 76);
    EXPECT_EQ(b[-1], 17);

    EXPECT_EQ(b[0], 42);
    EXPECT_EQ(b[1], 76);
    EXPECT_EQ(b[2], 17);

    EXPECT_EQ(b[3], 42);
    EXPECT_EQ(b[4], 76);
    EXPECT_EQ(b[5], 17);

    EXPECT_EQ(b[6], 42);
    EXPECT_EQ(b[7], 76);
    EXPECT_EQ(b[8], 17);
}
