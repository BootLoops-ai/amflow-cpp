// SPDX-License-Identifier: MIT
// GoogleTest entry point for amflow_tests.

#include <gtest/gtest.h>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
