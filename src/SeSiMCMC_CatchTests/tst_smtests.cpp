#include <catch2/catch_test_macros.hpp>
//testing Catch2 installation
TEST_CASE("My first test with Catch2", "[fancy]")
{
    REQUIRE(0 == 0);
}
#include "atgc.hpp"

TEST_CASE("atgc")
{
    REQUIRE(1 == Atgc::atgc2ushort('a'));
}
