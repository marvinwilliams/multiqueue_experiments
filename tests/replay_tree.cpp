#include "util/replay_tree.hpp"
#include "catch2/catch_test_macros.hpp"

#include <iostream>

struct extract_key {
    static int const& get(int const& i) {
        return i;
    }
};

TEST_CASE("replay_tree insert", "[replay_tree]") {
    ReplayTree<int, int, extract_key> replay_tree;
    REQUIRE(replay_tree.empty());
    for (int i = 0; i < 500'000; ++i) {
        replay_tree.insert(i);
    }
    for (int i = 0; i < 500'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == i);
        auto it = replay_tree.find(i);
        REQUIRE(it != replay_tree.end());
    }
    for (int i = 0; i < 500'000; ++i) {
        replay_tree.insert(1'000'000 - i - 1);
    }
    for (int i = 0; i < 1'000'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == i);
        auto it = replay_tree.find(i);
        REQUIRE(it != replay_tree.end());
    }
}

TEST_CASE("replay_tree no increase", "[replay_tree]") {
    ReplayTree<int, int, extract_key> replay_tree;
    REQUIRE(replay_tree.empty());
    replay_tree.insert(0);
    replay_tree.insert(1);
    replay_tree.insert(2);
    replay_tree.insert(3);
    replay_tree.insert(4);
    for (int i = 0; i < 5; ++i) {
        REQUIRE(replay_tree.get_rank(i) == 0);
        auto [success, rank, delay] = replay_tree.erase_val(i);
        REQUIRE(success);
        REQUIRE(delay == 0);
    }
    REQUIRE(replay_tree.empty());
    replay_tree.insert(4);
    replay_tree.insert(3);
    replay_tree.insert(2);
    replay_tree.insert(1);
    replay_tree.insert(0);
    for (int i = 0; i < 5; ++i) {
        REQUIRE(replay_tree.get_rank(4 - i) == 4 - i);
        auto [success, rank, delay] = replay_tree.erase_val(4 - i);
        REQUIRE(success);
        REQUIRE(delay == i);
    }
    REQUIRE(replay_tree.empty());
    for (int i = 0; i < 1'000'000; ++i) {
        replay_tree.insert(i);
    }
    for (int i = 0; i < 500'000; ++i) {
        REQUIRE(replay_tree.get_rank(500'000 - i - 1) == 500'000 - i - 1);
        auto [success, rank, delay] = replay_tree.erase_val(500'000 - i - 1);
        REQUIRE(success);
        REQUIRE(delay == i);
    }
    for (int i = 0; i < 500'000; ++i) {
        replay_tree.insert(i);
    }
    for (int i = 0; i < 1'000'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == i);
    }
    for (int i = 500'000; i < 1'000'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == 500'000);
        auto [success, rank, delay] = replay_tree.erase_val(i);
        REQUIRE(success);
        REQUIRE(delay == 0);
    }
    for (int i = 0; i < 500'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == 0);
        auto [success, rank, delay] = replay_tree.erase_val(i);
        REQUIRE(success);
        REQUIRE(delay == 500'000);
    }
}

TEST_CASE("replay_tree", "[replay_tree]") {
    ReplayTree<int, int, extract_key> replay_tree;
    REQUIRE(replay_tree.empty());
    replay_tree.insert(0);
    replay_tree.insert(1);
    replay_tree.insert(2);
    replay_tree.insert(3);
    replay_tree.insert(4);
    for (int i = 0; i < 5; ++i) {
        REQUIRE(replay_tree.get_rank(i) == 0);
        auto [success, rank, delay] = replay_tree.erase_val(i);
        REQUIRE(success);
        REQUIRE(delay == 0);
    }
    REQUIRE(replay_tree.empty());
    replay_tree.insert(4);
    replay_tree.insert(3);
    replay_tree.insert(2);
    replay_tree.insert(1);
    replay_tree.insert(0);
    replay_tree.increase_global_delay();
    for (int i = 0; i < 5; ++i) {
        REQUIRE(replay_tree.get_rank(4 - i) == 4 - i);
        auto [success, rank, delay] = replay_tree.erase_val(4 - i);
        REQUIRE(success);
        REQUIRE(delay == i + 1);
    }
    REQUIRE(replay_tree.empty());
    for (int i = 0; i < 1'000'000; ++i) {
        replay_tree.insert(i);
    }
    /* replay_tree.increase_delay(500); */
    /* for (int i = 0; i < 1000; ++i) { */
    /*   auto it = replay_tree.find(i); */
    /*   auto [success, delay] = replay_tree.erase(it); */
    /*   REQUIRE(success); */
    /*   std::cout << i << ' ' << delay << '\n'; */
    /* } */
    /* for (int i = 0; i < 1000; ++i) { */
    /*   auto it = replay_tree.find(i); */
    /*   auto [success, delay] = replay_tree.erase(it); */
    /*   REQUIRE(success); */
    /*   std::cout << i << ' ' << delay << '\n'; */
    /* } */
    for (int i = 0; i < 500'000; ++i) {
        REQUIRE(replay_tree.get_rank(500'000 - i - 1) == 500'000 - i - 1);
        auto [success, rank, delay] = replay_tree.erase_val(500'000 - i - 1);
        REQUIRE(success);
        REQUIRE(delay == i);
    }
    for (int i = 0; i < 500'000; ++i) {
        replay_tree.insert(i);
    }
    /* for (int i = 500; i < 1000; ++i) { */
    /*   auto it = replay_tree.find(i); */
    /*   auto [success, delay] = replay_tree.erase(it); */
    /*   REQUIRE(success); */
    /*   std::cout << i << ' ' << delay << '\n'; */
    /* } */
    for (int i = 500'000; i < 1'000'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == 500'000);
        auto [success, rank, delay] = replay_tree.erase_val(i);
        REQUIRE(success);
        REQUIRE(delay == 0);
    }
    for (int i = 0; i < 500'000; ++i) {
        REQUIRE(replay_tree.get_rank(i) == 0);
        auto [success, rank, delay] = replay_tree.erase_val(i);
        REQUIRE(success);
        REQUIRE(delay == 500'000);
    }
}
