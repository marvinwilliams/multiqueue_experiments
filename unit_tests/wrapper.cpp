#include "wrapper/capq.hpp"
#include "wrapper/klsm.hpp"
#include "wrapper/linden.hpp"
#include "wrapper/spraylist.hpp"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("capq", "[wrapper]") {
  using pq_t = wrapper::Capq<true, true, true>;
  pq_t capq{};
  auto h = capq.get_handle();
  pq_t::value_type retval;

  REQUIRE(!capq.try_delete_min(h, retval));

  retval.first = capq.min_valid_key;
  retval.second = 1;
  capq.push(h, retval);
  retval.first = capq.min_valid_key + 1;
  retval.second = 0;
  do {
  } while (!capq.try_delete_min(h, retval));
  REQUIRE(retval.first == capq.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = capq.max_valid_key;
  retval.second = 2;
  capq.push(h, retval);
  retval.first = capq.min_valid_key;
  retval.second = 0;
  do {
  } while (!capq.try_delete_min(h, retval));
  REQUIRE(retval.first == capq.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!capq.try_delete_min(h, retval));
}

TEST_CASE("klsm", "[wrapper]") {
  using pq_t = wrapper::Klsm<unsigned long, unsigned long, 1024>;
  pq_t klsm{};
  auto h = klsm.get_handle();
  pq_t::value_type retval;

  REQUIRE(!klsm.try_delete_min(h, retval));

  retval.first = klsm.min_valid_key;
  retval.second = 1;
  klsm.push(h, retval);
  retval.first = klsm.min_valid_key + 3;
  retval.second = 0;
  do {
  } while (!klsm.try_delete_min(h, retval));
  REQUIRE(retval.first == klsm.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = klsm.max_valid_key;
  retval.second = 2;
  klsm.push(h, retval);
  retval.first = klsm.min_valid_key;
  retval.second = 0;
  do {
  } while (!klsm.try_delete_min(h, retval));
  REQUIRE(retval.first == klsm.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!klsm.try_delete_min(h, retval));
}

TEST_CASE("linden", "[wrapper]") {
  using pq_t = wrapper::Linden;
  pq_t linden{};
  auto h = linden.get_handle();
  wrapper::Linden::value_type retval;

  REQUIRE(!linden.try_delete_min(h, retval));

  retval.first = linden.min_valid_key;
  retval.second = 1;
  linden.push(h, retval);
  retval.first = linden.min_valid_key + 1;
  retval.second = 0;
  do {
  } while (!linden.try_delete_min(h, retval));
  REQUIRE(retval.first == linden.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = linden.max_valid_key;
  retval.second = 2;
  linden.push(h, retval);
  retval.first = linden.min_valid_key;
  retval.second = 0;
  do {
  } while (!linden.try_delete_min(h, retval));
  REQUIRE(retval.first == linden.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!linden.try_delete_min(h, retval));
}

TEST_CASE("spraylist", "[wrapper]") {
  using pq_t = wrapper::Spraylist;
  pq_t spraylist{};
  spraylist.init_thread(1);
  auto h = spraylist.get_handle();
  wrapper::Spraylist::value_type retval;

  REQUIRE(!spraylist.try_delete_min(h, retval));

  retval.first = spraylist.min_valid_key;
  retval.second = 1;
  spraylist.push(h, retval);
  retval.first = spraylist.min_valid_key + 1;
  retval.second = 0;
  do {
  } while (!spraylist.try_delete_min(h, retval));
  REQUIRE(retval.first == spraylist.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = spraylist.max_valid_key;
  retval.second = 2;
  spraylist.push(h, retval);
  retval.first = spraylist.min_valid_key;
  retval.second = 0;
  do {
  } while (!spraylist.try_delete_min(h, retval));
  REQUIRE(retval.first == spraylist.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!spraylist.try_delete_min(h, retval));
}

