#include "wrapper/capq.hpp"
#include "wrapper/klsm.hpp"
#include "wrapper/linden.hpp"
#include "wrapper/spraylist.hpp"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("capq", "[wrapper]") {
  using pq_t = wrapper::Capq<true, true, true>;
  pq_t capq(1);
  auto h = capq.get_handle();
  pq_t::value_type retval;

  REQUIRE(!capq.try_pop(retval));

  retval.first = pq_t::min_valid_key;
  retval.second = 1;
  h.push(retval);
  retval.first = pq_t::min_valid_key + 1;
  retval.second = 0;
  do {
  } while (!capq.try_pop(retval));
  REQUIRE(retval.first == capq.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = pq_t::max_valid_key;
  retval.second = 2;
  h.push(retval);
  retval.first = pq_t::min_valid_key;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == capq.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!h.try_pop(retval));
}

TEST_CASE("klsm", "[wrapper]") {
  using pq_t = wrapper::Klsm<unsigned long, unsigned long, 1024>;
  pq_t klsm(1);
  auto h = klsm.get_handle();
  pq_t::value_type retval;

  REQUIRE(!h.try_pop(retval));

  retval.first = 0;
  retval.second = 1;
  h.push(retval);
  retval.first = 3;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == 0);
  REQUIRE(retval.second == 1);

  retval.first = 10;
  retval.second = 2;
  h.push(retval);
  retval.first = 0;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == 10);
  REQUIRE(retval.second == 2);

  REQUIRE(!h.try_pop(retval));
}

TEST_CASE("linden", "[wrapper]") {
  using pq_t = wrapper::Linden;
  pq_t linden(1);
  auto h = linden.get_handle();
  wrapper::Linden::value_type retval;

  REQUIRE(!h.try_pop(retval));

  retval.first = pq_t::min_valid_key;
  retval.second = 1;
  h.push(retval);
  retval.first = pq_t::min_valid_key + 1;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == linden.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = pq_t::max_valid_key;
  retval.second = 2;
  h.push(retval);
  retval.first = pq_t::min_valid_key;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == linden.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!h.try_pop(retval));
}

TEST_CASE("spraylist", "[wrapper]") {
  using pq_t = wrapper::Spraylist;
  pq_t spraylist(1);
  auto h = spraylist.get_handle();
  wrapper::Spraylist::value_type retval;

  REQUIRE(!h.try_pop(retval));

  retval.first = pq_t::min_valid_key;
  retval.second = 1;
  h.push(retval);
  retval.first = pq_t::min_valid_key + 1;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == spraylist.min_valid_key);
  REQUIRE(retval.second == 1);

  retval.first = pq_t::max_valid_key;
  retval.second = 2;
  h.push(retval);
  retval.first = pq_t::min_valid_key;
  retval.second = 0;
  do {
  } while (!h.try_pop(retval));
  REQUIRE(retval.first == spraylist.max_valid_key);
  REQUIRE(retval.second == 2);

  REQUIRE(!h.try_pop(retval));
}

