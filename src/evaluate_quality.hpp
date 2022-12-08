#pragma once

#include "quality.hpp"

#include <cstddef>
#include <filesystem>
#include <vector>

bool verify(PushLogType & push_log, PopLogType const& pop_log);
bool evaluate(PushLogType& push_log, PopLogType const& pop_log, std::filesystem::path const& rank_file,
              std::filesystem::path const& delay_file);
