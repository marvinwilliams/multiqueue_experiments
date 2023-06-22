#pragma once

#ifdef WITH_PAPI
#include <papi.h>
#endif

#include <sstream>
#include <string>

static std::string get_build_info() {
    std::stringstream ss;
    ss << "Built on " << __DATE__ << ' ' << __TIME__ << " with\n";
#if defined(__clang__)
    ss << "  Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
    ss << "  GCC " << __VERSION__ << '\n';
#else
    ss << "  Unknown compiler\n";
#endif
#ifdef NDEBUG
    ss << "  NDEBUG defined\n";
#else
    ss << "  NDEBUG not defined\n";
#endif
#if defined WITH_PAPI
    ss << "  PAPI " << PAPI_VER_CURRENT << '\n';
#else
    ss << "  PAPI unsupported\n";
#endif
    return ss.str();
}
