// SPDX-License-Identifier: MIT
// numeric::log — implementation.
//
// amflow::numeric::log.

#include "amflow/numeric/log.hpp"

#include <cstdlib>

namespace amflow::numeric::log {

bool trace_enabled(const char* env_name) {
    const char* v = std::getenv(env_name);
    return v != nullptr && v[0] != '\0';
}

}  // namespace amflow::numeric::log
