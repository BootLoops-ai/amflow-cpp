// SPDX-License-Identifier: MIT
// amflow::api::run_json — JSON-driven entry point.
//
// See run_json.cpp for the supported input schemas
// ("amflow", "solve_integrals", "black_box_amflow").

#ifndef AMFLOW_API_RUN_JSON_HPP
#define AMFLOW_API_RUN_JSON_HPP

#include <nlohmann/json.hpp>

namespace amflow::api {

// Apply input.options (if present), dispatch on input.mode, return the
// usual {result, options} envelope.  Throws std::runtime_error on any
// schema/parse/compute failure; callers (cli) are expected to catch.
nlohmann::json run_json(const nlohmann::json& input);

}  // namespace amflow::api

#endif  // AMFLOW_API_RUN_JSON_HPP
