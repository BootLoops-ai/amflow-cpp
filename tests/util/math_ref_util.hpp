// SPDX-License-Identifier: MIT
// Shared helpers for MMA-reference differential tests.
//
// The math_ref JSON files at tests/data/math_ref/*.json contain Laurent
// expansions emitted by Mathematica.  These helpers parse them and
// evaluate them for comparison with C++ pipeline output.

#ifndef AMFLOW_REFACTOR_TESTS_MATH_REF_UTIL_HPP
#define AMFLOW_REFACTOR_TESTS_MATH_REF_UTIL_HPP

#include <cctype>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::refactor_test {

namespace fs = std::filesystem;
using nlohmann::json;

inline bool kira_available() {
    return fs::exists("/usr/local/bin/kira") &&
           fs::exists("/usr/share/Ferl7/fer64");
}

inline fs::path source_root() {
    return fs::path(AMFLOW_SOURCE_DIR);
}

inline double acb_real_mid(const numeric::AcbValue& value) {
    return arf_get_d(arb_midref(acb_realref(value.raw())), ARF_RND_NEAR);
}

inline double acb_imag_mid(const numeric::AcbValue& value) {
    return arf_get_d(arb_midref(acb_imagref(value.raw())), ARF_RND_NEAR);
}

inline numeric::AcbValue acb_from_complex(const std::complex<double>& value) {
    numeric::AcbValue out;
    out.set_d_d(value.real(), value.imag());
    return out;
}

inline std::string strip_spaces(std::string text) {
    std::string out;
    out.reserve(text.size());
    for (unsigned char ch : text) {
        if (!std::isspace(ch)) out.push_back(static_cast<char>(ch));
    }
    return out;
}

inline std::string strip_mma_precision(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] != '`') {
            out.push_back(text[i]);
            continue;
        }
        ++i;
        while (i < text.size()) {
            const unsigned char ch = static_cast<unsigned char>(text[i]);
            if (!std::isdigit(ch) && text[i] != '.') break;
            ++i;
        }
        --i;
    }
    return out;
}

inline std::string mma_numeric_to_arb_string(std::string text) {
    text = strip_spaces(strip_mma_precision(std::move(text)));
    std::string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '*' && i + 1 < text.size() && text[i + 1] == '^') {
            out.push_back('e');
            ++i;
            continue;
        }
        out.push_back(text[i]);
    }
    if (!out.empty() && out.front() == '+') out.erase(out.begin());
    return out;
}

inline std::complex<double> parse_mma_complex(std::string text) {
    text = strip_spaces(strip_mma_precision(std::move(text)));
    if (text.empty()) {
        throw std::runtime_error("empty Mathematica numeric literal");
    }
    if (text.find('I') == std::string::npos) {
        return {std::stod(text), 0.0};
    }

    if (text.size() >= 2 && text.substr(text.size() - 2) == "*I") {
        text.erase(text.size() - 2);
    } else if (!text.empty() && text.back() == 'I') {
        text.pop_back();
    } else {
        throw std::runtime_error("unsupported Mathematica complex literal: " + text);
    }

    std::size_t sep = std::string::npos;
    for (std::size_t i = 1; i < text.size(); ++i) {
        if ((text[i] == '+' || text[i] == '-') &&
            text[i - 1] != 'e' && text[i - 1] != 'E') {
            sep = i;
        }
    }
    if (sep == std::string::npos) {
        return {0.0, std::stod(text)};
    }
    return {std::stod(text.substr(0, sep)), std::stod(text.substr(sep))};
}

inline std::complex<double> eval_laurent_json(const json& laurent, double eps) {
    std::complex<double> sum(0.0, 0.0);
    for (const auto& term : laurent) {
        const long order = term.at(0).get<long>();
        const auto coef = parse_mma_complex(term.at(1).get<std::string>());
        sum += coef * std::pow(eps, static_cast<double>(order));
    }
    return sum;
}

inline std::map<std::string, std::complex<double>>
load_math_ref_values(const fs::path& path, double eps) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open math_ref JSON: " + path.string());
    }
    json doc;
    in >> doc;
    std::map<std::string, std::complex<double>> out;
    for (const auto& entry : doc.at("results")) {
        out.emplace(entry.at("integral").get<std::string>(),
                     eval_laurent_json(entry.at("laurent"), eps));
    }
    return out;
}

inline std::map<long, std::complex<double>>
load_math_ref_laurent(const fs::path& path, const std::string& key) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open math_ref JSON: " + path.string());
    }
    json doc;
    in >> doc;
    for (const auto& entry : doc.at("results")) {
        if (entry.at("integral").get<std::string>() != key) continue;
        std::map<long, std::complex<double>> out;
        for (const auto& term : entry.at("laurent")) {
            out.emplace(term.at(0).get<long>(),
                         parse_mma_complex(term.at(1).get<std::string>()));
        }
        return out;
    }
    throw std::runtime_error("math_ref integral not found: " + key);
}

inline std::string j_to_math_ref_key(const qft::JIntegral& j) {
    std::string out = "j[" + j.family();
    for (long idx : j.indices()) {
        out += ", ";
        out += std::to_string(idx);
    }
    out += "]";
    return out;
}

inline const numeric::AcbValue&
lookup_solution_value(const pipeline::AMFSystem& sys,
                       const pipeline::AMFSystemSolution& sol,
                       const qft::JIntegral& target) {
    if (!sol.global_preferred.empty()) {
        for (std::size_t i = 0; i < sol.global_preferred.size(); ++i) {
            if (sol.global_preferred[i] == target) return sol.global_values[i];
        }
    }
    const auto& preferred = sys.preferred();
    for (std::size_t i = 0; i < preferred.size(); ++i) {
        if (preferred[i] == target) return sol.master_values[i];
    }
    throw std::runtime_error("target integral not found in AMFSystem solution: "
                              + target.to_string());
}

}  // namespace amflow::refactor_test

#endif  // AMFLOW_REFACTOR_TESTS_MATH_REF_UTIL_HPP
