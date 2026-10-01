// SPDX-License-Identifier: MIT
// qft::jintegral — implementation.
//

#include "amflow/qft/jintegral.hpp"

#include <algorithm>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace amflow::qft {

// ===========================================================================
//  JIntegral
// ===========================================================================

JIntegral::JIntegral(std::string family, std::vector<long> indices)
    : family_(std::move(family)), indices_(std::move(indices)) {
    if (family_.empty()) {
        throw std::invalid_argument("JIntegral: family must be non-empty");
    }
}

long JIntegral::n_props() const noexcept {
    long n = 0;
    for (long i : indices_) if (i > 0) ++n;
    return n;
}

long JIntegral::n_dots() const noexcept {
    long s = 0;
    for (long i : indices_) if (i > 0) s += (i - 1);
    return s;
}

long JIntegral::rank() const noexcept {
    long s = 0;
    for (long i : indices_) if (i < 0) s += -i;
    return s;
}

JIntegral JIntegral::sector() const {
    std::vector<long> ind(indices_.size());
    for (std::size_t i = 0; i < indices_.size(); ++i) {
        ind[i] = (indices_[i] > 0) ? 1 : 0;
    }
    return JIntegral(family_, std::move(ind));
}

std::vector<int> JIntegral::sector_pattern() const {
    std::vector<int> out(indices_.size());
    for (std::size_t i = 0; i < indices_.size(); ++i) {
        out[i] = (indices_[i] > 0) ? 1 : 0;
    }
    return out;
}

std::tuple<long, long, long, long> JIntegral::sort_weight() const noexcept {
    long min_idx = indices_.empty()
                 ? std::numeric_limits<long>::max()
                 : indices_[0];
    for (long i : indices_) if (i < min_idx) min_idx = i;
    return {-n_props(), -n_dots(), -rank(), min_idx};
}

bool JIntegral::operator==(const JIntegral& other) const noexcept {
    return family_ == other.family_ && indices_ == other.indices_;
}

bool JIntegral::operator<(const JIntegral& other) const {
    auto a = sort_weight();
    auto b = other.sort_weight();
    if (a != b) return a < b;
    if (family_ != other.family_) return family_ < other.family_;
    return indices_ < other.indices_;
}

std::string JIntegral::to_string() const {
    std::ostringstream oss;
    oss << "j[" << family_;
    for (long i : indices_) {
        oss << ", " << i;
    }
    oss << "]";
    return oss.str();
}

std::ostream& operator<<(std::ostream& os, const JIntegral& j) {
    return os << j.to_string();
}

// ===========================================================================
//  Sector geometry
// ===========================================================================

bool is_subsector(const std::vector<int>& sec1, const std::vector<int>& sec2) {
    if (sec1.size() != sec2.size()) return false;
    for (std::size_t i = 0; i < sec1.size(); ++i) {
        if (sec1[i] - sec2[i] > 0) return false;
    }
    return true;
}

bool is_top_sector(const std::vector<int>& sec,
                   const std::vector<std::vector<int>>& list) {
    long count = 0;
    for (const auto& other : list) {
        if (is_subsector(sec, other)) ++count;
    }
    return count == 1;
}

std::vector<int> get_top_sector(const std::vector<JIntegral>& jints) {
    if (jints.empty()) return {};
    std::size_t n = jints.front().n_indices();
    for (const auto& j : jints) {
        if (j.n_indices() != n) {
            throw std::invalid_argument(
                "get_top_sector: integrals have different n_indices");
        }
    }
    std::vector<int> out(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        long max_v = std::numeric_limits<long>::min();
        for (const auto& j : jints) {
            if (j.indices()[k] > max_v) max_v = j.indices()[k];
        }
        out[k] = (max_v > 0) ? 1 : 0;
    }
    return out;
}

std::vector<std::size_t>
get_top_position(const std::vector<JIntegral>& jints) {
    auto sec = get_top_sector(jints);
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < sec.size(); ++i) {
        if (sec[i] != 0) out.push_back(i);
    }
    return out;
}

std::vector<std::vector<int>>
get_top_sector_list(const std::vector<JIntegral>& jlist) {
    std::vector<std::vector<int>> sectors;
    for (const auto& j : jlist) {
        auto p = j.sector_pattern();
        bool seen = false;
        for (const auto& s : sectors) {
            if (s == p) { seen = true; break; }
        }
        if (!seen) sectors.push_back(std::move(p));
    }
    std::vector<std::vector<int>> out;
    for (const auto& s : sectors) {
        if (is_top_sector(s, sectors)) out.push_back(s);
    }
    return out;
}

std::vector<std::vector<JIntegral>>
split_target(const std::vector<JIntegral>& jlist) {
    auto top_sectors = get_top_sector_list(jlist);

    std::vector<bool> claimed(jlist.size(), false);
    std::vector<std::vector<JIntegral>> out;
    out.reserve(top_sectors.size());

    for (const auto& top : top_sectors) {
        std::vector<JIntegral> chunk;
        for (std::size_t k = 0; k < jlist.size(); ++k) {
            if (claimed[k]) continue;
            if (is_subsector(jlist[k].sector_pattern(), top)) {
                chunk.push_back(jlist[k]);
                claimed[k] = true;
            }
        }
        out.push_back(std::move(chunk));
    }

    for (std::size_t k = 0; k < jlist.size(); ++k) {
        if (!claimed[k]) {
            throw std::runtime_error(
                "split_target: integral "
                + jlist[k].to_string()
                + " could not be assigned to any top sector");
        }
    }
    return out;
}

std::vector<JIntegral>
sort_integrals(std::vector<JIntegral> jlist) {
    std::stable_sort(jlist.begin(), jlist.end());
    return jlist;
}

}  // namespace amflow::qft
