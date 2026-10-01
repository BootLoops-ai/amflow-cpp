// SPDX-License-Identifier: MIT
// qft::jintegral — integral identifier "j[family, n1, n2, ...]" and the
// sector algebra around it.
//
//
// Mirrors AMFlow.m lines 357-395:
//
//   ToJ / FromJ / JSector / JProp / JDot / JRank / IntegralWeight
//   SortIntegrals
//   GetTopSector / GetTopPosition
//   IntegralQ / CheckIntegrals
//   GetSector / SubSectorQ / TopSectorQ / GetTopSectorList / SplitTarget
//
// Mathematica uses the unevaluated head `j[fam, n1, n2, ...]` to represent
// an integral.  The C++ port stores the family name as a string and the
// index tuple as `vector<long>`.

#ifndef AMFLOW_QFT_JINTEGRAL_HPP
#define AMFLOW_QFT_JINTEGRAL_HPP

#include <cstddef>
#include <iosfwd>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace amflow::qft {

// ---------------------------------------------------------------------------
//  JIntegral
// ---------------------------------------------------------------------------

class JIntegral {
public:
    JIntegral() = default;
    JIntegral(std::string family, std::vector<long> indices);

    const std::string&        family()  const noexcept { return family_; }
    const std::vector<long>&  indices() const noexcept { return indices_; }
    std::size_t               n_indices() const noexcept { return indices_.size(); }

    long n_props() const noexcept;
    long n_dots()  const noexcept;
    long rank()    const noexcept;

    JIntegral sector() const;
    std::vector<int> sector_pattern() const;

    // (-n_props, -n_dots, -rank, min(indices))
    std::tuple<long, long, long, long> sort_weight() const noexcept;

    bool operator==(const JIntegral& other) const noexcept;
    bool operator!=(const JIntegral& other) const noexcept { return !(*this == other); }
    bool operator<(const JIntegral& other) const;

    std::string to_string() const;

private:
    std::string         family_;
    std::vector<long>   indices_;
};

std::ostream& operator<<(std::ostream& os, const JIntegral& j);

// ---------------------------------------------------------------------------
//  Sector geometry helpers
// ---------------------------------------------------------------------------

bool is_subsector(const std::vector<int>& sec1,
                  const std::vector<int>& sec2);

bool is_top_sector(const std::vector<int>& sec,
                   const std::vector<std::vector<int>>& list);

std::vector<int> get_top_sector(const std::vector<JIntegral>& jints);

std::vector<std::size_t> get_top_position(const std::vector<JIntegral>& jints);

std::vector<std::vector<int>>
get_top_sector_list(const std::vector<JIntegral>& jlist);

std::vector<std::vector<JIntegral>>
split_target(const std::vector<JIntegral>& jlist);

std::vector<JIntegral>
sort_integrals(std::vector<JIntegral> jlist);

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_JINTEGRAL_HPP
