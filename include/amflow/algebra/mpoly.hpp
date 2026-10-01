// SPDX-License-Identifier: MIT
// algebra::mpoly — multivariate polynomials and rational functions over ℚ
// with named-variable contexts.
//
//
// Three types:
//   * MpolyContext  — pinned, shared via shared_ptr<MpolyContext>.
//                      Variable order is set at construction.
//   * Mpoly         — wraps fmpz_mpoly_t (ℤ multivariate).
//   * Mfrac         — wraps fmpz_mpoly_q_t (ℚ multivariate; FLINT
//                      auto-reduces by gcd after every operation).
//
// All three are move-only with explicit clone(); cross-context operations
// throw on shared_ptr-identity mismatch.
//
// Sentinel-prefix variables (transparent on context narrowing — see
// algebra::project_mfrac_dropping in context_migration.hpp):
//   * __amf_*    internal substitution sentinels
//   * __feyn_*   Feynman parameters
//   * __zsq_*    squared-mass auxiliaries

#ifndef AMFLOW_ALGEBRA_MPOLY_HPP
#define AMFLOW_ALGEBRA_MPOLY_HPP

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/numeric/options.hpp"

namespace amflow::algebra {

// ---------------------------------------------------------------------------
//  MpolyContext  — pinned, shared via shared_ptr.
// ---------------------------------------------------------------------------

class MpolyContext {
public:
    explicit MpolyContext(std::vector<std::string> var_names);
    ~MpolyContext();

    // Pinned: FLINT's ctx struct holds internal pointers that copy/move
    // would invalidate.
    MpolyContext(const MpolyContext&)            = delete;
    MpolyContext& operator=(const MpolyContext&) = delete;
    MpolyContext(MpolyContext&&)                 = delete;
    MpolyContext& operator=(MpolyContext&&)      = delete;

    long              n_vars() const noexcept    { return static_cast<long>(var_names_.size()); }
    const std::vector<std::string>& var_names() const noexcept { return var_names_; }

    long              var_index(const std::string& name) const;   // -1 if not found
    const std::string& var_name(long idx) const;

    fmpz_mpoly_ctx_struct*       raw()       noexcept { return ctx_; }
    const fmpz_mpoly_ctx_struct* raw() const noexcept { return ctx_; }

    // C-style variable name array for FLINT routines.
    const char* const* var_names_cstr() const noexcept { return var_cstr_.data(); }

private:
    std::vector<std::string>  var_names_;
    std::vector<const char*>  var_cstr_;          // points into var_names_
    fmpz_mpoly_ctx_struct*    ctx_;
};

// ---------------------------------------------------------------------------
//  Mpoly  — ℤ[x1, ..., xn].  Move-only.
// ---------------------------------------------------------------------------

class Mpoly {
public:
    Mpoly() noexcept;                                      // null (no ctx)
    explicit Mpoly(std::shared_ptr<MpolyContext> ctx);     // zero polynomial

    Mpoly(const Mpoly&)            = delete;
    Mpoly& operator=(const Mpoly&) = delete;

    Mpoly(Mpoly&& other) noexcept;
    Mpoly& operator=(Mpoly&& other) noexcept;
    ~Mpoly();

    Mpoly clone() const;

    // ---- factory helpers ----
    static Mpoly zero(std::shared_ptr<MpolyContext> ctx);
    static Mpoly one (std::shared_ptr<MpolyContext> ctx);
    static Mpoly constant(std::shared_ptr<MpolyContext> ctx, long c);
    static Mpoly constant(std::shared_ptr<MpolyContext> ctx, const fmpz_t c);
    static Mpoly variable(std::shared_ptr<MpolyContext> ctx, long var_idx);
    static Mpoly variable(std::shared_ptr<MpolyContext> ctx, const std::string& name);
    static Mpoly monomial(std::shared_ptr<MpolyContext> ctx,
                          const fmpz_t c,
                          const std::vector<unsigned long>& exponents);

    // Parse from a string like "1+2*eta-3*eta^2+s*msq".  Variable names must
    // be in ctx.  Whitespace ignored.
    //
    // Internally accepts fractional coefficients via fmpq_mpoly parser then
    // converts to integer.  When `out_scale != nullptr`, fractional input is
    // accepted and the common denominator is written to `out_scale`; when
    // null, throws on fractional coefficient.
    static Mpoly from_string(std::shared_ptr<MpolyContext> ctx,
                             const std::string& expr,
                             fmpz_t out_scale = nullptr);

    // ---- accessors ----
    const std::shared_ptr<MpolyContext>& ctx() const noexcept { return ctx_; }
    fmpz_mpoly_struct*       raw()        noexcept { return p_; }
    const fmpz_mpoly_struct* raw() const noexcept { return p_; }

    bool is_zero() const noexcept;
    bool is_one()  const noexcept;
    bool is_constant() const noexcept;

    long total_degree() const noexcept;          // -1 for zero
    long degree(long var_idx) const noexcept;    // -1 for zero

    // ---- arithmetic ----
    Mpoly  operator-() const;
    Mpoly& operator+=(const Mpoly& other);
    Mpoly& operator-=(const Mpoly& other);
    Mpoly& operator*=(const Mpoly& other);

    // Returns true and writes a/b to out iff b divides a exactly.  Required
    // for Bareiss elimination in MpolyMatrix::det.  All three Mpolys share ctx.
    static bool exact_divide(Mpoly& out, const Mpoly& a, const Mpoly& b);

    // Returns gcd(a, b).  All three Mpolys share ctx.  Used by the
    // batched-lcm-sum path in `ibp::diffeq` matrix assembly.
    static Mpoly gcd(const Mpoly& a, const Mpoly& b);

    Mpoly& multiply_by_si(long c);
    Mpoly& multiply_by_fmpz(const fmpz_t c);

    Mpoly derivative(long var_idx) const;

    // Coefficient of x_{vars[0]}^{exps[0]} * x_{vars[1]}^{exps[1]} * ...
    // (Mathematica's Coefficient[expr, Loop[[j]] * Loop[[k]]] etc.)
    Mpoly coeff_of(const std::vector<long>& vars,
                   const std::vector<unsigned long>& exps) const;
    Mpoly coeff_of(long var, unsigned long power) const;

    // x_var -> integer value.
    Mpoly substitute(long var_idx, const fmpz_t value) const;
    Mpoly substitute(long var_idx, long value) const;

    // Evaluate at fully-numeric integer inputs.  values[i] points to the
    // fmpz to substitute for variable i.  Returns 1 on success; 0 on FLINT
    // failure.
    int evaluate(fmpz * const * values, fmpz_t out) const;

    // Evaluate to acb (numerical).  values is a contiguous acb array of
    // length n_vars.  out must be pre-initialised.
    void evaluate_acb(acb_ptr out,
                      acb_srcptr values,
                      long prec = numeric::working_prec_bits()) const;

    std::string to_string() const;

    bool operator==(const Mpoly& other) const;
    bool operator!=(const Mpoly& other) const { return !(*this == other); }

private:
    void check_same_ctx(const Mpoly& other) const;
    std::shared_ptr<MpolyContext> ctx_;
    fmpz_mpoly_struct*            p_;
};

Mpoly operator+(const Mpoly& a, const Mpoly& b);
Mpoly operator-(const Mpoly& a, const Mpoly& b);
Mpoly operator*(const Mpoly& a, const Mpoly& b);

std::ostream& operator<<(std::ostream& os, const Mpoly& p);

// ---------------------------------------------------------------------------
//  Mfrac  — ℚ(x1, ..., xn).  Move-only.  Auto-reduced after every op.
// ---------------------------------------------------------------------------

class Mfrac {
public:
    Mfrac() noexcept;                                      // null (no ctx)
    explicit Mfrac(std::shared_ptr<MpolyContext> ctx);     // zero
    Mfrac(Mpoly num, Mpoly den);                           // throws on zero den

    Mfrac(const Mfrac&)            = delete;
    Mfrac& operator=(const Mfrac&) = delete;

    Mfrac(Mfrac&& other) noexcept;
    Mfrac& operator=(Mfrac&& other) noexcept;
    ~Mfrac();

    Mfrac clone() const;

    static Mfrac zero(std::shared_ptr<MpolyContext> ctx);
    static Mfrac one (std::shared_ptr<MpolyContext> ctx);
    static Mfrac from_mpoly(Mpoly p);
    static Mfrac from_si(std::shared_ptr<MpolyContext> ctx, long c);
    static Mfrac from_fmpq(std::shared_ptr<MpolyContext> ctx, const fmpq_t q);

    // ---- accessors ----
    const std::shared_ptr<MpolyContext>& ctx() const noexcept { return ctx_; }
    fmpz_mpoly_q_struct*        raw()        noexcept { return q_; }
    const fmpz_mpoly_q_struct*  raw() const  noexcept { return q_; }

    Mpoly numerator()   const;
    Mpoly denominator() const;

    bool is_zero() const noexcept;
    bool is_one()  const noexcept;

    // True iff the (canonical) denominator is a constant polynomial.  In
    // other words this Mfrac is just a polynomial (with possibly rational
    // coefficients).
    bool is_polynomial() const noexcept;

    // ---- arithmetic ----
    Mfrac  operator-() const;
    Mfrac& operator+=(const Mfrac& other);
    Mfrac& operator-=(const Mfrac& other);
    Mfrac& operator*=(const Mfrac& other);
    Mfrac& operator/=(const Mfrac& other);

    Mfrac& multiply_by_mpoly(const Mpoly& p);

    Mfrac derivative(long var_idx) const;

    Mfrac substitute(long var_idx, const fmpq_t value) const;
    Mfrac substitute(long var_idx, long value) const;

    // Evaluate fully at rational inputs.  values is an array of length
    // n_vars.  Returns false (and leaves out unchanged) if the denominator
    // vanishes; true on success.
    bool evaluate(fmpq * const * values, fmpq_t out) const;

    // Evaluate to acb (numerical).  values length = n_vars.  out must be
    // pre-initialised.
    void evaluate_acb(acb_ptr out,
                      acb_srcptr values,
                      long prec = numeric::working_prec_bits()) const;

    std::string to_string() const;

    bool operator==(const Mfrac& other) const;
    bool operator!=(const Mfrac& other) const { return !(*this == other); }

private:
    void check_same_ctx(const Mfrac& other) const;
    std::shared_ptr<MpolyContext> ctx_;
    fmpz_mpoly_q_struct*          q_;
};

Mfrac operator+(const Mfrac& a, const Mfrac& b);
Mfrac operator-(const Mfrac& a, const Mfrac& b);
Mfrac operator*(const Mfrac& a, const Mfrac& b);
Mfrac operator/(const Mfrac& a, const Mfrac& b);

std::ostream& operator<<(std::ostream& os, const Mfrac& f);

}  // namespace amflow::algebra

#endif  // AMFLOW_ALGEBRA_MPOLY_HPP
