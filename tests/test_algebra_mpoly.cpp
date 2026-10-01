// SPDX-License-Identifier: MIT
// Tests for amflow::algebra MpolyContext / Mpoly / Mfrac.

#include <gtest/gtest.h>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace alg = amflow::algebra;

namespace {

std::shared_ptr<alg::MpolyContext> make_ctx(std::vector<std::string> names) {
    return std::make_shared<alg::MpolyContext>(std::move(names));
}

}  // namespace

// ---------------------------------------------------------------------------
//  MpolyContext
// ---------------------------------------------------------------------------

TEST(MpolyContext, ConstructionAndLookup) {
    auto ctx = make_ctx({"x", "y", "z"});
    EXPECT_EQ(ctx->n_vars(), 3);
    EXPECT_EQ(ctx->var_index("x"), 0);
    EXPECT_EQ(ctx->var_index("y"), 1);
    EXPECT_EQ(ctx->var_index("z"), 2);
    EXPECT_EQ(ctx->var_index("missing"), -1);
    EXPECT_EQ(ctx->var_name(0), "x");
    EXPECT_EQ(ctx->var_name(2), "z");
    EXPECT_THROW(ctx->var_name(3), std::out_of_range);
    EXPECT_THROW(ctx->var_name(-1), std::out_of_range);
}

TEST(MpolyContext, EmptyNamesThrows) {
    auto bad = []() {
        alg::MpolyContext c{std::vector<std::string>{}};
        (void)c;
    };
    EXPECT_THROW(bad(), std::invalid_argument);
}

TEST(MpolyContext, DuplicateNamesThrows) {
    auto bad = []() {
        alg::MpolyContext c{std::vector<std::string>{"x", "y", "x"}};
        (void)c;
    };
    EXPECT_THROW(bad(), std::invalid_argument);
}

// ---------------------------------------------------------------------------
//  Mpoly factories
// ---------------------------------------------------------------------------

TEST(Mpoly, ZeroAndOne) {
    auto ctx = make_ctx({"x"});
    auto z = alg::Mpoly::zero(ctx);
    EXPECT_TRUE(z.is_zero());
    EXPECT_FALSE(z.is_one());
    EXPECT_EQ(z.total_degree(), -1);

    auto o = alg::Mpoly::one(ctx);
    EXPECT_FALSE(o.is_zero());
    EXPECT_TRUE(o.is_one());
    EXPECT_TRUE(o.is_constant());
    EXPECT_EQ(o.total_degree(), 0);
}

TEST(Mpoly, VariableByIndexAndName) {
    auto ctx = make_ctx({"x", "y"});
    auto x_idx = alg::Mpoly::variable(ctx, 0);
    auto x_name = alg::Mpoly::variable(ctx, "x");
    EXPECT_EQ(x_idx, x_name);
    EXPECT_EQ(x_idx.degree(0), 1);
    EXPECT_EQ(x_idx.degree(1), 0);
}

TEST(Mpoly, VariableUnknownNameThrows) {
    auto ctx = make_ctx({"x"});
    EXPECT_THROW(alg::Mpoly::variable(ctx, "nope"), std::invalid_argument);
    EXPECT_THROW(alg::Mpoly::variable(ctx, 5), std::out_of_range);
}

TEST(Mpoly, ConstantFromInteger) {
    auto ctx = make_ctx({"x"});
    auto c = alg::Mpoly::constant(ctx, 42);
    EXPECT_TRUE(c.is_constant());
    EXPECT_FALSE(c.is_zero());
}

TEST(Mpoly, MonomialBuildsCorrectly) {
    auto ctx = make_ctx({"x", "y"});
    fmpz_t two;
    fmpz_init(two);
    fmpz_set_si(two, 2);
    // 2 * x^3 * y^1
    auto m = alg::Mpoly::monomial(ctx, two, {3u, 1u});
    EXPECT_EQ(m.degree(0), 3);
    EXPECT_EQ(m.degree(1), 1);
    EXPECT_EQ(m.total_degree(), 4);
    fmpz_clear(two);
}

// ---------------------------------------------------------------------------
//  Mpoly::from_string
// ---------------------------------------------------------------------------

TEST(MpolyFromString, IntegerCoefficients) {
    auto ctx = make_ctx({"x"});
    auto p = alg::Mpoly::from_string(ctx, "x^2 + 2*x + 1");
    EXPECT_EQ(p.degree(0), 2);
    EXPECT_FALSE(p.is_zero());
}

TEST(MpolyFromString, UnknownVarThrows) {
    auto ctx = make_ctx({"x"});
    EXPECT_THROW(alg::Mpoly::from_string(ctx, "x + foo"),
                 std::invalid_argument);
}

TEST(MpolyFromString, RationalCoeffWithoutScaleThrows) {
    auto ctx = make_ctx({"x"});
    EXPECT_THROW(alg::Mpoly::from_string(ctx, "1/2*x"),
                 std::invalid_argument);
}

TEST(MpolyFromString, RationalCoeffWithScaleSucceeds) {
    auto ctx = make_ctx({"x"});
    fmpz_t scale;
    fmpz_init(scale);
    auto p = alg::Mpoly::from_string(ctx, "1/2*x", scale);
    // Scale is 2; numerator is x.
    EXPECT_EQ(fmpz_get_si(scale), 2);
    EXPECT_EQ(p.degree(0), 1);
    fmpz_clear(scale);
}

// ---------------------------------------------------------------------------
//  Mpoly arithmetic
// ---------------------------------------------------------------------------

TEST(MpolyArithmetic, Addition) {
    auto ctx = make_ctx({"x"});
    auto x = alg::Mpoly::variable(ctx, 0);
    auto sum = x + x;   // 2x
    auto two_x = alg::Mpoly::from_string(ctx, "2*x");
    EXPECT_EQ(sum, two_x);
}

TEST(MpolyArithmetic, Multiplication) {
    auto ctx = make_ctx({"x", "y"});
    auto x = alg::Mpoly::variable(ctx, 0);
    auto y = alg::Mpoly::variable(ctx, 1);
    auto xy = x * y;
    EXPECT_EQ(xy.degree(0), 1);
    EXPECT_EQ(xy.degree(1), 1);
    EXPECT_EQ(xy.total_degree(), 2);
}

TEST(MpolyArithmetic, ContextMismatchThrows) {
    auto ctx_a = make_ctx({"x"});
    auto ctx_b = make_ctx({"x"});
    auto a = alg::Mpoly::variable(ctx_a, 0);
    auto b = alg::Mpoly::variable(ctx_b, 0);
    EXPECT_THROW(a + b, std::invalid_argument);
}

TEST(MpolyArithmetic, MoveSemanticsClearSource) {
    auto ctx = make_ctx({"x"});
    auto a = alg::Mpoly::variable(ctx, 0);
    auto b = std::move(a);
    EXPECT_EQ(b.degree(0), 1);
}

TEST(MpolyArithmetic, CloneIsDeepCopy) {
    auto ctx = make_ctx({"x"});
    auto a = alg::Mpoly::from_string(ctx, "2*x + 3");
    auto b = a.clone();
    auto x = alg::Mpoly::variable(ctx, 0);
    a += x;     // a = 3x + 3; b unchanged.
    EXPECT_NE(a, b);
}

TEST(MpolyArithmetic, ExactDivideSuccess) {
    auto ctx = make_ctx({"x"});
    auto x_plus_1 = alg::Mpoly::from_string(ctx, "x + 1");
    auto x_squared_minus_1 = alg::Mpoly::from_string(ctx, "x^2 - 1");
    alg::Mpoly q(ctx);
    EXPECT_TRUE(alg::Mpoly::exact_divide(q, x_squared_minus_1, x_plus_1));
    auto x_minus_1 = alg::Mpoly::from_string(ctx, "x - 1");
    EXPECT_EQ(q, x_minus_1);
}

TEST(MpolyArithmetic, ExactDivideFailureReturnsFalse) {
    auto ctx = make_ctx({"x"});
    auto x_plus_2 = alg::Mpoly::from_string(ctx, "x + 2");
    auto x_squared_minus_1 = alg::Mpoly::from_string(ctx, "x^2 - 1");
    alg::Mpoly q(ctx);
    EXPECT_FALSE(alg::Mpoly::exact_divide(q, x_squared_minus_1, x_plus_2));
}

// ---------------------------------------------------------------------------
//  Mpoly::derivative / coeff_of / substitute
// ---------------------------------------------------------------------------

TEST(Mpoly, DerivativeOfPolynomial) {
    auto ctx = make_ctx({"x"});
    auto p = alg::Mpoly::from_string(ctx, "x^3 - 2*x^2 + x");
    auto d = p.derivative(0);   // 3x^2 - 4x + 1
    auto expected = alg::Mpoly::from_string(ctx, "3*x^2 - 4*x + 1");
    EXPECT_EQ(d, expected);
}

TEST(Mpoly, CoefficientOfMonomial) {
    auto ctx = make_ctx({"x", "y"});
    // 5 * x^2 * y + 3 * x * y + 7
    auto p = alg::Mpoly::from_string(ctx, "5*x^2*y + 3*x*y + 7");
    auto coef_x_y = p.coeff_of({0, 1}, {1, 1});   // coefficient of x*y
    auto expected = alg::Mpoly::from_string(ctx, "3");
    EXPECT_EQ(coef_x_y, expected);
}

TEST(Mpoly, CoefficientDuplicateVarsThrows) {
    auto ctx = make_ctx({"x"});
    auto p = alg::Mpoly::from_string(ctx, "x^2 + x");
    EXPECT_THROW(p.coeff_of({0, 0}, {1, 1}), std::invalid_argument);
}

TEST(Mpoly, SubstituteRationalIntoVariable) {
    auto ctx = make_ctx({"x", "y"});
    auto p = alg::Mpoly::from_string(ctx, "x^2 + y");
    auto p_x_eq_3 = p.substitute(0, 3);   // 9 + y
    auto expected = alg::Mpoly::from_string(ctx, "y + 9");
    EXPECT_EQ(p_x_eq_3, expected);
}

TEST(Mpoly, EvaluateAcb) {
    auto ctx = make_ctx({"x"});
    auto p = alg::Mpoly::from_string(ctx, "x^2 + 1");

    amflow::numeric::AcbValue x_val;
    x_val.set_si(2);
    amflow::numeric::AcbValue out;
    p.evaluate_acb(out.raw(), x_val.raw());

    amflow::numeric::AcbValue expected;
    expected.set_si(5);
    EXPECT_NE(acb_eq(out.raw(), expected.raw()), 0);
}

// ---------------------------------------------------------------------------
//  Mfrac
// ---------------------------------------------------------------------------

TEST(Mfrac, ZeroAndOne) {
    auto ctx = make_ctx({"x"});
    auto z = alg::Mfrac::zero(ctx);
    EXPECT_TRUE(z.is_zero());
    EXPECT_FALSE(z.is_one());
    EXPECT_TRUE(z.is_polynomial());

    auto o = alg::Mfrac::one(ctx);
    EXPECT_FALSE(o.is_zero());
    EXPECT_TRUE(o.is_one());
}

TEST(Mfrac, ZeroDenominatorThrows) {
    auto ctx = make_ctx({"x"});
    auto num = alg::Mpoly::variable(ctx, 0);
    auto den = alg::Mpoly::zero(ctx);
    EXPECT_THROW(alg::Mfrac(std::move(num), std::move(den)),
                 std::invalid_argument);
}

TEST(Mfrac, AutoReductionAfterCanonicalConstruction) {
    auto ctx = make_ctx({"x"});
    // 2*x / (4*x) should reduce to 1/2 — canonical form has denominator
    // monic-ish so equality with explicitly-canonical form works.
    auto num = alg::Mpoly::from_string(ctx, "2*x");
    auto den = alg::Mpoly::from_string(ctx, "4*x");
    alg::Mfrac r(std::move(num), std::move(den));

    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);
    auto half_frac = alg::Mfrac::from_fmpq(ctx, half);
    EXPECT_EQ(r, half_frac);
    fmpq_clear(half);
}

TEST(Mfrac, AdditionAutoReduces) {
    auto ctx = make_ctx({"x"});
    // 1/(x+1) + 1/(x-1) = ((x-1) + (x+1)) / ((x+1)(x-1)) = 2x/(x^2 - 1)
    auto a = alg::Mfrac(alg::Mpoly::one(ctx),
                         alg::Mpoly::from_string(ctx, "x + 1"));
    auto b = alg::Mfrac(alg::Mpoly::one(ctx),
                         alg::Mpoly::from_string(ctx, "x - 1"));
    auto sum = a + b;
    auto expected = alg::Mfrac(alg::Mpoly::from_string(ctx, "2*x"),
                                 alg::Mpoly::from_string(ctx, "x^2 - 1"));
    EXPECT_EQ(sum, expected);
}

TEST(Mfrac, MultiplicationCancels) {
    auto ctx = make_ctx({"x"});
    // (x+1)/x * x/(x-1) = (x+1)/(x-1)
    auto a = alg::Mfrac(alg::Mpoly::from_string(ctx, "x + 1"),
                         alg::Mpoly::variable(ctx, 0));
    auto b = alg::Mfrac(alg::Mpoly::variable(ctx, 0),
                         alg::Mpoly::from_string(ctx, "x - 1"));
    auto prod = a * b;
    auto expected = alg::Mfrac(alg::Mpoly::from_string(ctx, "x + 1"),
                                 alg::Mpoly::from_string(ctx, "x - 1"));
    EXPECT_EQ(prod, expected);
}

TEST(Mfrac, DivideByZeroThrows) {
    auto ctx = make_ctx({"x"});
    auto a = alg::Mfrac::one(ctx);
    auto b = alg::Mfrac::zero(ctx);
    EXPECT_THROW(a / b, std::domain_error);
}

TEST(Mfrac, IsPolynomial) {
    auto ctx = make_ctx({"x"});
    auto poly = alg::Mfrac::from_mpoly(alg::Mpoly::from_string(ctx, "x^2 + 1"));
    EXPECT_TRUE(poly.is_polynomial());

    auto rat = alg::Mfrac(alg::Mpoly::one(ctx),
                           alg::Mpoly::from_string(ctx, "x"));
    EXPECT_FALSE(rat.is_polynomial());
}

TEST(Mfrac, DerivativeQuotientRule) {
    auto ctx = make_ctx({"x"});
    // d/dx [ x / (x^2 + 1) ] = (1*(x^2+1) - x*2x) / (x^2+1)^2 = (1-x^2)/(x^2+1)^2
    auto r = alg::Mfrac(alg::Mpoly::variable(ctx, 0),
                          alg::Mpoly::from_string(ctx, "x^2 + 1"));
    auto d = r.derivative(0);

    auto expected = alg::Mfrac(alg::Mpoly::from_string(ctx, "1 - x^2"),
                                 alg::Mpoly::from_string(ctx, "x^4 + 2*x^2 + 1"));
    EXPECT_EQ(d, expected);
}

TEST(Mfrac, SubstituteRationalIntoVariable) {
    auto ctx = make_ctx({"x"});
    auto r = alg::Mfrac(alg::Mpoly::from_string(ctx, "x + 1"),
                          alg::Mpoly::from_string(ctx, "x - 1"));
    fmpq_t two;
    fmpq_init(two);
    fmpq_set_si(two, 2, 1);
    auto sub = r.substitute(0, two);   // (2+1)/(2-1) = 3
    auto expected = alg::Mfrac::from_si(ctx, 3);
    EXPECT_EQ(sub, expected);
    fmpq_clear(two);
}

TEST(Mfrac, EvaluateAtRationals) {
    auto ctx = make_ctx({"x", "y"});
    // r = (x + y) / (x - y); at x=3, y=1: 4/2 = 2
    auto r = alg::Mfrac(alg::Mpoly::from_string(ctx, "x + y"),
                          alg::Mpoly::from_string(ctx, "x - y"));

    fmpq_t three, one;
    fmpq_init(three); fmpq_init(one);
    fmpq_set_si(three, 3, 1);
    fmpq_set_si(one, 1, 1);
    fmpq* values[] = {three, one};
    fmpq_t out;
    fmpq_init(out);
    EXPECT_TRUE(r.evaluate(values, out));

    fmpq_t expected;
    fmpq_init(expected);
    fmpq_set_si(expected, 2, 1);
    EXPECT_NE(fmpq_equal(out, expected), 0);

    fmpq_clear(three); fmpq_clear(one);
    fmpq_clear(out);
    fmpq_clear(expected);
}

TEST(Mfrac, EvaluateAcbWithAcbInputs) {
    auto ctx = make_ctx({"x"});
    auto r = alg::Mfrac(alg::Mpoly::from_string(ctx, "x + 1"),
                          alg::Mpoly::from_string(ctx, "x"));
    amflow::numeric::AcbVector vals(1);
    acb_set_si(vals.at(0), 4);

    amflow::numeric::AcbValue out;
    r.evaluate_acb(out.raw(), vals.raw());

    fmpq_t five_quarters;
    fmpq_init(five_quarters);
    fmpq_set_si(five_quarters, 5, 4);
    amflow::numeric::AcbValue expected;
    expected.set_fmpq(five_quarters);
    EXPECT_NE(acb_overlaps(out.raw(), expected.raw()), 0);
    fmpq_clear(five_quarters);
}

TEST(Mfrac, MoveSemantics) {
    auto ctx = make_ctx({"x"});
    auto a = alg::Mfrac::from_si(ctx, 7);
    auto b = std::move(a);
    auto expected = alg::Mfrac::from_si(ctx, 7);
    EXPECT_EQ(b, expected);
}

TEST(Mfrac, CloneIsDeepCopy) {
    auto ctx = make_ctx({"x"});
    auto a = alg::Mfrac::from_si(ctx, 5);
    auto b = a.clone();
    a *= alg::Mfrac::from_si(ctx, 2);   // a = 10
    auto five = alg::Mfrac::from_si(ctx, 5);
    auto ten  = alg::Mfrac::from_si(ctx, 10);
    EXPECT_EQ(b, five);
    EXPECT_EQ(a, ten);
}
