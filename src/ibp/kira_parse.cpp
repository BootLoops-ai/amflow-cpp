// SPDX-License-Identifier: MIT
// ibp::kira_parse_expression — Mathematica-style algebraic parser.
//

#include "amflow/ibp/kira.hpp"

#include <cctype>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>

#include <flint/fmpq.h>
#include <flint/fmpz.h>

namespace amflow::ibp {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

namespace {

class Parser {
public:
    Parser(const std::shared_ptr<MpolyContext>& ctx, const std::string& src)
        : ctx_(ctx), src_(src), pos_(0) {}

    Mfrac parse() {
        skip_ws();
        Mfrac r = parse_add();
        skip_ws();
        if (pos_ != src_.size()) {
            throw std::runtime_error(
                "kira_parse_expression: trailing input at position "
                + std::to_string(pos_) + " (\"" + src_.substr(pos_) + "\")");
        }
        return r;
    }

private:
    void skip_ws() {
        while (pos_ < src_.size()
                && std::isspace((unsigned char)src_[pos_])) ++pos_;
    }

    Mfrac parse_add() {
        Mfrac left = parse_mul();
        while (true) {
            skip_ws();
            if (pos_ < src_.size() && src_[pos_] == '+') {
                ++pos_;
                Mfrac r = parse_mul();
                left += r;
            } else if (pos_ < src_.size() && src_[pos_] == '-') {
                ++pos_;
                Mfrac r = parse_mul();
                left -= r;
            } else break;
        }
        return left;
    }

    Mfrac parse_mul() {
        Mfrac left = parse_unary();
        while (true) {
            skip_ws();
            if (pos_ < src_.size() && src_[pos_] == '*') {
                ++pos_;
                Mfrac r = parse_unary();
                left *= r;
            } else if (pos_ < src_.size() && src_[pos_] == '/') {
                ++pos_;
                Mfrac r = parse_unary();
                left /= r;
            } else break;
        }
        return left;
    }

    Mfrac parse_unary() {
        skip_ws();
        if (pos_ < src_.size() && src_[pos_] == '+') { ++pos_; return parse_unary(); }
        if (pos_ < src_.size() && src_[pos_] == '-') {
            ++pos_;
            Mfrac r = parse_unary();
            return -r;
        }
        return parse_pow();
    }

    Mfrac parse_pow() {
        Mfrac base = parse_atom();
        skip_ws();
        if (pos_ < src_.size() && src_[pos_] == '^') {
            ++pos_;
            skip_ws();
            long sign = 1;
            if (pos_ < src_.size() && src_[pos_] == '-') { sign = -1; ++pos_; }
            else if (pos_ < src_.size() && src_[pos_] == '+') { ++pos_; }
            std::string num;
            while (pos_ < src_.size() &&
                    std::isdigit((unsigned char)src_[pos_])) {
                num += src_[pos_]; ++pos_;
            }
            if (num.empty()) {
                throw std::runtime_error(
                    "kira_parse_expression: '^' must be followed by an integer");
            }
            long e = sign * std::stol(num);
            if (e == 0) {
                return Mfrac::one(ctx_);
            }
            Mfrac result = Mfrac::one(ctx_);
            Mfrac b = (e > 0) ? base.clone() : (Mfrac::one(ctx_) / base);
            long n = std::abs(e);
            for (long i = 0; i < n; ++i) result *= b;
            return result;
        }
        return base;
    }

    Mfrac parse_atom() {
        skip_ws();
        if (pos_ >= src_.size()) {
            throw std::runtime_error(
                "kira_parse_expression: unexpected end of input");
        }
        char c = src_[pos_];
        if (c == '(') {
            ++pos_;
            Mfrac r = parse_add();
            skip_ws();
            if (pos_ >= src_.size() || src_[pos_] != ')') {
                throw std::runtime_error(
                    "kira_parse_expression: missing ')'");
            }
            ++pos_;
            return r;
        }
        if (std::isdigit((unsigned char)c)) {
            std::string num;
            while (pos_ < src_.size() &&
                    std::isdigit((unsigned char)src_[pos_])) {
                num += src_[pos_]; ++pos_;
            }
            fmpz_t big;
            fmpz_init(big);
            int rc = fmpz_set_str(big, num.c_str(), 10);
            if (rc != 0) {
                fmpz_clear(big);
                throw std::runtime_error(
                    "kira_parse: cannot parse integer '" + num + "'");
            }
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_fmpz(q, big);
            Mfrac result = Mfrac::from_fmpq(ctx_, q);
            fmpq_clear(q);
            fmpz_clear(big);
            return result;
        }
        if (std::isalpha((unsigned char)c) || c == '_') {
            std::string ident;
            while (pos_ < src_.size() &&
                    (std::isalnum((unsigned char)src_[pos_])
                        || src_[pos_] == '_')) {
                ident += src_[pos_]; ++pos_;
            }
            long var = ctx_->var_index(ident);
            if (var < 0) {
                throw std::runtime_error(
                    "kira_parse_expression: unknown identifier '" + ident + "'");
            }
            return Mfrac::from_mpoly(Mpoly::variable(ctx_, var));
        }
        throw std::runtime_error(
            "kira_parse_expression: unexpected character '"
            + std::string(1, c) + "' at position "
            + std::to_string(pos_));
    }

    std::shared_ptr<MpolyContext> ctx_;
    const std::string& src_;
    std::size_t pos_;
};

}  // namespace

Mfrac kira_parse_expression(const std::shared_ptr<MpolyContext>& ctx,
                              const std::string& expr) {
    Parser p(ctx, expr);
    return p.parse();
}

}  // namespace amflow::ibp
