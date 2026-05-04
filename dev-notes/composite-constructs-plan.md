# brms SEM Extension — Implementation Plan (FC-SEM)

**Branch:** `claude/brms-composite-constructs-7hU8W`
**Working dir:** `/home/user/brms`
**Reference paper:** Schamberger, Schuberth, Henseler & Rosseel (2026),
*Factor- and Composite-Based Structural Equation Modeling*, arXiv 2508.06112v2,
copy at `refs/2508.06112v2.pdf` (extracted text at `refs/2508.06112v2.txt`).

**Long-term goal:** Evaluate whether **brms** can subsume **lavaan** (frequentist
SEM) and **blavaan** (Bayesian SEM). If yes, ship upstream PRs to
`paul-buerkner/brms` covering at least the FC-SEM subset.

---

## 0. Resource access

- ✅ brms source: local at `/home/user/brms`.
- ✅ lavaan source on GitHub. The `<~` operator is real
  (`R/lav_syntax_parser.R`, `real_operators <- c("=~", "<~", ...)`).
- ✅ blavaan source on GitHub (Stan in `inst/stan/`, entry `bsem`/`bcfa`).
- ✅ Paper PDF/text in `refs/`.

### 0.1 Paper findings that drive the design

The paper presents **FC-SEM** (factor- and composite-based SEM). Concretely:

1. **Two indicator equations** (Eqs. 1–2):
   `y_f = Λ_f η_f + ε_f` (reflective / common-factor block)
   `η_c = W' y_c` (composite block — *no* error term).
2. **Composite-indicator covariances are free.** The block `T` (Eq. 3)
   contains the unrestricted `cov(y_j, y_k)` for indicators forming the
   *same* composite (off-block elements are zero). These are model
   parameters.
3. **Composite loadings** are derived, not free
   (Eq. 7): `Λ_c = T W (W' T W)^{-1}`. They are the OLS regression
   coefficients of composites on their indicators.
4. **Model-implied covariance** (Eq. 13):
   `Σ(θ) = Λ V(η) Λ' + Θ`, with
   `Λ = blockdiag(Λ_f, Λ_c)`, `Θ = blockdiag(Θ_f, Θ_c)`,
   `Θ_c = T − Λ_c W' T W Λ_c'`.
5. **Structural-error variances for composite-η are not free.** The
   diagonal of `V(η)` for composites must equal the diagonal of `W' T W`
   (Eq. 10 + the constraint just below it).
6. **Estimator.** ML under `y ∼ N(0, Σ(θ))` — the standard
   Wishart/sufficient-statistic fit function (Eq. 15):
   `log L = −(n/2)[log|Σ(θ)| + tr(S Σ(θ)^{−1})] − (NP/2) log(2π)`,
   minimised in the discrepancy form
   `F = log|Σ(θ)| + tr(S Σ(θ)^{−1}) − log|S| − P` (Eq. 16).
7. **Default scaling.** Fix `w_1 = 1` per composite (matches lavaan's
   default `=~` scaling). Alternatives: fix composite variance, fix
   weight sum.
8. **Identification.** Each composite must be related to at least one
   variable that is *not* one of its own indicators (Dijkstra 2017).
9. **lavaan's current FC-SEM** (≥ 0.6-20) requires
   `optim.gradient = "numerical"` and does not yet support mean
   structure or scaling-to-unit-variance during estimation. Footnote 4:
   the `<~` operator in earlier lavaan versions had different semantics
   ("formative latent with disturbance variance fixed to zero"); we
   target the new (≥ 0.6-20) FC-SEM semantics.
10. **Reference R code.** Authors host an OSF replication package at
    `https://osf.io/9tm2y/` — usable as a parity benchmark.

### 0.2 Implications of the paper for the plan

FC-SEM is composites + reflective + structural model, fitted by ML on a
multivariate-normal likelihood through the model-implied Σ(θ). The user
has, however, asked that the brms-side implementation reuse the existing
`mi()` / `me()` mechanism (`R/formula-sp.R`) wherever possible. That
constraint changes the picture substantially:

- **Reflective `=~` adds no new term class.** brms already supports
  latent-as-fully-missing via the standard idiom
  `bf(eta | mi() ~ 1) + bf(y_k ~ mi(eta), ...)` (see
  `?brms::mi` and the `brms_multivariate`/missing-data vignettes). A
  string-level desugarer is enough: `eta =~ y1 + y2 + y3` becomes the
  multivariate `bf(...) + bf(...) + bf(...) + bf(eta | mi() ~ 1)` with
  the first loading fixed via `set_prior("constant(1)", ...)` for
  identification. **No Stan additions, no new R class.**
- **Composites need a new constructor**, but shaped exactly like `mi()`:
  a tiny `composite()` helper in `R/formula-sp.R` (or a new
  `R/formula-cc.R`) that captures unevaluated arguments, returns a
  `cc_term`/`sp_term` object, and is parsed by `brmsterms()` in the
  same dispatch already used by `mi_term`/`me_term`/`mo_term`. `<~` is
  *only* sugar that emits a `composite()` call.
- **MVN likelihood is independent of the term types.** Without it, brms
  with `composite()` and `mi()` already gives a valid Bayesian posterior
  that *replaces blavaan* (whose Stan code also evaluates per-row MVN).
  *Replacing lavaan's ML* is what needs the sufficient-statistic
  Wishart form. We therefore separate PR-B (composites) from PR-C (MVN
  likelihood); the latter is only required for exact lavaan ML parity.
- **`algorithm = "optimize"` (PR-A) is still first** — smallest patch,
  and it unblocks lavaan parity in the harness.

### 0.3 Upstream context — issue #304 and the `brms3` branch

This subsection rebaselines the plan against ongoing upstream work
discovered after the first draft.

**Resource access (this turn):**
- ✅ `upstream/brms3` and `upstream/generalize-mi` fetched locally
  (added as remote `upstream`; SHAs `4d3b308a` / `b9184f1c`).
  `brms3` is `Version: 2.99.9002` — the brms-3.0 development branch.
- ❌ **GitHub issue
  [paul-buerkner/brms#304](https://github.com/paul-buerkner/brms/issues/304)**
  — both WebFetch and the unauthenticated GitHub API hit rate limits
  in this turn. The thread is referenced by every external SEM-in-brms
  discussion but I have *not* read its current contents in this
  revision. **Action for the user:** paste the issue body and the
  maintainer's most recent comments (or the relevant excerpts) under
  `dev-notes/refs/issue-304/` so the next pass can ground the plan in
  Paul's stated direction. The plan below is conservative w.r.t.
  whatever Paul has already proposed in #304 — if he has a concrete
  syntax in mind, ours should yield to his.

**What `brms3` already ships that is directly relevant:**

1. **PR #1733 (`generalize-mi`, merged into `brms3`).** `mi(x, idx)`
   now accepts **non-unique** indices, and the `NEWS.md` entry says
   verbatim:
   > "Extend `mi` addition terms to handle non-unique indexes via
   > argument `idx`. This allows to express that multiple observations
   > share the same latent missing value."
   This is *exactly* the upstream-blessed primitive for reflective
   measurement: a single latent value `eta_i` per person, shared
   across that person's indicator observations, in a long-format
   multivariate response. **The reflective desugarer in our PR-B does
   not need to invent any latent-sharing mechanism — it desugars to
   `mi(eta, idx = person_id)` plus `bf(eta | mi(idx = person_id) ~ 1)`
   verbatim**, and brms3 handles the rest.
2. **PR #1687, the new `re()` predictor.** brms3 exports a new
   special term that lets *group-level effects defined in one part of
   the model be used as predictors in another part*. Stan-side wiring
   already exists (`stan_re()` + `Jsub_<id>` index data, see
   `R/stan-predictor.R` in brms3). For our purposes:
   - `re()` is **another viable backbone** for `=~`-style reflective
     models — the latent value can be a group-level intercept that
     `re()` then supplies as a predictor. We will **not** rely on
     `re()` for PR-B (the `mi(idx)` path is more direct for single
     latent variables) but we will document it in the vignette as the
     idiomatic way to express *latent-variable-as-predictor* once a
     latent has been declared.
3. **No upstream composite (`<~`) work.** A grep over `upstream/brms3`
   for `composite|formative|<~|FC-SEM` returns nothing. `composite()`
   remains net-new in PR-B.
4. **No upstream `algorithm = "optimize"`.** brms3
   `R/backends.R` `algorithm_choices()` is unchanged
   (`sampling/meanfield/fullrank/pathfinder/laplace/fixed_param`).
   PR-A remains net-new.
5. **No upstream sufficient-statistic MVN likelihood.** PR-C remains
   net-new.

**What is now redundant in our plan:**

- Any plumbing for "latent value shared across multiple observations"
  — fully covered by `mi(idx)` in brms3. Drop any language in §3.2.1
  suggesting we add such plumbing.
- Manual emission of `mi()` `Yl_*` predictors and the latent
  imputation block — also handled by brms3.
- The vague "use brms's existing latent-as-fully-missing idiom"
  language in §3.2.1 is now specific: *"desugar `=~` to a long-format
  multivariate `bf(...)` set that uses `mi(eta, idx = person_id)` and
  `bf(eta | mi(idx = person_id) ~ 1)`, exactly as in the brms3 `?mi`
  example with non-unique `idx`."*

**Implications for upstream targeting and merge strategy:**

- **Target `upstream/brms3`, not `upstream/master`.** The reflective
  sugar in PR-B depends on the generalized `mi(idx)` shipped only in
  brms3. Submitting PR-B against current master would force us to
  back-port #1733 ourselves, which Paul will (rightly) reject.
- **PR ordering is unchanged but rebased:** PR-A → PR-B → PR-C, all
  branched from `upstream/brms3`.
- **Snapshot baseline must be regenerated against brms3.** Stan code
  emitted by brms3 differs from master in `formula-sp.R`,
  `formula-ad.R`, `prepare_predictions.R`, and `stan-predictor.R`.
  CP-0 below should snapshot brms3, not master.
- **CRAN/release timing.** brms3 will ship as brms 3.0; if Paul is
  close to releasing, our PRs land in 3.0 or 3.1. If brms3 is months
  away, our PRs effectively delay until then. Ask on issue #304
  before any PR is opened.
- **Working branch.** Recommend creating a working branch *off*
  `upstream/brms3` rather than continuing on top of `upstream/master`:
  ```
  git checkout -b sem-on-brms3 upstream/brms3
  git cherry-pick <plan-commit>          # bring just the plan over
  ```
  All future PR-A/B/C commits land on `sem-on-brms3`. The
  `claude/brms-composite-constructs-7hU8W` branch is preserved as the
  planning branch.

**Open question for the user:** if issue #304 contains a concrete
syntax proposal from Paul (e.g. an `fa()` or `latent()` constructor
he prefers), our `composite()` constructor should adopt that
naming/shape rather than introducing a parallel one. The next plan
revision should fold in those preferences.

---

## 1. Adversarial-collaboration team layout

| Role | Responsibility | Performed by |
|---|---|---|
| **Architect** | Translate FC-SEM to brms; choose minimum-touch insertion points; own this file. | Claude |
| **Implementer** | R + Stan code; brms style. | Claude |
| **Adversary** | Before each checkpoint advances, write a 3-bullet "attack list" answering *Why does this not ship upstream?* | Claude (separate prompt or `/review` skill) |
| **Executor** | Runs all R verification scripts (this environment cannot run R). | User |

No checkpoint advances on a red verification or an unanswered Adversary
bullet.

---

## 2. Constraints

1. **Minimum diff.** New behaviour in new files
   (`R/formula-cc.R`, `R/stan-cc.R`, `R/lav-string.R`,
   `R/likelihood-mvn.R`, `R/optimize-backend.R`). Existing files
   touched only at single, additive dispatch points behind feature
   switches. Reflective measurement (`=~`) reuses the existing `mi()`
   machinery and adds **no new R class and no new Stan code**.
1a. **Target branch.** All PRs branch from `upstream/brms3`
   (Version 2.99.9002), not `upstream/master`. PR-B's reflective
   sugar depends on the generalized `mi(idx)` from PR #1733 which is
   merged only into brms3. See §0.3.
2. **Backwards compatibility.** All existing brms3 tests pass unmodified.
   New behaviour activates only when the user opts in
   (`composite()`, `lav("...")`, `algorithm = "optimize"`,
   `set_likelihood("mvn")`).
3. **Style.** 2-space indent; snake_case; one new S3 class only
   (`cc_term`/`sp_term`, the same dual-class pattern `mi_term` uses);
   roxygen2 docs; no new hard `Imports:`. Suggests-only deps
   (`lavaan`, `blavaan`, `MASS`) for the parity harness.
4. **No code vendored from lavaan/blavaan.** Re-implement under brms's
   GPL-2 with credit in `NEWS.md`.
5. **Three upstream PRs**, in this order:
   - **PR-A: `algorithm = "optimize"`** — pure backend addition.
   - **PR-B: `composite()` helper + `<~`/`=~` lavaan-string sugar**
     — modelled on `mi()` / `me()` in `R/formula-sp.R`. Reflective `=~`
     desugars into existing `mi()`-based brms idioms (no new
     parameters); composite `<~` desugars into a new `composite()`
     constructor that adds weight parameters and a transformed-parameter
     latent. **No new likelihood mode required.**
   - **PR-C: `set_likelihood("mvn")`** — optional sufficient-statistic
     ML for `mvbf(...) + set_rescor(TRUE)` models. Required only for
     *exact* parity with lavaan's ML fit function on FC-SEM models;
     does not change posteriors under MCMC.

---

## 3. Technical design

### 3.1 `algorithm = "optimize"` (PR-A)

Add `"optimize"` to `algorithm_choices()` in `R/backends.R`. New
private functions `.fit_model_rstan_optimize()` /
`.fit_model_cmdstanr_optimize()` call `rstan::optimizing()` /
`cmdstanr_model$optimize()` and return a thin subclass
`brmsfit_optimize` of `brmsfit`.

- Methods: `summary`, `print`, `fixef`, `coef`, `vcov` (from Hessian).
- `posterior_*` methods error with "no draws under
  `algorithm = 'optimize'`".
- This is **MAP** unless priors are flat; document loudly.
  `prior = empty_prior()` gives true ML. Alternatively expose a switch
  `optimize(... , flat_prior = TRUE)` that internally substitutes
  `target += 0` for all `*_lpdf`/`*_lprior` priors. Default off.

### 3.2 `composite()` helper + `<~`/`=~` sugar (PR-B)

This is the headline feature. It must look and feel exactly like
`mi()`/`me()` so that brms users discover it the same way.

#### 3.2.1 Reflective measurement (`=~`) reuses `mi()`

A reflective measurement model is, in brms3's vocabulary, a
**latent-as-fully-missing** variable shared across an indicator's
observations via the new non-unique `mi(idx)` mechanism (PR #1733,
shipped in `brms3`; see §0.3). For `eta =~ y1 + y2 + y3` over `N`
people, the canonical desugaring is:

```r
# wide -> long: stack y1, y2, y3 into one column with an
# 'indicator' factor and a 'person_id' grouping column.
# eta is N rows long; the long-format response is 3*N rows long.

bform <- bf(y_long | mi() ~ 0 + indicator +
                              indicator:mi(eta, idx = person_id)) +
         bf(eta    | mi(idx = person_id) ~ 1) +
         set_rescor(FALSE)

# pin first loading at 1 for identification:
prior <- set_prior("constant(1)", coef = "indicatory1:mieta",
                   resp = "ylong")
```

Equivalent wide-form (kept as an option for users who prefer it; brms3
also accepts this and the `idx` is implicit because each row is one
person):

```r
bform <- bf(y1 ~ 0 + mi(eta)) +
         bf(y2 ~ 0 + mi(eta)) +
         bf(y3 ~ 0 + mi(eta)) +
         bf(eta | mi() ~ 1) +
         set_rescor(FALSE)
```

Therefore **the reflective operator requires no new R class and no new
Stan code in brms3**. PR-B contributes only:

- A small parser in `R/lav-string.R` that recognises `=~` inside
  `lav("...")` strings and emits the multi-`bf(...)` call above.
- A small helper `as_brms_reflective(latent, indicators, ...)`
  (internal, not exported) that returns the same multi-formula object
  directly without going through the string parser. Used by
  `lav_to_brms()` and by users who prefer programmatic construction.
  Documented as *"sugar over `mi()`"*; the docstring points at `?mi`.
  No new exported user-facing function for reflective.
- A vignette section showing the equivalence.

#### 3.2.2 Composite measurement (`<~`) — new `composite()` helper

Composites cannot be expressed as `mi()` because they are deterministic
functions of *observed* (not latent) indicators with free
indicator-indicator covariances. A new `sp_term`-style helper is needed,
modelled byte-for-byte on `mi()`:

```r
# R/formula-sp.R (or new R/formula-cc.R if Paul prefers a separate file)
#'
#' Composite (Formative) Constructs in \pkg{brms} Models
#'
#' Specify a composite latent variable formed as a (weighted) linear
#' combination of observed indicators. The function does not evaluate
#' its arguments — it exists purely to help set up a model, in the
#' same way as \code{\link{mi}} and \code{\link{me}}.
#'
#' @param formula A two-sided formula whose LHS is the latent name and
#'   whose RHS is a sum of indicators, e.g. \code{eta ~ y1 + y2 + y3}.
#' @param identification One of \code{"first"} (default; pin
#'   \code{w_1 = 1}), \code{"var"} (pin \code{Var(eta) = 1}), or
#'   \code{"sum"} (constrain \code{sum(w) = 1}).
#' @param weights Optional named numeric vector of fixed weights.
#' @param fixed_T Logical. If \code{TRUE} (default, matches lavaan
#'   0.6-20 FC-SEM), composite-indicator (co)variances are fixed at
#'   their sample values; if \code{FALSE}, they become free parameters
#'   gathered in \code{T}.
#' @return An object of class \code{c("cc_term", "sp_term")} that can
#'   be added to a \code{brmsformula} with \code{+}.
#' @export
composite <- function(formula,
                      identification = c("first", "var", "sum"),
                      weights = NULL, fixed_T = TRUE) {
  identification <- match.arg(identification)
  lhs <- deparse0(formula[[2L]])
  rhs <- all.vars(formula[[3L]])
  label <- deparse0(match.call())
  out <- nlist(latent = lhs, indicators = rhs,
               identification, weights, fixed_T, label)
  class(out) <- c("cc_term", "sp_term")
  out
}
```

Wiring (mirrors what `mi_term` does today):

1. `R/brmsformula.R` accepts `composite()` objects in the `+` chain
   exactly the way it accepts `set_rescor()` / `set_mecor()` results
   today. (One new branch in the existing dispatch; ~3 lines.)
2. `R/brmsterms.R` extracts `cc_term`s alongside the existing
   `sp_term` extraction loop and stores them on the parsed object as
   `bterms$composites` (parallels the existing `bterms$sp` lists).
3. `R/stan-cc.R` (new) emits, for each composite, the parameters and
   transformed-parameter blocks shown in §3.3.1 below.
4. `R/stan-predictor.R` resolves predictor names that match a composite
   latent to the corresponding transformed parameter — single new lookup
   in the existing predictor-substitution loop.

The user-facing API ends up matching brms's existing `mi()` look:

```r
brm(
  bf(z ~ x + composite(eta ~ y1 + y2 + y3)),
  data = d
)

# or, lavaan-style sugar:
brm(lav("z ~ x + eta
         eta <~ y1 + y2 + y3"), data = d)
```

#### 3.2.3 What the `composite()` Stan code looks like

For each composite η_c with K indicators `y_c[, 1:K]`:

```stan
parameters {
  vector[K - 1] w_eta_free;       // identification = "first" by default
  // Free intra-composite (co)variance block, only if fixed_T = FALSE.
  // Otherwise T is precomputed in transformed data from sample stats.
}
transformed parameters {
  vector[K] w_eta;
  w_eta[1] = 1;
  w_eta[2:K] = w_eta_free;
  vector[N] eta = Y_eta_indicators * w_eta;     // N x K times K
}
```

Composite-loading recovery `Λ_c = T W (W' T W)^{-1}` (paper Eq. 7) is
only needed when **fitting via the model-implied Σ(θ)** (PR-C). Under
brms's standard per-row likelihood (PR-B alone), `eta` enters the
predictor as a transformed parameter and the rest of the formula
proceeds as usual; this is the path the user asked for.

### 3.3 `set_likelihood("mvn")` (PR-C)

Adds the sufficient-statistic ML form for `mvbf(...) + set_rescor(TRUE)`
models, optionally with a `composite()` term active. Required only for
exact lavaan FC-SEM ML parity; under MCMC the per-row form gives the
same posterior. Today brms's MVN model with `rescor = TRUE` already uses
`multi_normal_cholesky_lpdf` per row. Add the option to evaluate the
**equivalent sufficient-statistic form** when:

- All responses in `mvbf` are Gaussian.
- `set_rescor(TRUE)`.
- The mean structure is constant across rows for *all* responses (no
  predictors), or row-varying means are explicitly handled (see below).
- Data are complete (no missing). Otherwise fall back with one-line
  message.

Stan emitted (sketch):

```stan
// transformed data
matrix[N, P] Y = ...;          // existing
vector[P]   ybar = ...;        // new, computed once
matrix[P, P] S   = (Y - rep_matrix(ybar', N))' * (...) / N;

// model
target += -0.5 * N * (log_determinant_spd(Sigma)
                      + trace_quad_form(inverse_spd(Sigma), S)
                      + (ybar - mu0)' * mdivide_left_spd(Sigma, ybar - mu0)
                      + P * log(2 * pi()));
```

- Use Cholesky factor of `Sigma` (already tracked by brms) and
  `mdivide_left_tri_low` rather than `inverse_spd` for stability.
- For row-varying mean (predictors present) the speedup is lost; in that
  case keep the per-row form regardless.
- Same posterior, same MAP, just faster and matches the form used by
  lavaan/blavaan. This is the "SEM-style likelihood" the user asked for.

### 3.4 Bidirectional translator (PR-B, secondary)

Two pure functions and a checker, all in `R/lav-string.R`:

```r
# lavaan model string  ->  brmsformula object
lav_to_brms(model, data = NULL,
            check = c("strict", "lenient", "none")) -> brmsformula

# brmsformula object   ->  lavaan model string
brms_to_lav(bform,
            check = c("strict", "lenient", "none")) -> character

# pure validity checker; called automatically inside the two above
lav_brms_translatable(x) -> list(ok = logical, reason = character)
```

Translation rules (translatable subset):

| lavaan | brms |
|---|---|
| `eta =~ y1 + y2 + y3` | `bf(y1 ~ 0+mi(eta)) + bf(y2 ~ 0+mi(eta)) + bf(y3 ~ 0+mi(eta)) + bf(eta\|mi() ~ 1)` with `prior = constant(1)` on first loading |
| `eta <~ y1 + y2 + y3` | `composite(eta ~ y1 + y2 + y3)` |
| `y ~ x1 + x2`         | `bf(y ~ x1 + x2)` |
| `y1 ~~ y2`            | `set_rescor(TRUE)` (only at the model level; per-pair residual covariances are not yet expressible in brms — see "untranslatable" below) |

`lav_brms_translatable()` returns `ok = FALSE` with a `reason` for any of:

- multi-group syntax (`group = "..."`) — out of scope.
- equality / inequality constraints (`a == b`, `a > b`) — needs brms
  custom Stan code, not first-class.
- ordinal indicators with thresholds.
- pairwise residual covariances `y1 ~~ y2` other than via global
  `set_rescor`.
- mean-structure intercepts on latents (lavaan FC-SEM does not yet
  support these either, per paper § "Implementation in lavaan").
- exogenous variances/means (`y1 ~~ y1`, `y1 ~ 1` on its own).

`check = "strict"` (default) errors on any untranslatable element;
`"lenient"` warns and drops them; `"none"` translates the recognised
subset silently. Round-trip test:
`brms_to_lav(lav_to_brms(s)) == canonicalised(s)` for every model in
`dev-notes/parity/models/`. The canonicaliser strips whitespace, sorts
RHS terms, and lower-cases operators.

The translator lives entirely in `R/lav-string.R` and is unit-tested in
`tests/testthat/tests.lav-string.R`. It does NOT call lavaan; it is a
self-contained string operation.

### 3.5 Files added vs. touched

**New files** (all carry full unit tests):

- `R/formula-cc.R`          — `composite()` and `cc_term`. Mirrors the
                              shape of `mi()` in `R/formula-sp.R`.
- `R/stan-cc.R`             — composite-block Stan helpers.
- `R/optimize-backend.R`    — `algorithm = "optimize"` plumbing.
- `R/lav-string.R`          — bidirectional translator
                              (`lav_to_brms`, `brms_to_lav`,
                              `lav_brms_translatable`) and `<~`/`=~`
                              recognition. No new term class for
                              reflective; reflective desugars to `mi()`.
- `R/likelihood-mvn.R`      — `set_likelihood("mvn")` data + model
                              code. PR-C only.
- `inst/chunks/fun_sem_implied.stan` — only if helper functions are
                              reused (e.g. composite-loadings
                              recovery for PR-C); otherwise inline.
- `tests/testthat/tests.optimize-backend.R`
- `tests/testthat/tests.formula-cc.R`
- `tests/testthat/tests.lav-string.R`
- `tests/testthat/tests.likelihood-mvn.R`
- `vignettes/brms_sem.Rmd`  — short worked example. Long-form
                              parity write-up lives in `dev-notes/`.
- `dev-notes/parity/`       — comparison harness; in `.Rbuildignore`.

**Touched files** (single additive insertion each):

- `R/formula-sp.R` (or `R/brmsformula.R`) — register `composite()` as
  a recognised special, alongside the existing `mi()`/`me()`/`mo()`
  registration. Three lines.
- `R/brmsterms.R`           — extend the existing `sp_term` extraction
  loop to also collect `cc_term`s onto `bterms$composites`.
- `R/stan-predictor.R`      — when a predictor name matches a
  composite latent, substitute the transformed-parameter symbol
  written by `R/stan-cc.R`. Single new lookup.
- `R/stan-likelihood.R`     — branch on `likelihood = "mvn"` (PR-C).
- `R/stan-data.R`           — emit `S`, `ybar` when MVN mode is on
  (PR-C).
- `R/backends.R`            — append `"optimize"` to
  `algorithm_choices()`.
- `R/brm.R`                 — pipe `algorithm`, `likelihood`
  arguments.
- `NAMESPACE`, `NEWS.md`, `DESCRIPTION` (Suggests-only bump for the
  parity harness: `lavaan`, `blavaan`, `MASS`).

---

## 4. Phased plan with checkpoints

The cadence is: Implementer writes → Executor runs the **CP-N script**
→ pastes output → Adversary review → next phase.

### CP-0 — Baseline and snapshots (≈ 30 min user time)

```r
setwd("/home/user/brms")
library(devtools)

# 0a. Baseline: existing tests pass on the brms3-rebased branch.
#     Confirm the working tree is rebased onto upstream/brms3 first.
stopifnot(grepl("brms3",
  system("git merge-base --is-ancestor upstream/brms3 HEAD && echo OK",
         intern = TRUE),
  fixed = TRUE) ||
  identical(
    system("git rev-parse upstream/brms3", intern = TRUE),
    system("git merge-base HEAD upstream/brms3", intern = TRUE)))
stopifnot(packageVersion("brms") >= "2.99.9002")
res <- devtools::test(reporter = "summary")
stopifnot(all(as.data.frame(res)$failed == 0))

# 0b. Stan-code snapshots that all later checkpoints must preserve byte-
#     for-byte for non-SEM models.
library(brms)
dir.create("dev-notes/parity/snapshots", recursive = TRUE,
           showWarnings = FALSE)

snap_uni <- brms::make_stancode(
  bf(mpg ~ wt + hp), data = mtcars, family = gaussian())
writeLines(snap_uni, "dev-notes/parity/snapshots/uni_gaussian.stan")

snap_mv <- brms::make_stancode(
  mvbf(mpg ~ wt, hp ~ wt) + set_rescor(TRUE),
  data = mtcars, family = gaussian())
writeLines(snap_mv, "dev-notes/parity/snapshots/mv_gaussian_rescor.stan")

# 0c. Confirm parity-harness deps. lavaan must be ≥ 0.6-20 for FC-SEM.
have_lavaan  <- requireNamespace("lavaan", quietly = TRUE)
stopifnot(have_lavaan,
          utils::packageVersion("lavaan") >= "0.6.20")
have_blavaan <- requireNamespace("blavaan", quietly = TRUE)
cat("lavaan:",  as.character(utils::packageVersion("lavaan")),
    " blavaan:", if (have_blavaan)
      as.character(utils::packageVersion("blavaan")) else "MISSING",
    "\n")
```

**Pass:** all blocks green; lavaan ≥ 0.6-20. Install blavaan before CP-4.

---

### CP-1 — `algorithm = "optimize"` (PR-A scope)

Implementer delivers §3.1.

```r
setwd("/home/user/brms"); devtools::load_all()
library(testthat); library(brms)

# 1a. New tests pass.
test_file("tests/testthat/tests.optimize-backend.R")

# 1b. No regressions.
test_dir("tests/testthat",
         filter = "^(?!tests\\.optimize-backend).*",
         perl = TRUE, reporter = "summary")

# 1c. Snapshot byte-equality.
old_uni <- readLines("dev-notes/parity/snapshots/uni_gaussian.stan")
new_uni <- strsplit(brms::make_stancode(
   bf(mpg ~ wt + hp), data = mtcars, family = gaussian()), "\n")[[1]]
stopifnot(identical(new_uni, old_uni))

# 1d. ML parity with lavaan on a saturated regression.
set.seed(1); n <- 500
d <- data.frame(x = rnorm(n)); d$y <- 0.7 * d$x + rnorm(n, 0, 0.5)

lav_fit <- lavaan::sem("y ~ x", data = d, meanstructure = TRUE)
brm_fit <- brm(y ~ x, data = d, family = gaussian(),
               algorithm = "optimize",
               prior = brms::empty_prior(),
               refresh = 0)

lav_b <- lavaan::coef(lav_fit)["y~x"]
brm_b <- fixef(brm_fit)["x", "Estimate"]
stopifnot(abs(lav_b - brm_b) < 1e-3)
cat(sprintf("lavaan b = %.6f, brms-ML b = %.6f\n", lav_b, brm_b))
```

**Adversary checklist for CP-1:**
- Does `summary.brmsfit_optimize` invent MCMC diagnostics that don't
  exist?
- What happens with `algorithm = "optimize"` and group-level effects
  (random intercepts)? Define behaviour explicitly — error or warn.
- Hessian-based standard errors: which transformations are applied? Are
  they on the constrained or unconstrained scale?

---

### CP-2 — `composite()` + `<~`/`=~` sugar + translator (PR-B scope)

Implementer delivers §3.2 and §3.4. **No** `set_likelihood("mvn")` yet
— this checkpoint validates that per-row brms with `composite()` and
`mi()`-based reflective sugar matches lavaan and blavaan within
posterior tolerances. (PR-C will tighten ML parity at CP-3.)

```r
setwd("/home/user/brms"); devtools::load_all()
library(testthat); library(brms); library(lavaan)

# 2a. Unit tests.
test_file("tests/testthat/tests.formula-cc.R")
test_file("tests/testthat/tests.lav-string.R")

# 2b. Translator round-trip for every parity model. Strict mode must
#     refuse anything not in the translatable subset; round-trips of
#     supported models must be exact (after canonicalisation).
src_models <- c(
  "y ~ x",
  "eta =~ y1 + y2 + y3",
  "eta <~ y1 + y2 + y3
   z ~ eta",
  "f =~ y1 + y2 + y3
   g =~ y4 + y5 + y6
   c <~ y7 + y8 + y9
   g ~ f + c"
)
for (s in src_models) {
  bform <- brms::lav_to_brms(s, check = "strict")
  back  <- brms::brms_to_lav(bform, check = "strict")
  stopifnot(brms:::canonicalise_lav(back) ==
            brms:::canonicalise_lav(s))
}

# 2c. Strict-mode rejections (untranslatable elements MUST error).
expect_error(brms::lav_to_brms("y1 ~~ y2",      check = "strict"))
expect_error(brms::lav_to_brms("a == b",        check = "strict"))
expect_error(brms::lav_to_brms("y1 + y2 ~ x",   check = "strict"))

# 2d. Reflective parity via mi()-based sugar. Single common factor,
#     three indicators. Compare brms posterior mean to lavaan ML.
set.seed(3); n <- 800
eta_true <- rnorm(n)
d <- data.frame(
  y1 = 1.0 * eta_true + rnorm(n, 0, sqrt(0.5)),
  y2 = 0.8 * eta_true + rnorm(n, 0, sqrt(0.5)),
  y3 = 0.9 * eta_true + rnorm(n, 0, sqrt(0.5))
)
lav_cf <- lavaan::sem("eta =~ y1 + y2 + y3", data = d)
brm_cf <- brm(brms::lav_to_brms("eta =~ y1 + y2 + y3"), data = d,
              chains = 2, iter = 1000, warmup = 500, refresh = 0)
lav_l <- lavaan::coef(lav_cf)[c("eta=~y2", "eta=~y3")]
brm_l <- fixef(brm_cf)[c("y2_mieta", "y3_mieta"), "Estimate"]   # naming TBD
stopifnot(max(abs(lav_l - brm_l) / abs(lav_l)) < 0.05)

# 2e. Composite parity, three indicators + Dijkstra-required outside link.
set.seed(4); n <- 1000
T_chol <- chol(matrix(c(6, 2, 1.2, 2, 5, 1.5, 1.2, 1.5, 2), 3))
y_c <- matrix(rnorm(n * 3), n, 3) %*% T_chol
z   <- 0.4 * (1.0 * y_c[, 1] + 0.4 * y_c[, 2] + 0.6 * y_c[, 3]) +
       rnorm(n, 0, 1)
d2 <- data.frame(y7 = y_c[, 1], y8 = y_c[, 2], y9 = y_c[, 3], z = z)

lav_c <- lavaan::sem("eta <~ y7 + y8 + y9
                      z   ~ eta",
                     data = d2, optim.gradient = "numerical")
brm_c <- brm(
  bf(z ~ eta) + composite(eta ~ y7 + y8 + y9, identification = "first"),
  data = d2,
  chains = 2, iter = 1000, warmup = 500, refresh = 0
)
lav_w <- lavaan::coef(lav_c)[c("eta<~y8", "eta<~y9")]
brm_w <- fixef(brm_c)[c("w_eta_2", "w_eta_3"), "Estimate"]
stopifnot(max(abs(lav_w - brm_w)) < 0.10)

# 2f. blavaan parity for the same composite model (MCMC vs. MCMC).
if (requireNamespace("blavaan", quietly = TRUE)) {
  blav_c <- blavaan::bsem("eta <~ y7 + y8 + y9; z ~ eta",
                          data = d2, n.chains = 2,
                          burnin = 500, sample = 500)
  blav_w <- blavaan::parameterEstimates(blav_c) |>
    subset(lhs == "eta" & op == "<~" & rhs %in% c("y8", "y9"),
           select = est)
  stopifnot(max(abs(blav_w$est - brm_w)) < 0.15)
}

# 2g. Snapshot byte-equality for non-SEM models.
stopifnot(identical(
  readLines("dev-notes/parity/snapshots/uni_gaussian.stan"),
  strsplit(brms::make_stancode(
    bf(mpg ~ wt + hp), data = mtcars, family = gaussian()), "\n")[[1]]))
```

**Adversary checklist for CP-2:**
- Does `composite()` cleanly reuse the `sp_term` machinery, or did we
  introduce a parallel parser? Parallel parsers are an upstream-PR
  killer.
- Reflective sugar via `mi()`: how are user-supplied priors propagated?
  In particular, the `prior = constant(1)` for first-loading
  identification must be auto-injected by the translator; confirm it
  appears in `prior_summary(brm_cf)`.
- Composite-indicator covariances `T`: under brms's per-row likelihood
  these are not modelled at all (the indicators just enter the
  composite as a deterministic linear combination). That is *the
  difference* between brms-PR-B and lavaan FC-SEM. Document explicitly:
  "PR-B gives the same posterior as blavaan; PR-C is required for ML
  parity with lavaan."
- Identification: under `identification = "first"`, sign of η is
  pinned; under `"var"`, multiple optima exist. brms must detect the
  latter and either pin sign or warn.
- Translator: confirm `lav_brms_translatable()` rejects every item in
  the untranslatable list (§3.4) under `check = "strict"`.

---

### CP-3 — `set_likelihood("mvn")` (PR-C scope)

Implementer delivers §3.3. This checkpoint adds the
sufficient-statistic ML form so that we can claim *exact* lavaan ML
parity on FC-SEM models.

```r
setwd("/home/user/brms"); devtools::load_all()
library(testthat); library(brms); library(lavaan)
test_file("tests/testthat/tests.likelihood-mvn.R")

# 3a. Per-row vs. sufficient-statistic forms must give identical MAP
#     for an intercept-only bivariate Gaussian.
set.seed(2); n <- 1000
S <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
Y <- MASS::mvrnorm(n, c(0, 0), S)
d <- data.frame(y1 = Y[, 1], y2 = Y[, 2])

base <- mvbf(y1 ~ 1, y2 ~ 1) + set_rescor(TRUE)
f1 <- brm(base,                         data = d, family = gaussian(),
          algorithm = "optimize", prior = brms::empty_prior(), refresh = 0)
f2 <- brm(base + set_likelihood("mvn"), data = d, family = gaussian(),
          algorithm = "optimize", prior = brms::empty_prior(), refresh = 0)
stopifnot(max(abs(fixef(f1)[, "Estimate"] -
                  fixef(f2)[, "Estimate"])) < 1e-5)

# 3b. lavaan ML parity on the same model.
lav <- lavaan::sem("y1 ~ 1; y2 ~ 1; y1 ~~ y2",
                   data = d, meanstructure = TRUE)
cov_lav <- lavaan::coef(lav)["y1~~y2"]
cov_brm <- VarCorr(f2)$residual__$cov[1, 2]   # field name TBD
stopifnot(abs(cov_lav - cov_brm) / abs(cov_lav) < 1e-2)

# 3c. lavaan ML parity for a one-factor reflective model fitted via
#     PR-B sugar + PR-C likelihood.
lav_cf <- lavaan::sem("eta =~ y1 + y2 + y3", data = d)   # reuse d2 from CP-2 here
brm_cf <- brm(brms::lav_to_brms("eta =~ y1 + y2 + y3") + set_likelihood("mvn"),
              data = d, algorithm = "optimize",
              prior = brms::empty_prior(), refresh = 0)
stopifnot(abs(lavaan::coef(lav_cf)["eta=~y2"] -
              fixef(brm_cf)["y2_mieta", "Estimate"]) < 1e-3)

# 3d. All earlier snapshots still byte-identical.
stopifnot(identical(
  readLines("dev-notes/parity/snapshots/uni_gaussian.stan"),
  strsplit(brms::make_stancode(
    bf(mpg ~ wt + hp), data = mtcars, family = gaussian()), "\n")[[1]]))
```

**Adversary checklist for CP-3:**
- Cholesky stability: any remaining `inverse_spd` in the emitted Stan?
  Replace with `mdivide_left_*` patterns.
- Row-varying mean (predictors present) silently degrades the speedup
  to zero; the emitter must notice and either fall back loudly or skip
  the sufficient-stat form with a printed message.
- `loo()`/`waic()` need per-row log-lik. Either also compute per-row
  log-lik in `generated quantities` (recommended; cost negligible) or
  error on `loo()` under `mvn` mode.
- Composite-indicator covariances `T`: under PR-C with FC-SEM
  semantics, `T` enters Σ(θ) explicitly. Verify the `T`
  block matches Eq. 3 of the paper (only intra-composite covariances
  are non-zero) and that `Λ_c = T W (W' T W)^{-1}` is computed
  per-composite, not jointly across composites.

---

### CP-4 — Paper replication (integration)

Reproduce the paper's empirical example (American Customer Satisfaction
Index, four composites + three reflective + one single-indicator factor).
Data are available via the OSF link in §0.1 and Hwang et al. (2023). If
data are unavailable in this environment, fall back to the population
example in paper §"Scenario Analysis" (paper Eqs. 19–25), which is fully
specified.

```r
setwd("/home/user/brms"); devtools::load_all()
library(brms); library(lavaan)

# 4a. Population scenario from paper §"Scenario Analysis".
#     Two common factors (eta1, eta2), two composites (eta3, eta4).
source("dev-notes/parity/models/05_paper_scenario.R")
# expects this script to (a) simulate from Eqs. 19-25, (b) fit with
# lavaan via FC-SEM, (c) fit with brms FC-SEM, (d) return an
# 'all_close' logical and a side-by-side table.

# 4b. blavaan parity (MCMC), if installed.
if (requireNamespace("blavaan", quietly = TRUE)) {
  source("dev-notes/parity/models/05_paper_scenario_mcmc.R")
}
```

**Pass criteria:**
- 4a: brms FC-SEM ML and lavaan FC-SEM ML agree to **rel. err. < 1 %**
  on every parameter the paper reports in Tables 1–3 of the simulation
  section (loadings, weights, structural betas, σ²(ζ_2), V(η_1)).
- 4b: brms-MCMC posterior mean within **10 %** of blavaan posterior
  mean on the same set, under matched priors and chain length.

**Adversary checklist for CP-4:**
- Are the same scaling defaults being used on both sides? (Paper says
  lavaan default is `w_1 = 1`; brms must match.)
- Does any standardisation step happen implicitly inside lavaan that
  brms is not doing?
- Identifiability of the simulation: does the structural model satisfy
  Dijkstra's rule for both composites?

---

### CP-5 — Documentation, vignette, R CMD check

```r
setwd("/home/user/brms")
devtools::document()
devtools::check(args = "--as-cran")    # 0 errors, 0 warnings, 0 notes
devtools::test()
covr::package_coverage()               # informational
```

**Pass:** R CMD check fully clean; new code ≥ 90 % line coverage;
`vignettes/brms_sem.Rmd` builds and reproduces the paper's
"Sample-Perfectly-Representing-Population" results table.

---

## 5. Parity harness

```
dev-notes/parity/
  snapshots/                   # CP-0 snapshots; never overwrite
  models/
    01_simple_regression.R          # CP-1 ML parity
    02_intercept_only_mvn.R         # CP-2 MVN parity
    03_one_factor.R                 # CP-3 reflective parity
    04_one_composite.R              # CP-3 composite parity
    05_paper_scenario.R             # CP-4 FC-SEM ML parity
    05_paper_scenario_mcmc.R        # CP-4 blavaan parity
  results/                     # CSVs of (model, estimator, parameter,
                               # estimate, se, time) per run
  run_parity.R                 # iterates models/, writes results.csv
```

Each `models/*.R` declares `NAME` and produces a `results` data frame
with one row per (parameter × estimator). The runner can be invoked at
any checkpoint; the Adversary uses it before signing off.

---

## 6. Feasibility verdict (after CP-4)

Architect writes `dev-notes/feasibility-verdict.md` answering:

1. Does brms reproduce **lavaan FC-SEM ML** within 1 % on every model
   in `dev-notes/parity/models/`?
2. Does brms reproduce **blavaan posterior means** within 10 % under
   matched priors and ≥ 1000 post-warmup draws?
3. Is brms runtime within 3× of lavaan (ML) and within 3× of blavaan
   (MCMC)?
4. Lavaan/blavaan features that brms still cannot express even with
   these extensions (multi-group, ordinal indicators with thresholds,
   growth-curve models with phantom variables, MIMIC, mediation with
   bootstrap CIs, FIML for missing data). Each is a follow-up PR.

**Recommendation rule.** If 1–3 are yes and 4 is a finite ranked list,
recommend submitting **PR-A → PR-B → PR-C** to upstream brms. If any
of 1–3 fails, ship the subset that did pass and keep the rest on this
branch.

---

## 7. Risks and open questions

- **Mean structure** is not implemented in lavaan FC-SEM (paper §
  "Implementation in lavaan"). Skipping it in brms keeps parity easy;
  enabling it later is a small additive change to the implied-Σ block
  (add `μ(θ)` and use `multi_normal` instead of mean-zero MVN).
- **Free vs. fixed `T`.** Paper says lavaan defaults to fixing
  composite-indicator covariances to sample values. brms should match
  this default but expose `composite(..., fixed_T = FALSE)` to free
  them. With free `T`, identification gets harder; document.
- **`algorithm = "optimize"` is MAP** unless `prior = empty_prior()`.
  Always print a one-line note in `summary()` explaining which case the
  user is in.
- **MVN sufficient-stat form and `loo()`.** Per-row log-lik must still
  be computed in `generated quantities` for `loo()`/`waic()` to work
  (see Adversary CP-2).
- **Identification under `identification = "var"`** does not pin sign;
  multiple optima exist. brms should detect this and either pin sign or
  warn.
- **The paper's reference R code** lives behind an OSF anonymised link
  — copy results into `dev-notes/parity/expected/` once retrieved so
  the harness can diff against numerical ground truth without
  re-running lavaan every time.
- **Out of scope for this branch but in scope for "replace lavaan +
  blavaan"**: multi-group, ordinal/threshold indicators, FIML for
  missing data, growth curves, bootstrap CIs. Tracked in §6 item 4.

---

## 8. Concrete next action for the user

1. **Paste excerpts of issue
   [paul-buerkner/brms#304](https://github.com/paul-buerkner/brms/issues/304)**
   under `dev-notes/refs/issue-304/`. Both WebFetch and the GitHub API
   were rate-limited in this turn (§0.3). At minimum we need: the
   maintainer's most recent comment on syntax, any explicit naming
   he prefers (`fa`, `latent`, `composite`, ...), and any "this is
   out of scope" statements. The plan will yield to whatever Paul
   has already proposed.
2. **Rebase the working branch onto `upstream/brms3`:**
   ```
   git checkout -b sem-on-brms3 upstream/brms3
   git cherry-pick <plan-commit-sha>
   ```
   Reason: PR-B requires the generalized `mi(idx)` from #1733
   (brms3 only). See §0.3.
3. Run **CP-0 verification** (§4 above) on the rebased branch and
   paste the output back.
4. Confirm whether `vignettes/brms_sem.Rmd` should reproduce the
   paper's full empirical example (American Customer Satisfaction
   Index — needs the OSF data) or only the analytic scenario (paper §
   "Scenario Analysis", self-contained).
5. Confirm three-PR split for upstream (PR-A → PR-B → PR-C, §2 item 5),
   all targeted at `upstream/brms3`.
