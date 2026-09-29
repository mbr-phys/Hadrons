# Flowed stochastic fermion observables in Hadrons

This branch extends Hadrons' gradient-flow machinery in three closely related directions:

1. forward flow of arbitrary collections of `FermionField` and `PropagatorField` objects;
2. stochastic estimators for volume-summed flowed fermion bilinears, including an adjoint-flow formulation; and
3. a factorised stochastic estimator for connected meson two-point functions.

The principal physics application is the measurement and normalization of fermion observables at positive flow time. In particular, the one-point infrastructure supplies scalar and kinetic traces used for flowed scalar bilinears and ringed fermion fields, while the ordinary `Meson` module remains the reference route for deterministic doubly-flowed two-point functions. This note describes the code currently on the `test/AdjointFlow` branch. It does not promote the new estimators to production status: their basic module and smoke-test infrastructure exists, whereas comprehensive numerical validation remains to be done.

## Fermion gradient flow

Let $V_\tau$ denote the gauge field obtained by Wilson flow and let $\Delta_\tau=D_\mu[V_\tau]D_\mu[V_\tau]$ be the covariant Laplacian. The fermion flow used here is

$$
\partial_\tau\chi_\tau=\Delta_\tau\chi_\tau,\qquad \chi_0=\psi,\qquad \chi_\tau=K_\tau\psi .
$$

The corresponding row field carries the adjoint kernel, $\bar\chi_\tau=\bar\psi K_\tau^\dagger$. `MGradientFlow::FermionFlow` implements the existing third-order Runge--Kutta discretisation of this evolution. The branch generalises the module so that a single invocation can evolve both `FermionField` and `PropagatorField` inputs. `propTypes` gives their types individually, or `defaultType` assigns one type to the entire `props` list.

At every step the module advances the gauge field once, constructs the three gauge Runge--Kutta stages, and applies the same stages to every listed field. This is important for stochastic estimators: a noise vector and its Dirac solution can be flowed together without duplicating gauge evolution or allowing their discrete flow trajectories to differ. Results are stored at `meas_interval` (and always at the final step), either with the default `<field>_t<time>` names or a supplied `outProps` list. The final flowed gauge field is available as `<module>_U`. The implementation accepts temporal boundary condition `bc=+1` or `-1` only.

The shared flow action for `FermionFlow` and `WilsonFlow` is `WilsonAction`, implemented as the plaquette-plus-rectangle Grid action with vanishing rectangle coefficient. Thus the gauge stages used to flow fermions are compatible with the Wilson-flow trajectory used by the adjoint route below.

## Positive-flow stochastic traces

For a four-dimensional Dirac matrix $D_f$ at flow time zero, take a noise field satisfying $\mathbb E_\eta[\eta\eta^\dagger]=1$, solve

$$
D_f\phi=\eta,
$$

and forward-flow both members of the pair:

$$
\eta_\tau=K_\tau\eta,\qquad \phi_\tau=K_\tau D_f^{-1}\eta.
$$

For a volume-summed local operator $A_\tau$, built with the gauge field at the same flow time, the Grassmann contraction is

$$
\left\langle\sum_x\bar\chi_\tau(x)A_\tau\chi_\tau(x)\right\rangle_f
=-\operatorname{Tr}\left[K_\tau^\dagger A_\tau K_\tau D_f^{-1}\right]
=-\mathbb E_\eta\left[(\eta_\tau,A_\tau\phi_\tau)\right].
$$

The closed-fermion-loop sign is included by `MContraction::StochasticCondensate{Fermion,Propagator}`. These two module registrations have the same interface and select whether their two input objects are a `FermionField` or a `PropagatorField`. The `gammas` parameter is a space-separated list of Hadrons gamma labels, or `all`; it produces one volume-summed result per insertion. Its scalar result, for example, is

$$
S_f(\tau)=-\sum_x\eta_\tau(x)^\dagger\phi_\tau(x)
+c_{\rm fl}\sum_x\eta_\tau(x)^\dagger\eta_\tau(x).
$$

The $c_{\rm fl}$ counterterm is deliberately applied only to the identity gamma channel; the module records the coefficient actually applied. Its choice is action- and improvement-scheme dependent. In particular, it is not a generic correction for pseudoscalar, vector, or kinetic traces, and the branch makes no claim to provide non-perturbative improvement coefficients.

`MUtilities::DslashField{Fermion,Propagator}` constructs the symmetric gauge-covariant lattice operator

$$
\not\!D\,q(x)=\frac12\sum_\mu\gamma_\mu
\left[V_\mu(\tau,x)q(x+\hat\mu)
-V_\mu(\tau,x-\hat\mu)^\dagger q(x-\hat\mu)\right].
$$

Applying it to $\phi_\tau$ and contracting the output against $\eta_\tau$ through `StochasticCondensate` gives the one-sided kinetic trace $-\langle\eta_\tau,\not\!D[V_\tau]\phi_\tau\rangle$. This is the raw lattice quantity relevant to a ringed-field normalization analysis. The conversion to a chosen $Z_\chi$ convention, including its perturbative factors and normalizations, is intentionally outside Hadrons. Operators requiring the explicit antisymmetrised derivative $\bar\chi\Gamma(D_\mu-\overleftarrow D_\mu)\chi$ require derivatives on both stochastic legs; that full operator is not implemented by the present single `DslashField` path.

Coordinate-diluted sources can be produced by the established `SparseSpinColorDiagonal -> Z2Diluted -> PropagatorVectorUnpack` workflow. For `nSrcs` base noises and sparse factor `nsparse`, this gives $n_{\rm Srcs}n_{\rm sparse}^4$ components. Those components partition each base-noise estimator, so analysis must combine all components and divide by the number of base noises, not additionally by the dilution factor.

## Adjoint flow formulation

The alternative estimator represents the flowed row field by transporting a noise source from the measurement time $\tau$ back to zero flow time. Define

$$
(\partial_s+\Delta_s)\xi(\tau;s)=0,\qquad
\xi(\tau;\tau)=\eta,
$$

so that $\xi(\tau;0)=K_\tau^\dagger\eta$. For the scalar trace,

$$
-\mathbb E_\eta\left[(\xi(\tau;0),D_f^{-1}\xi(\tau;0))\right]
=-\operatorname{Tr}\left[K_\tau^\dagger K_\tau D_f^{-1}\right],
$$

which is the same flowed bilinear as the positive-flow expression. Gauge flow is never reversed: the backward evolution is only the discrete adjoint of the fermion update and must use gauge stages generated on the forward trajectory.

`Evolution::adjoint_laplace_flow` implements the reverse Hermitian adjoint of the forward third-order RK map. For a step from $s$ to $s+\epsilon$, using the forward gauge stages $W_0,W_1,W_2$, its update is

$$
\begin{aligned}
\lambda_3&=\xi_{s+\epsilon},\\
\lambda_2&=\tfrac{3\epsilon}{4}\Delta(W_2)\lambda_3,\\
\lambda_1&=\lambda_3+\tfrac{8\epsilon}{9}\Delta(W_1)\lambda_2,\\
\xi_s&=\lambda_1+\lambda_2+
\tfrac{\epsilon}{4}\Delta(W_0)
\left(\lambda_1-\tfrac89\lambda_2\right).
\end{aligned}
$$

`MGradientFlow::AdjointFermionFlow` applies exactly one such step to a list of `FermionField` and/or `PropagatorField` objects. A full reverse trajectory is therefore represented explicitly in the Hadrons module graph, with modules ordered from the largest flow time to zero. This intentional restriction makes each dependency on an earlier gauge configuration visible and avoids an implicit, potentially memory-heavy checkpointing scheme.

The module either reconstructs $W_0,W_1,W_2$ from its supplied earlier-time gauge field or accepts saved $W_1,W_2$ stages. For the latter route, `WilsonFlow` must be run at fixed step size with `save_history=true`, `save_rk_stages=true`, and `meas_interval=1`. It then publishes $U_0,U_{\epsilon},\ldots$ as `<name>_U_t<time>` and the raw middle stages as `<name>_W1_t<start>` and `<name>_W2_t<start>`. Retaining stages avoids their reconstruction during the reverse evolution but costs two extra gauge fields per step. Gauge history is not exposed for adaptive flow, and both history options default to false.

The central numerical identity is

$$
(\xi(\tau;s),\chi_s)=(\eta,\chi_\tau),
$$

for a forward-flowed test field $\chi$. `MUtilities::NormCheck` was added to record norms, inner products, and a relative field difference for either fermion or propagator fields. `Test_adjoint_flow` constructs the one-step identity data and compares the saved-stage and reconstructed-stage paths in its output. It records diagnostics rather than enforcing a numerical tolerance, and it is not yet a multi-step precision validation.

## Stochastic connected mesons

`MContraction::StochasticMeson` provides a distinct factorised estimator for connected two-point functions. It takes two source--solution pairs, $(\eta_A,\phi_A=D^{-1}\eta_A)$ and $(\eta_B,\phi_B=D^{-1}\eta_B)$, and requires independent noise fields unless `detSrc=1` explicitly declares a deterministic-source use case. With sink and source gamma matrices $\Gamma_{\rm snk}$ and $\Gamma_{\rm src}$, the implemented contraction is arranged as

$$
\left[\sum_{\mathbf y,\,y_0=t_0}
\eta_A(y)^\dagger\Gamma_{\rm src}^\dagger\gamma_5\eta_B(y)\right]
\left[\phi_B(x)^\dagger\gamma_5\Gamma_{\rm snk}\phi_A(x)\right],
$$

followed by the specified `MSink` momentum projection. This placement puts both noise factors at the source wall and both solutions at the sink, with the second quark line oriented through $\gamma_5$ hermiticity. The result stores the absolute source time, gamma pair, noise labels, flow-time label and the correlator over sink times.

The module is agnostic about flow: supplying both pairs at a common positive time produces the intended doubly-flowed stochastic construction; supplying one pair at zero time and one at positive time produces a mixed-flow correlator instead. Its `flowTime` parameter labels output only and does not evolve any field. Independent noise ensembles and a direct comparison to the deterministic `Meson` contraction are essential before using this as a precision correlator measurement. In particular, equality of gamma labels does not make a mixed-flow stochastic result equivalent to an equal-time, doubly-flowed `Meson` result.

## Current scope and validation requirements

The branch adds `Test_propagator_flow`, `Test_fermion_flow`, `Test_adjoint_flow`, and `Test_stochastic_meson` examples to the test build. They exercise the principal data-flow patterns: homogeneous propagator flow, heterogeneous/batched stochastic flow, a one-step adjoint smoke workflow, and the wall-noise meson contraction. They are examples and diagnostics, not a completed production-validation suite.

Before physics production, the following checks are still required:

- establish the adjoint inner-product identity to the expected RK and floating-point precision, including multi-step trajectories and both stage-storage modes;
- test step-size convergence, gauge covariance, temporal boundary-condition handling, and agreement of positive- and adjoint-flow estimators at matched statistics;
- validate stochastic normalization and dilution aggregation, including the $\tau\to0$ limit where applicable;
- validate the stochastic meson estimator against controlled deterministic sources and the conventional `Meson` contraction; and
- perform the physical analysis in a flow-time window $a^2\ll\tau\ll\Lambda_{\rm QCD}^{-2}$, taking continuum and small-flow-time limits with the required matching, improvement, finite-volume, and renormalization analysis.

The positive-flow route is operationally simple and gives all requested flow times after the zero-time inversions, whereas the adjoint route moves the kernel onto the stochastic source and requires forward gauge history. Their relative cost and variance are observable- and ensemble-dependent; the code is intended to make a controlled comparison possible rather than presuming one route is universally preferable.
