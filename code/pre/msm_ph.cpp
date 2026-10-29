#include <TMB.hpp>

// Pre-marriage block with a phase after debut.
//
//   0 V   never had sex, never married
//   1 X1  had sex recently, union formation under way
//   2 X2  had sex, no union forming
//   3 M   in first union (absorbing here; the post-marriage block is separate)
//
// Rates
//   q0(a)  V  -> X1  debut
//   q1(a)  V  -> M   marriage in the year of debut
//   h1(a)  X1 -> M   marriage while the first phase lasts
//   lam    X1 -> X2  constant rate of leaving the first phase
//   h2(a)  X2 -> M   marriage on the age schedule
//
// Why a phase and not a clock. A hazard indexed on duration since debut is read
// off a debut date, and inside a V episode the chain is given the row's reported
// debut date, so the clock forms apply the duration-zero rate to every woman at
// risk. Under a phase, which of the two states after debut a woman occupies is
// what represents the time since her debut, and she enters X1 at her own debut
// whenever the chain takes her there. Among women
// still unmarried the marriage intensity declines from h1 towards h2 with a time
// constant 1/lam.
//
// The likelihood is taken per WOMAN, not per episode, because the phase is not
// observed at an episode boundary. The path is fixed by (afs, afm, aai):
//   type 1  afs = 0, afm = 0      in V to the interview
//   type 2  afs > 0, afm = 0      V, debut year, then X to the interview
//   type 3  afm > 0, afs >= afm   V, then marriage in the year of debut
//   type 4  0 < afs < afm         V, debut year, X, then marriage
// Within a span the forbidden states are zeroed each year, which is the same
// product of interval probabilities the episode form uses when the state after
// debut is a single state.
#define eps CppAD::numeric_limits<double>::epsilon()

template <class Type>
Eigen::SparseMatrix<Type> makeQ_ph(Type q0, Type q1, Type h1, Type lam, Type h2)
{
  Eigen::SparseMatrix<Type> Q(4, 4);
  Q.coeffRef(0, 1) = q0;  Q.coeffRef(0, 3) = q1;
  Q.coeffRef(1, 2) = lam; Q.coeffRef(1, 3) = h1;
  Q.coeffRef(2, 3) = h2;
  Q.coeffRef(0, 0) = -(q0 + q1);
  Q.coeffRef(1, 1) = -(lam + h1);
  Q.coeffRef(2, 2) = -h2;
  Q.coeffRef(3, 3) = Type(0);
  Q.makeCompressed();
  return Q;
}

template <class Type>
Type logSASLH(Type a, Type mu_pk, Type sigma, Type del, Type gam,
              Type p, Type tau)
{
  // the sinh-arcsinh log hazard of code/msm_t.cpp, power form, p held at one
  Type la = log(a / mu_pk);
  Type z  = (exp(p * la) - Type(1)) / (p * sigma);
  Type u  = z + sinh(del);
  Type r  = sqrt(u * u + Type(1));
  Type w  = CppAD::CondExpGt(u, Type(0), u + r, Type(1) / (r - u));
  Type wt = exp(tau * log(w));
  Type sh = Type(0.5) * (wt - Type(1) / wt);
  Type ch = Type(0.5) * (wt + Type(1) / wt);
  Type s  = sh * cosh(tau * del) - ch * sinh(tau * del);
  return -gam * s * s;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type dll = 0;

  DATA_IVECTOR(afs);     // age at first sex, 0 if none
  DATA_IVECTOR(afm);     // age at first marriage, 0 if none
  DATA_IVECTOR(aai);     // age at interview
  DATA_IVECTOR(ic);      // index of this woman's cohort in coh_u, zero based
  DATA_VECTOR(coh_u);    // the distinct birth cohorts, centred, in decades
  DATA_INTEGER(a_max);   // highest age used
  DATA_VECTOR(n);        // weight
  DATA_INTEGER(nopen);   // 1 switches every prior off

  PARAMETER_VECTOR(mu);          // log peak age, one per transition
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(lh_del);
  PARAMETER_VECTOR(log_c2);      // log curvature of log h at its peak, per transition
  PARAMETER_VECTOR(log_lh_tau);
  PARAMETER_VECTOR(b_coh);
  PARAMETER_VECTOR(itc);
  PARAMETER(d1);                 // log multiplier of the X1 marriage intensity
  PARAMETER(log_lambda);         // log rate of leaving the first phase

  vector<Type> sigma = exp(log_sigma), tau = exp(log_lh_tau);
  // The sharpness gam and the tail exponent tau are nearly the same parameter
  // twice: both set how fast log h falls away from its peak. Pinning gam at a
  // common value across transitions was measured on 2026-10-21 and it lifts the
  // debut tail, because the fitted sharpness differs by two orders of magnitude
  // between transitions. Estimating the PEAK CURVATURE instead separates them,
  // since tau then only sets the far tail power a^(2 tau).
  //   d2 log h / da2 at a = peak  =  -2 gam tau^2 / (cosh(del)^2 peak^2 sigma^2)
  vector<Type> gam(3);
  for (int z = 0; z < 3; z++)
    gam[z] = exp(log_c2[z]) * pow(cosh(lh_del[z]), 2) * exp(Type(2) * mu[z]) *
             sigma[z] * sigma[z] / (Type(2) * tau[z] * tau[z]);
  ADREPORT(gam);
  Type lam = exp(log_lambda);
  ADREPORT(lam);
  ADREPORT(d1);

  if (!nopen) {
    dll -= dnorm(mu[0], Type(3), Type(1), true);
    dll -= dnorm(mu[1], Type(3), Type(1), true);
    dll -= dnorm(mu[2], Type(3), Type(1), true);
    dll -= dnorm(sigma, Type(0), Type(1), true).sum() + log_sigma.sum();
    dll -= dnorm(b_coh, Type(0), Type(1), true).sum();
    dll -= dnorm(lh_del, Type(0), Type(3), true).sum();
    dll -= dnorm(log_c2, Type(0), Type(5), true).sum();
    dll -= dnorm(log_lh_tau, Type(-1.386), Type(1.5), true).sum();
    for (int z = 0; z < 3; z++) dll -= dnorm(itc[z], Type(-2), Type(3), true);
    dll -= dnorm(d1, Type(0), Type(3), true);
    dll -= dnorm(log_lambda, Type(0), Type(1.5), true);
  }

  // The annual transition matrix depends on age and cohort alone, so it is built
  // once for each pair and reused across women. Building it per woman instead
  // tapes the same matrix exponential tens of thousands of times over.
  int nC = coh_u.size();
  std::vector< matrix<Type> > PA(a_max * nC);
  matrix<Type> RS(a_max, nC);          // the V exit rate, for the survival product
  vector<Type> rr(3);
  for (int a = 0; a < a_max; a++) for (int cj = 0; cj < nC; cj++) {
    for (int z = 0; z < 3; z++)
      rr[z] = exp(itc[z] + logSASLH(Type(a + 1),
                                    exp(mu[z] + b_coh[z] * coh_u[cj]),
                                    sigma[z], lh_del[z], gam[z], Type(1), tau[z]));
    matrix<Type> Qd = makeQ_ph(rr[0], rr[1], rr[2] * exp(d1), lam, rr[2]).toDense();
    PA[a * nC + cj] = atomic::expm(Qd);
    RS(a, cj) = rr[0] + rr[1];
  }
  // the matrix for age a (1 based) and cohort column cj
  auto Pm = [&](int a, int cj) -> const matrix<Type>& {
    int aa = a < 1 ? 1 : (a > a_max ? a_max : a);
    return PA[(aa - 1) * nC + cj];
  };

  // mask for the years known to be spent in X, either phase
  vector<Type> mX(4);
  mX << Type(0), Type(1), Type(1), Type(0);

  for (int i = 0; i < afs.size(); i++) {
    int married = afm[i] > 0;
    int debut_first = afs[i] > 0 && (!married || afs[i] < afm[i]);
    int a_leave = debut_first ? afs[i] : (married ? afm[i] : aai[i]);
    int cj = ic[i];

    // The V span, ages 1 to a_leave - 1. V is a single state and its only exits
    // leave it, so the masked propagation is the survival product.
    Type log_s = Type(0);
    for (int a = 1; a < a_leave; a++) {
      int aa = a > a_max ? a_max : a;
      log_s -= RS(aa - 1, cj);
    }
    vector<Type> vp(4); vp.setZero(); vp(0) = exp(log_s);
    if (!married && !debut_first) {                     // type 1
      dll -= n[i] * log(vp(0) + eps);
      continue;
    }
    if (married && !debut_first) {                      // type 3
      vp = (vp.matrix().transpose() * Pm(a_leave, cj)).array();
      dll -= n[i] * log(vp(3) + eps);
      continue;
    }
    // the debut year, then the X span, with the phase split kept between years
    int a_end = married ? afm[i] : (aai[i] > a_leave ? aai[i] : a_leave + 1);
    for (int a = a_leave; a < a_end; a++) {
      vp = (vp.matrix().transpose() * Pm(a, cj)).array();
      vp = vp * mX;
    }
    if (!married) {                                     // type 2
      dll -= n[i] * log(vp(1) + vp(2) + eps);
      continue;
    }
    vp = (vp.matrix().transpose() * Pm(a_end, cj)).array();   // type 4
    dll -= n[i] * log(vp(3) + eps);
  }

  SIMULATE {
    // annual transition matrices and rates by age for the reference cohort
    array<Type> PP(4, 4, a_max);
    matrix<Type> QQ(a_max, 4);
    PP.setZero(); QQ.setZero();
    for (int a = 1; a <= a_max; a++) {
      for (int z = 0; z < 3; z++)
        rr[z] = exp(itc[z] + logSASLH(Type(a), exp(mu[z]), sigma[z], lh_del[z],
                                      gam[z], Type(1), tau[z]));
      QQ(a - 1, 0) = rr[0]; QQ(a - 1, 1) = rr[1];
      QQ(a - 1, 2) = rr[2] * exp(d1); QQ(a - 1, 3) = rr[2];
      matrix<Type> Qd = makeQ_ph(rr[0], rr[1], rr[2] * exp(d1), lam, rr[2]).toDense();
      PP.col(a - 1) = atomic::expm(Qd);
    }
    REPORT(PP);
    REPORT(QQ);
  }

  return dll;
}
