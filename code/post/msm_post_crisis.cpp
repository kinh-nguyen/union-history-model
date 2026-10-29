// Variant msm_post_crisis: msm_post_par with a richer calendar-time crisis term - per-rate centre and
// width, a one-sided width ratio (skew), and an optional second bump. Nests msm_post_par exactly when
// useskew = usebump2 = 0 and the centre/width are tied across rates by the parameter map.
// Base variant msm_post_2d: black-box penalised tensor surface (calendar year x years since first marriage)
// on each log-hazard, plus a penalised smooth in age at marriage. Otherwise identical to msm_post.cpp.
// Post-marriage block of the eight-state model, fitted on its own (it shares no
// parameters with the pre-marriage block). Episodes start in M at age at first
// marriage; the state at interview is current-status data.
//
// Speed: episodes are grouped by path (age at marriage x cohort). Each group's
// chain is run once and every interview duration is read off the same pass,
// and the group (not the row) is the unit handed to the parallel accumulator.
//
// States: 0 M, 1 D1, 2 W1, 3..3+K later union R (K = 0: one state; K > 0:
// tunnel states R_0..R_K = completed years in the later union), 4+K D2, 5+K W2.
#include <TMB.hpp>

template<class Type> Type softplus(Type x) { return log(Type(1) + exp(x)); }

// peak-anchored hyperbola in duration (sinh-arcsinh log-hazard with tau = 1/2), 0 at the peak
template<class Type> Type hyper(Type d, vector<Type> dh, Type wmin) {
  Type L = exp(dh(1)), R = exp(dh(2)), w = wmin + exp(dh(3));
  d = d - exp(dh(0));
  Type d0 = w * (L - R) / (Type(2) * sqrt(L * R));
  return Type(0.5) * (L - R) * d - Type(0.5) * (L + R) * (sqrt(w * w + (d + d0) * (d + d0)) - sqrt(w * w + d0 * d0));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // groups
  DATA_VECTOR(gstart);   // age at marriage (episode start)
  DATA_VECTOR(gcoh);     // cohort, decades from 1975
  DATA_IVECTOR(gdmax);   // longest interview duration in the group
  DATA_IVECTOR(gy0);     // row of PB for calendar year at the start of the episode
  DATA_IVECTOR(ga0);     // row of AB for age at the start of the episode
  DATA_IVECTOR(roff);    // rows of group g are roff(g) .. roff(g+1)-1, sorted by duration
  // rows
  DATA_IVECTOR(rdur); DATA_IVECTOR(rZ); DATA_VECTOR(rn);
  // bases (linear parts removed), computed outside so they match the prototype exactly
  DATA_MATRIX(PB); DATA_MATRIX(AB);
  // switches
  DATA_INTEGER(hump); DATA_INTEGER(frail); DATA_INTEGER(use_per); DATA_INTEGER(per_random);
  DATA_SCALAR(lam); DATA_INTEGER(ties); DATA_INTEGER(K); DATA_INTEGER(agespl);
  DATA_SCALAR(wmin); DATA_INTEGER(use_kappa); DATA_INTEGER(use_coh); DATA_INTEGER(tq);
  DATA_INTEGER(nmax);    // cap on uniformization terms (TMB default 100)
  DATA_MATRIX(PF);       // [1, PB]: period basis of the surface, rows as PB
  DATA_MATRIX(DF);       // [1, dlin, DBt]: duration basis, row j = completed years since first marriage
  DATA_MATRIX(AF);       // age-at-marriage basis (linear part removed), row = round(age at marriage)
  DATA_INTEGER(surf); DATA_SCALAR(lamP); DATA_SCALAR(lamD); DATA_INTEGER(afs); DATA_SCALAR(lamA);
  DATA_SCALAR(ybase);    // calendar year of PB row 0 (1930)
  DATA_MATRIX(MB); DATA_INTEGER(usem); DATA_SCALAR(lamM);    // marriage-year smooth, rows as PB, indexed by gy0
  DATA_INTEGER(usebump); DATA_IVECTOR(bumpadd);              // calendar-time bump; bumpadd(z)=1: added to the hazard
  DATA_INTEGER(useskew); DATA_INTEGER(usebump2);             // one-sided width ratio; second bump
  DATA_INTEGER(useq); DATA_INTEGER(useed);                   // quadratic in year of marriage; early-divorce excess

  PARAMETER_VECTOR(itc);        // divorce, widowhood, remarriage at age 25, duration 10
  PARAMETER_VECTOR(gp_b);       // age slopes
  PARAMETER_VECTOR(b_coh);      // cohort slopes
  PARAMETER(lkap);              // log kappa: later-union divorce multiplier
  PARAMETER_VECTOR(b_tm);       // duration slopes: (shared) or (widowhood, remarriage) with hump
  PARAMETER_VECTOR(dh);         // divorce duration curve: log peak, log rise, log fall, log width
  PARAMETER_VECTOR(fr);         // fragile class: logit pi, log mult divorce, log mult remarriage
  PARAMETER_MATRIX(per);        // period spline coefficients, 3 x NPB
  PARAMETER_VECTOR(log_sd_per); // sd of second differences (random-effect mode)
  PARAMETER_VECTOR(tie);        // log kappa_w (R->W2 vs M->W1), log rho2 (D2,W2->R vs D1->R)
  PARAMETER_VECTOR(ras);        // remarriage: nonlinear age (3 df)
  PARAMETER_MATRIX(A0); PARAMETER_MATRIX(A1); PARAMETER_MATRIX(A2);  // surfaces: divorce, widowhood, remarriage
  PARAMETER_MATRIX(fa);         // 3 x KA age-at-marriage smooths
  PARAMETER_MATRIX(fm);         // 3 x NM marriage-year smooths
  PARAMETER_MATRIX(bmp);        // 3 x 4 by rate: centre - 1990, log sd, amplitude (log where bumpadd), log width ratio
  PARAMETER_MATRIX(bmp2);       // 3 x 3 by rate: second bump, same columns 0-2
  PARAMETER_VECTOR(qc);         // quadratic marriage cohort, ((y0 - 1990)/10)^2, per rate
  PARAMETER_VECTOR(ed);         // early divorce: log A at y0 = 2000, zeta per decade of y0, log theta

  int G = gstart.size(), S = K + 6, NPB = PB.cols(), nk = frail ? 2 : 1;
  int iD2 = 4 + K, iW2 = 5 + K, nY = PB.rows(), nA = AB.rows();
  vector<int> obs(S);
  obs(0) = 0; obs(1) = 1; obs(2) = 2; for (int k = 0; k <= K; k++) obs(3 + k) = 3; obs(iD2) = 4; obs(iW2) = 5;

  Type nll = 0;
  if (use_per) {
    for (int z = 0; z < 3; z++) {
      for (int i = 0; i + 2 < NPB; i++) {
        Type d2 = per(z, i + 2) - Type(2) * per(z, i + 1) + per(z, i);
        if (per_random) nll -= dnorm(d2, Type(0), exp(log_sd_per(z)), true);
        else            nll += lam * d2 * d2;
      }
      for (int i = 0; i < NPB; i++) {
        if (per_random) nll -= dnorm(per(z, i), Type(0), Type(22.360679775), true);  // = ridge 1e-3
        else            nll += Type(1e-3) * per(z, i) * per(z, i);
      }
    }
  }

  std::vector<matrix<Type> > SS(3);
  if (surf) {
    for (int z = 0; z < 3; z++) {
      matrix<Type> A = (z == 0) ? A0 : ((z == 1) ? A1 : A2);
      int nr = A.rows(), nc = A.cols();
      for (int cc = 0; cc < nc; cc++) for (int r = 1; r + 2 < nr; r++) {
        Type d2 = A(r + 2, cc) - Type(2) * A(r + 1, cc) + A(r, cc); nll += lamP * d2 * d2; }
      for (int r = 0; r < nr; r++) for (int cc = 2; cc + 2 < nc; cc++) {
        Type d2 = A(r, cc + 2) - Type(2) * A(r, cc + 1) + A(r, cc); nll += lamD * d2 * d2; }
      for (int r = 0; r < nr; r++) for (int cc = 0; cc < nc; cc++) nll += Type(1e-3) * A(r, cc) * A(r, cc);
      matrix<Type> DT = DF.transpose();
      matrix<Type> PA = PF * A; SS[z] = PA * DT;
    }
  }
  if (afs) for (int z = 0; z < 3; z++) {
    for (int i = 0; i + 2 < fa.cols(); i++) { Type d2 = fa(z, i + 2) - Type(2) * fa(z, i + 1) + fa(z, i); nll += lamA * d2 * d2; }
    for (int i = 0; i < fa.cols(); i++) nll += Type(1e-3) * fa(z, i) * fa(z, i);
  }
  if (usem) for (int z = 0; z < 3; z++) {
    for (int i = 0; i + 2 < fm.cols(); i++) { Type d2 = fm(z, i + 2) - Type(2) * fm(z, i + 1) + fm(z, i); nll += lamM * d2 * d2; }
    for (int i = 0; i < fm.cols(); i++) nll += Type(1e-3) * fm(z, i) * fm(z, i);
  }
  Type kap = use_kappa ? exp(lkap) : Type(1);
  Type kw = ties ? exp(tie(0)) : Type(1), r2 = ties ? exp(tie(1)) : Type(1);
  Type pif = frail ? invlogit(fr(0)) : Type(0);
  sparse_matrix_exponential::config<Type> cfg; cfg.trace = false; cfg.Nmax = nmax;

  parallel_accumulator<Type> dll(this);
  matrix<Type> PR(rdur.size(), 6); PR.setZero();   // predicted cell probabilities per row (REPORT)
  for (int g = 0; g < G; g++) {
    matrix<Type> V(nk, S); V.setZero();
    for (int k = 0; k < nk; k++) V(k, 0) = Type(1);
    Type lg = 0; int r = roff(g);
    for (int j = 0; j < gdmax(g); j++) {
      Type age = gstart(g) + Type(j), tm = softplus(Type(j));
      int yi = gy0(g) + j; yi = yi < 0 ? 0 : (yi >= nY ? nY - 1 : yi);
      int ai = ga0(g) + j; ai = ai < 0 ? 0 : (ai >= nA ? nA - 1 : ai);
      for (int k = 0; k < nk; k++) {
        Type bd = itc(0) + gp_b(0) * (age - Type(25));
        Type ew = itc(1) + gp_b(1) * (age - Type(25));
        Type er = itc(2) + gp_b(2) * (age - Type(25));
        if (use_coh) { bd += b_coh(0) * gcoh(g); ew += b_coh(1) * gcoh(g); er += b_coh(2) * gcoh(g); }
        Type btw = hump ? b_tm(0) : b_tm(0), btr = hump ? b_tm(1) : b_tm(0);
        ew += btw * (tm - Type(10)); er += btr * (tm - Type(10));
        if (agespl) for (int i = 0; i < 3; i++) er += AB(ai, i) * ras(i);
        if (use_per) for (int i = 0; i < NPB; i++) {
          bd += PB(yi, i) * per(0, i); ew += PB(yi, i) * per(1, i); er += PB(yi, i) * per(2, i);
        }
        if (surf) { int jj = j < DF.rows() ? j : DF.rows() - 1;
          bd += SS[0](yi, jj); ew += SS[1](yi, jj); er += SS[2](yi, jj); }
        if (afs) { int ia = (int) std::floor(asDouble(gstart(g)) + 0.5); ia = ia < 0 ? 0 : (ia >= AF.rows() ? AF.rows() - 1 : ia);
          for (int i = 0; i < AF.cols(); i++) { bd += AF(ia, i) * fa(0, i); ew += AF(ia, i) * fa(1, i); er += AF(ia, i) * fa(2, i); } }
        Type y0 = ybase + Type(gy0(g));
        Type exD = 0, exW = 0, exR = 0, edm = 1;
        if (usem) { int im = gy0(g) < 0 ? 0 : (gy0(g) >= MB.rows() ? MB.rows() - 1 : gy0(g));
          for (int i = 0; i < MB.cols(); i++) { bd += MB(im, i) * fm(0, i); ew += MB(im, i) * fm(1, i); er += MB(im, i) * fm(2, i); } }
        if (useq) { Type u = (y0 - Type(1990)) / Type(10); bd += qc(0) * u * u; ew += qc(1) * u * u; er += qc(2) * u * u; }
        if (usebump) {
          Type tc = y0 + Type(j) + Type(0.5);
          Type* lin[3] = {&bd, &ew, &er}; Type* exc[3] = {&exD, &exW, &exR};
          for (int z = 0; z < 3; z++) {
            Type x = (tc - Type(1990) - bmp(z, 0)) / exp(bmp(z, 1));
            if (useskew) { Type k = exp(bmp(z, 3)); x = CppAD::CondExpLt(x, Type(0), x * k, x / k); }
            Type bb = exp(Type(-0.5) * x * x);
            if (usebump2) { Type x2 = (tc - Type(1990) - bmp2(z, 0)) / exp(bmp2(z, 1));
              Type b2 = exp(Type(-0.5) * x2 * x2);
              if (bumpadd(z)) *exc[z] += exp(bmp2(z, 2)) * b2; else *lin[z] += bmp2(z, 2) * b2; }
            if (bumpadd(z)) *exc[z] += exp(bmp(z, 2)) * bb; else *lin[z] += bmp(z, 2) * bb;
          }
        }
        if (useed) edm = Type(1) + exp(ed(0) + ed(1) * (y0 - Type(2000)) / Type(10)) * exp(-(Type(j) + Type(0.5)) / exp(ed(2)));
        if (frail && k == 1) { bd += fr(1); er += fr(2); }
        Type qdM = exp(bd + (hump ? hyper(tm, dh, wmin) : b_tm(0) * (tm - Type(10)))) * edm + exD;
        Type qw = exp(ew) + exW, qr = exp(er) + exR;

        std::vector<Eigen::Triplet<Type> > tl; vector<Type> rs(S); rs.setZero();
        auto add = [&](int a, int b, Type v) {
          if (tq) tl.push_back(Eigen::Triplet<Type>(b, a, v)); else tl.push_back(Eigen::Triplet<Type>(a, b, v));
          rs(a) += v;
        };
        add(0, 1, qdM); add(0, 2, qw); add(1, 3, qr); add(2, 3, qr);
        for (int kk = 0; kk <= K; kk++) {
          Type qdR = qdM;
          if (K > 0) qdR = exp(bd + (hump ? hyper(softplus(Type(kk)), dh, wmin) : b_tm(0) * (softplus(Type(kk)) - Type(10)))) * edm + exD;
          add(3 + kk, iD2, kap * qdR); add(3 + kk, iW2, kw * qw);
        }
        add(iD2, 3, r2 * qr); add(iW2, 3, r2 * qr);
        for (int s = 0; s < S; s++) tl.push_back(Eigen::Triplet<Type>(s, s, -rs(s)));
        Eigen::SparseMatrix<Type> Q(S, S); Q.setFromTriplets(tl.begin(), tl.end());
        sparse_matrix_exponential::expm_generator<Type> pg(Q, cfg);
        vector<Type> v = V.row(k); v = pg(v);
        if (K > 0) {                     // one more completed year in the later union
          vector<Type> o = v;
          v(3) = Type(0);
          for (int kk = 1; kk < K; kk++) v(3 + kk) = o(3 + kk - 1);
          v(3 + K) = o(3 + K - 1) + o(3 + K);
        }
        V.row(k) = v;
      }
      while (r < roff(g + 1) && rdur(r) == j + 1) {
        int z = rZ(r); Type pz = 0;
        for (int k = 0; k < nk; k++) {
          Type wk = frail ? (k == 0 ? Type(1) - pif : pif) : Type(1);
          for (int s = 0; s < S; s++) PR(r, obs(s)) += wk * V(k, s);
        }
        pz = PR(r, z);
        lg -= rn(r) * log(pz + Type(1e-300));
        r++;
      }
    }
    dll += lg;                           // one accumulation per group
  }
  dll += nll;                            // penalty / prior, added once
  REPORT(PR);
  ADREPORT(kap);
  return dll;
}
