

FV const& p1 = ps.p1();
FV const& p2 = ps.p2();
FV const& k1 = ps.k1();
FV const& k2 = ps.k2();
FV const& p3 = ps.p3();
FV const& s1 = ps.s1();
FV const& s2 = ps.s2();

double s_tt = 2.0+2.0*sp(k1,k2);
double beta_tt2 = 1.0-4.0/s_tt;
