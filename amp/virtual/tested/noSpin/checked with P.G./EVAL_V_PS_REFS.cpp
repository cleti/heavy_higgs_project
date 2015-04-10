

FV const& p1 = ps.p1();
FV const& p2 = ps.p2();
FV const& k1 = ps.k1();
FV const& k2 = ps.k2();

double s = 2.0*sp(p1,p2);
double beta_y = 2.0*(sp(p1,k2)-sp(p1,k1))/s;
double beta = sqrt(1.0-4.0/s);
double y = beta_y/beta;
