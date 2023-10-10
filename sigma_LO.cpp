namespace sigma_LO_k_tmp{
  double k_tmp, xp_tmp;
}


double sigma_LO_r(double r, double xp){
  using namespace params;
  if(with_Nc||with_CF){
    return xf(xp,mu2)*Sr_interp_1D(r);
  }
  else if(with_gg){
    return xf(xp,mu2)*Sq(Sr_interp_1D(r));
  }
  return 0;
}


double integrand_sigma_LO_k(double r){
  using namespace sigma_LO_k_tmp;
  return r*gsl_sf_bessel_J0(k_tmp*r)*sigma_LO_r(r,xp_tmp);
}


double sigma_LO_k(double k, double xp){
  using namespace sigma_LO_k_tmp;
  k_tmp=k;
  xp_tmp=xp;
  double integral, error;
  intdeo(integrand_sigma_LO_k,0,k,epsrel_intde,&integral,&error);
  return integral/(2*M_PI);
}

/*
double integrand_sigma_LO_p(double z){
    using namespace params;
    using namespace sigma_LO_p_tmp;
    
    // Compute FF(z)
    fragini_.fini=0;
    double FF;
    int ih_tmp = ih;
    int ic_tmp = ic;
    int io_tmp = io;
    int parton_tmp = parton;
    double scalesqr_tmp = scalesqr;
    fdss_(ih_tmp, ic_tmp, io_tmp, z, scalesqr_tmp, parton_tmp, FF);
    
    // Compute sigma_LO_k(p/z, tau/z)
    double k = p_tmp/z;
    mu2=Sq(k);
    double xp=(k/SQRTS)*exp(yh);
    double xg=(k/SQRTS)*exp(-yh);
    init_interp(xp,xg,k);
    double mult = sigma_LO_k(k,xp);
    clear_interp();
    
    return mult*FF/Sq(z);
}


double sigma_LO_p(double p){
    using namespace sigma_LO_p_tmp;
    p_tmp = p;
    double integral, error;
    double zmin = 0.05;
    intde(integrand_sigma_LO_p, zmin, 1, epsrel_intde, &integral, &error);
    return integral;
}
*/
