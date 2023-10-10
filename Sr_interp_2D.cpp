namespace interp_Sr_2D{
  gsl_spline2d *spline;
  gsl_interp_accel *xacc, *yacc;
}


double Sr_interp_2D(double r, double y){
  using namespace amplitude;
  using namespace interp_Sr_2D;
  using namespace params;
  if(r<minr) return 1;
  if(r>maxr) return 0;
  
  /*
  if(col.compare("pA") == 0){
    return exp(-0.125*sigma0*Anucleus*TA*pow(Sq(r)*Qs02,gamm)*log(1/(r*LambdaQCD)+ec*M_E));
  }

  return exp(-0.25*pow(Sq(r)*Qs02,gamm)*log(1/(r*LambdaQCD)+ec*M_E));
  */

  return gsl_spline2d_eval(spline,r,max(0.0,y-std::log(1/x0)),xacc,yacc);        // ResumBK or KCBK
  //return gsl_spline2d_eval(spline,r,max(0.0,y-std::log(1/(x0*min(1.0,Sq(r))))),xacc,yacc);    // TBK
}
 

 
/*
double Sr_interp_2D(double r, double y){
    using namespace amplitude;
    double answer;
    double Y0 = std::exp(-x0threshold);
    if(y <= Y0) y = 0;
    else y -= Y0;
    return N.DipoleAmplitude(r,y);
}
 */


void init_Sr_interp_2D(){
  using namespace amplitude;
  using namespace interp_Sr_2D;
  double *Srvals=new double[rpoints*ypoints];
  spline=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,ypoints);
  xacc=gsl_interp_accel_alloc();
  yacc=gsl_interp_accel_alloc();
  for(int i=0; i<ypoints; i++){
    for(int j=0; j<rpoints; j++){
      double tmp=max(0.,Sr[i][j]);
      gsl_spline2d_set(spline,Srvals,j,i,tmp);
    }
  }
  gsl_spline2d_init(spline,rvals,yvals,Srvals,rpoints,ypoints);
  delete[] Srvals;
}





void clear_Sr_interp_2D(){
  using namespace interp_Sr_2D;
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}
