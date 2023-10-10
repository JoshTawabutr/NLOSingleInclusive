namespace interp_Sr_1D{
  gsl_spline *spline;
  gsl_interp_accel *acc;
#pragma omp threadprivate(spline,acc)
}


double Sr_interp_1D(double r){
  using namespace amplitude;
  using namespace interp_Sr_1D;
  if(r<minr) return 1;
  if(r>maxr) return 0;
  return gsl_spline_eval(spline,r,acc);
}


void init_Sr_interp_1D(double y){
  using namespace amplitude;
  using namespace interp_Sr_1D;
  double *Srvals=new double[rpoints];
  for(int i=0; i<rpoints; i++){
    double r=rvals[i];
    Srvals[i]=max(0.,Sr_interp_2D(r,y));
  }
  acc=gsl_interp_accel_alloc();
  spline=gsl_spline_alloc(gsl_interp_cspline,rpoints);
  gsl_spline_init(spline,rvals,Srvals,rpoints);
  delete[] Srvals;
}


void clear_Sr_interp_1D(){
  using namespace interp_Sr_1D;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
