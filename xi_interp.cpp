namespace interp_xi{
  gsl_spline2d *spline_I2, *spline_J, *spline_J1, *spline_K3, *spline_H2, *spline_H3, *spline_H5, *spline_Jv, *spline_Jv2, *spline_JJv_xi1;
  gsl_interp_accel *xacc_I2, *yacc_I2, *xacc_J, *yacc_J, *xacc_J1, *yacc_J1, *xacc_K3, *yacc_K3, *xacc_H2, *yacc_H2, *xacc_H3, *yacc_H3, *xacc_H5, *yacc_H5, *xacc_Jv, *yacc_Jv, *xacc_Jv2, *yacc_Jv2,
    *xacc_JJv_xi1, *yacc_JJv_xi1;
}


double I2_interp(double r, double y){
  using namespace interp_xi;
  using namespace amplitude;
  double res=gsl_spline2d_eval(spline_I2,r,max(y,std::log(1/x0)),xacc_I2,yacc_I2);
  return res;
}


double J1_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_J1,r,max(y,std::log(1/x0)),xacc_J1,yacc_J1);
}


double H2_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_H2,r,max(y,std::log(1/x0)),xacc_H2,yacc_H2);
}


double H3_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_H3,r,max(y,std::log(1/x0)),xacc_H3,yacc_H3);
}


double H5_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_H5,r,max(y,std::log(1/x0)),xacc_H5,yacc_H5);
}


double K3_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_K3,r,max(y,std::log(1/x0)),xacc_J,yacc_J);
}


double J_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_J,r,max(y,std::log(1/x0)),xacc_J,yacc_J);
}


double Jv_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_Jv,r,max(y,std::log(1/x0)),xacc_Jv,yacc_Jv);
}


double Jv2_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_Jv2,r,max(y,std::log(1/x0)),xacc_Jv2,yacc_Jv2);
}


double JJv_xi1_interp(double r, double y){
  using namespace interp_xi;
    using namespace amplitude;
  return gsl_spline2d_eval(spline_JJv_xi1,r,max(y,std::log(1/x0)),xacc_JJv_xi1,yacc_JJv_xi1);
}


void init_xi_interp(double xg){
  using namespace amplitude;
  using namespace interp_xi;

  
  double yg=log(1/xg);

  int realypoints=0;
  for(int i=0; i<ypoints; i++){
    realypoints++;
    if(1-(yvals[i]+std::log(1/x0))/yg<1e-5) break;
  }
  //cout << realypoints << endl;    

  double *yvals_tmp=new double[realypoints];
  for(int i=0; i<realypoints; i++){
    yvals_tmp[i]=min(yvals[i]+std::log(1/x0),yg);
  }

  double *I2_vals=new double[rpoints*realypoints];
  spline_I2=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_I2=gsl_interp_accel_alloc();
  yacc_I2=gsl_interp_accel_alloc();

  double *J_vals=new double[rpoints*realypoints];
  spline_J=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_J=gsl_interp_accel_alloc();
  yacc_J=gsl_interp_accel_alloc();
    
  double *J1_vals=new double[rpoints*realypoints];
  spline_J1=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_J1=gsl_interp_accel_alloc();
  yacc_J1=gsl_interp_accel_alloc();
    
  double *H2_vals=new double[rpoints*realypoints];
  spline_H2=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_H2=gsl_interp_accel_alloc();
  yacc_H2=gsl_interp_accel_alloc();
    
  double *H3_vals=new double[rpoints*realypoints];
  spline_H3=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_H3=gsl_interp_accel_alloc();
  yacc_H3=gsl_interp_accel_alloc();
    
  double *H5_vals=new double[rpoints*realypoints];
  spline_H5=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_H5=gsl_interp_accel_alloc();
  yacc_H5=gsl_interp_accel_alloc();
    
  double *K3_vals=new double[rpoints*realypoints];
  spline_K3=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_K3=gsl_interp_accel_alloc();
  yacc_K3=gsl_interp_accel_alloc();

  double *Jv_vals=new double[rpoints*realypoints];
  spline_Jv=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_Jv=gsl_interp_accel_alloc();
  yacc_Jv=gsl_interp_accel_alloc();
    
  double *Jv2_vals=new double[rpoints*realypoints];
  spline_Jv2=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_Jv2=gsl_interp_accel_alloc();
  yacc_Jv2=gsl_interp_accel_alloc();

  double *JJv_xi1_vals=new double[rpoints*realypoints];
  spline_JJv_xi1=gsl_spline2d_alloc(gsl_interp2d_bicubic,rpoints,realypoints);
  xacc_JJv_xi1=gsl_interp_accel_alloc();
  yacc_JJv_xi1=gsl_interp_accel_alloc();
    
   // cout << "check" << endl;

  for(int i=0; i<realypoints; i++){
    double y=yvals_tmp[i];
      //cout << "y = " << y << endl;
    double xi=1-exp(y-yg);
      //cout << "xi = " << xi << endl;
#pragma omp parallel
    {
      init_Sr_interp_1D(y);
      func_tmp::w_x=gsl_integration_workspace_alloc(gsl_maxpoints);
      func_tmp::w_phi=gsl_integration_workspace_alloc(gsl_maxpoints);
    }
#pragma omp parallel for schedule(dynamic)
    for(int j=0; j<rpoints; j++){
      double r=rvals[j];
      double I2_tmp=(with_CF ? I2(r,xi) : 0);
      double J_tmp=(with_Nc ? J(r,xi) : 0);
      double J1_tmp=(with_gl ? J1(r,xi) : 0);
      double H2_tmp=(with_gg ? H2(r,xi) : 0);
      double H3_tmp=(with_gg ? H3(r,xi) : 0);
      double H5_tmp=(with_gg ? H5(r,xi) : 0);
      double K3_tmp=(with_gq ? K3(r,xi) : 0);
      double Jv_tmp=(with_Nc ? Jv(r,xi) : 0);
      double Jv2_tmp=(with_Nc ? Jv2(r,xi) : 0);
      double JJv_xi1_tmp=(with_xi1 ? JJv_xi1(r) : 0);
      gsl_spline2d_set(spline_I2,I2_vals,j,i,I2_tmp);
      gsl_spline2d_set(spline_J,J_vals,j,i,J_tmp);
      gsl_spline2d_set(spline_J1,J1_vals,j,i,J1_tmp);
      gsl_spline2d_set(spline_H2,H2_vals,j,i,H2_tmp);
      gsl_spline2d_set(spline_H3,H3_vals,j,i,H3_tmp);
      gsl_spline2d_set(spline_H5,H5_vals,j,i,H5_tmp);
      gsl_spline2d_set(spline_K3,K3_vals,j,i,K3_tmp);
      gsl_spline2d_set(spline_Jv,Jv_vals,j,i,Jv_tmp);
      gsl_spline2d_set(spline_Jv2,Jv2_vals,j,i,Jv2_tmp);
      gsl_spline2d_set(spline_JJv_xi1,JJv_xi1_vals,j,i,JJv_xi1_tmp);

	//cout << i << "," << j << "," << J1_tmp << endl;

    }
#pragma omp parallel
    {
      clear_Sr_interp_1D();
      gsl_integration_workspace_free(func_tmp::w_x);
      gsl_integration_workspace_free(func_tmp::w_phi);
    }
  }

  gsl_spline2d_init(spline_I2,rvals,yvals_tmp,I2_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_J,rvals,yvals_tmp,J_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_J1,rvals,yvals_tmp,J1_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_H2,rvals,yvals_tmp,H2_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_H3,rvals,yvals_tmp,H3_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_H5,rvals,yvals_tmp,H5_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_K3,rvals,yvals_tmp,K3_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_Jv,rvals,yvals_tmp,Jv_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_Jv2,rvals,yvals_tmp,Jv2_vals,rpoints,realypoints);
  gsl_spline2d_init(spline_JJv_xi1,rvals,yvals_tmp,JJv_xi1_vals,
                    rpoints,realypoints);
    

  delete[] yvals_tmp;
  delete[] I2_vals;
  delete[] J_vals;
  delete[] J1_vals;
  delete[] H2_vals;
  delete[] H3_vals;
  delete[] H5_vals;
  delete[] K3_vals;
  delete[] Jv_vals;
  delete[] Jv2_vals;
  delete[] JJv_xi1_vals;
}


void clear_xi_interp(){
  using namespace interp_xi;

  gsl_spline2d_free(spline_I2);
  gsl_interp_accel_free(xacc_I2);
  gsl_interp_accel_free(yacc_I2);

  gsl_spline2d_free(spline_J);
  gsl_interp_accel_free(xacc_J);
  gsl_interp_accel_free(yacc_J);
    
  gsl_spline2d_free(spline_J1);
  gsl_interp_accel_free(xacc_J1);
  gsl_interp_accel_free(yacc_J1);
    
  gsl_spline2d_free(spline_H2);
  gsl_interp_accel_free(xacc_H2);
  gsl_interp_accel_free(yacc_H2);
    
  gsl_spline2d_free(spline_H3);
  gsl_interp_accel_free(xacc_H3);
  gsl_interp_accel_free(yacc_H3);
    
  gsl_spline2d_free(spline_H5);
  gsl_interp_accel_free(xacc_H5);
  gsl_interp_accel_free(yacc_H5);
    
  gsl_spline2d_free(spline_K3);
  gsl_interp_accel_free(xacc_K3);
  gsl_interp_accel_free(yacc_K3);

  gsl_spline2d_free(spline_Jv);
  gsl_interp_accel_free(xacc_Jv);
  gsl_interp_accel_free(yacc_Jv);
    
  gsl_spline2d_free(spline_Jv2);
  gsl_interp_accel_free(xacc_Jv2);
  gsl_interp_accel_free(yacc_Jv2);

  gsl_spline2d_free(spline_JJv_xi1);
  gsl_interp_accel_free(xacc_JJv_xi1);
  gsl_interp_accel_free(yacc_JJv_xi1);
}
