namespace interp_xi_int{
  gsl_spline *spline;
  gsl_interp_accel *acc;
}


double xi_int_interp(double r){
  using namespace amplitude;
  using namespace interp_xi_int;
  if(r<minr) r=minr;
  if(r>maxr) return 0;
  return gsl_spline_eval(spline,r,acc);
}


double I1(double r, double y){
  using namespace params;  
  double res=(Sr_interp_2D(r,y)*(2*log(c0/r)-log(mu2)))/(2*M_PI);
  if(alpha_s_running==PARENT || alpha_s_running==DAUGHTER || alpha_s_running==SMALLEST){
    return alpha_s_pos(r)*res;
  }else{
    return res;
  }
}


double H1(double r, double y){
  using namespace params;
  double res=(Sq(Sr_interp_2D(r,y))*(2*log(c0/r)-log(mu2)))/M_PI;
  if(alpha_s_running==PARENT || alpha_s_running==DAUGHTER || alpha_s_running==SMALLEST){
    return alpha_s_pos(r)*res;
  }else{
    return res;
  }
}


double integrand_xi(double xi, void *userdata){
  using namespace params;
  double r=((double *)userdata)[0];
  double xp=((double *)userdata)[1];
  double xg=((double *)userdata)[2];
  double k=((double *)userdata)[3];
  //
  double X=xg/(1-xi);
  double y=log(1/X); if(y < 0) y = 0;
  //
  double xi2=Sq(xi);
  double xq_xi=(xi>xp ? xf(xp/xi,mu2) : 0);
  double xq_xi1=xf(xp,mu2);
  double xg_xi=(xi>xp ? xf(xp/xi,mu2) : 0);
  double xg_xi1=xf(xp,mu2);
  //
  double res=0;
  if(with_CF){
    double real_CF=0;
    double virt_CF=-M_1_PI*(2*log(k) - log(mu2) + 2*log(1-xi))*Sr_interp_2D(r,y);
    if(alpha_s_running==PARENT || alpha_s_running==DAUGHTER || alpha_s_running==SMALLEST){
        virt_CF*=alpha_s_pos(r);
        real_CF+=(xi>xp ? I1(r,y)+I1(xi*r,y) - 4*I2_interp(r,y) : 0);
    }
    else if(alpha_s_running==MIXED || alpha_s_running==MIXEDBD){
        virt_CF*=alpha_s_mom(k);
        real_CF+=(xi>xp ? (I1(r,y)+I1(xi*r,y) - 4*I2_interp(r,y))*alpha_s_mom(k) : 0);
    }
    else{
        real_CF+=(xi>xp ? I1(r,y) + I1(xi*r,y) - 4*I2_interp(r,y) : 0);
    }
    res+=(CF*(1+xi2)*(real_CF*xq_xi+virt_CF*xq_xi1))/(1-xi);
  }
  if(with_Nc){
    double real_Nc=(xi>xp ? J_interp(r,y) : 0);
    double virt_Nc=-Jv_interp(r,y);
    double resNc=real_Nc*xq_xi+virt_Nc*xq_xi1;
    res+=(Nc*(1+xi2)*resNc)/(1-xi);
  }
  if(with_gl){
    double resgl1 = (xi>xp ? 0.25*H1(r,y) + 0.5*I1(xi*r,y) : 0);
    double resgl2 = (xi>xp ? -J1_interp(r,y) : 0);
    if(alpha_s_running==PARENT){
        resgl2*=alpha_s_pos(r);
    }
    else if(alpha_s_running==MIXED || alpha_s_running==MIXEDBD){
        resgl1*=alpha_s_pos(k);
        resgl2*=alpha_s_pos(k);
    }
    res += ((resgl1+resgl2)*Nc*xq_xi*(1+Sq(1-xi)))/xi;
}
  if(with_gq){
    double resgq1 = (xi>xp ? 0.25*H1(xi*r,y) + 0.5*I1(r,y) : 0);
    double resgq2 = (xi>xp ? - K3_interp(r,y) : 0);
    if(alpha_s_running==MIXED || alpha_s_running==MIXEDBD){
        resgq1 *= alpha_s_mom(k);
        if(alpha_s_running == MIXEDBD){
            resgq2 *= alpha_s_mom(k);
        }
    }
    double resgq = resgq1 + resgq2;
    //resgq *= Nf;               // Multiply by Nf to normalize. Technically, no multiplication but sum over flavors with fragmentation fn.
    res += (resgq*xg_xi*(Sq(xi)+Sq(1-xi)));
  }
  if(with_gg){
    double gg11 = (xi>xp ? H1(r,y) + H1(xi*r,y) : 0);
    double gg12 = (xi>xp ? - H2_interp(r,y) : 0);
    double gg2 = H3_interp(r,y) - M_1_PI*(2*log(k)-log(mu2)+2*log(1-xi))*Sq(Sr_interp_2D(r,y));
    double gg3 = H5_interp(r,y) - 0.5*M_1_PI*(2*log(k)-log(mu2)+2*log(1-xi))*Sq(Sr_interp_2D(r,y));
    //  double gg2 = H3_interp(r,y) - 2*H4_interp(r,y) - M_1_PI*(2*log(k)-log(mu2)+2*log(1-xi))*Sq(Sr_interp_2D(r,y));
    //  double gg3 = H5_interp(r,y) - H4_interp(r,y) - 0.5*M_1_PI*(2*log(k)-log(mu2)+2*log(1-xi))*Sq(Sr_interp_2D(r,y));
    if(alpha_s_running == PARENT){
        gg2 *= alpha_s_pos(r);
        gg3 *= alpha_s_pos(r);
    }
    else if(alpha_s_running == MIXED || alpha_s_running == MIXEDBD){
        gg2 *= alpha_s_mom(k);
        gg3 *= alpha_s_mom(k);
        gg11 *= alpha_s_mom(k);
        if (alpha_s_running == MIXEDBD){
            gg12 *= alpha_s_mom(k);
        }
    }
    double resgg = (Nc*xg_xi*(gg11+gg12)*Sq(1-(xi*(1-xi)))/(xi*(1-xi))) + (Nc*xg_xi1*gg2*((2*xi/(1-xi))+(xi*(1-xi))))  + (Nf*xg_xi1*gg3*(Sq(xi)+Sq(1-xi)));
    res += resgg;
  }
  if(with_xi1){
    double sub=Nc*2*JJv_xi1_interp(r,y)*xq_xi1 / (1-xi);
    if(alpha_s_running==PARENT) sub*=alpha_s_pos(r);
    res+=(with_Nc ? -sub : sub);
  }
  //res/=1-xi;
  if(gsl_finite(res)==1){
    return res;
  }else{
    return 0;
  }
}


double integral_xi(double r, double xp, double xg, double k){
  using namespace amplitude;
  double userdata[4];
  userdata[0]=r;
  userdata[1]=xp;
  userdata[2]=xg;
  userdata[3]=k;
  gsl_function F;
  F.function=&integrand_xi;
  F.params=userdata;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(gsl_maxpoints);
  double result, error;
    
    gsl_integration_qags(&F,0,1-xg,0,epsrel_gsl,gsl_maxpoints,w,&result,&error);
    //gsl_integration_qags(&F,0,1-xg/x0,0,epsrel_gsl,gsl_maxpoints,w,&result,&error);
    //gsl_integration_qags(&F,0,1-xg/(std::exp(-Y0threshold)),0,epsrel_gsl,gsl_maxpoints,w,&result,&error);
    //gsl_integration_qags(&F,1-xg/(std::exp(-Y0threshold)),1-xg/x0,0,epsrel_gsl,gsl_maxpoints,w,&result,&error);
    
  gsl_integration_workspace_free(w);
     
  return result;
}


void init_xi_int_interp(double xp, double xg, double k){
  using namespace amplitude;
  using namespace interp_xi_int;
  double *intvals=new double[rpoints];
  for(int i=0; i<rpoints; i++){
    double r=rvals[i];
    intvals[i]=integral_xi(r,xp,xg,k);
  }
  acc=gsl_interp_accel_alloc();
  spline=gsl_spline_alloc(gsl_interp_cspline,rpoints);
  gsl_spline_init(spline,rvals,intvals,rpoints);
  delete[] intvals;
}


void clear_xi_int_interp(){
  using namespace interp_xi_int;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
