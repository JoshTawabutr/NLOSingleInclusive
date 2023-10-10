namespace func_tmp{
  gsl_integration_workspace *w_x, *w_phi;
#pragma omp threadprivate(w_x,w_phi)
}


double integrand_x(double lnx, void *userdata){
  double r=((double *)userdata)[0];
  double xi=((double *)userdata)[1];
  double phi=((double *)userdata)[2];
  double flag=((double *)userdata)[3];
  // flag=1: J, flag=-1: Jv, flag=-2: Jv2, flag=0: J-Jv(xi=1), flag=2: I2, flag=3: J1,
  // flag=4: H2, flag=5: H3, flag=6: H5, flag=7: K3

  double x=exp(lnx), x2=Sq(x), r2=Sq(r), sprx=r*x*cos(phi), rpx2=r2+x2+2*sprx, rmx2=r2+x2-2*sprx;

  double res;
  if(flag>6.5){
    // K3
    double dip1 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    double K = (x2+sprx)/(x2*rpx2);
    res = 2 * exp(2*lnx) * K * Sr_interp_1D(xi*x) * Sr_interp_1D(dip1);
    if(alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    else if(alpha_s_running==MIXED){
      res *= alpha_s_pos(min(xi*x, xi*sqrt(rpx2)));
    }
    else if(alpha_s_running==PARENT){
      res *= alpha_s_pos(r);
    }
  }
  else if(flag>5.5){
    // H5
    double dip1 = sqrt(r2 + Sq(xi)*x2 - 2*xi*sprx);
    double dip2 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    res = 2 * (Sr_interp_1D(dip2) * Sr_interp_1D(dip1) - Sq(Sr_interp_1D(dip2)));
    if(alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
  }
  else if(flag>4.5){
    // H3
    double dip1 = sqrt(r2 + Sq(xi)*x2 + 2*xi*sprx);
    double dip2 = sqrt(r2 + Sq(1-xi)*x2 - 2*(1-xi)*sprx);
    double dip3 = sqrt(r2 + Sq(1-xi)*x2 + 2*(1-xi)*sprx);
    res = 4 * (Sr_interp_1D(x) * Sr_interp_1D(dip1) * Sr_interp_1D(dip2) - Sq(Sr_interp_1D(dip3)));
    if(alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
  }
  else if(flag>3.5){
    // H2
    double dip1 = sqrt(Sq(xi)*rpx2);
    double dip2 = sqrt(Sq(xi)*r2 + Sq(1-xi)*x2 - 2*xi*(1-xi)*sprx);
    double K = (x2+sprx)/(x2*rpx2);
    res = 8 * exp(2*lnx) * K * Sr_interp_1D(x) * Sr_interp_1D(dip1) * Sr_interp_1D(dip2);
    if(alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
    else if(alpha_s_running==PARENT){
        res *= alpha_s_pos(r);
    }
    else if(alpha_s_running==MIXED){
        res *= alpha_s_pos(min(xi*x, xi*sqrt(rpx2)));
    }
  }
  else if(flag>2.5){
    // J1
    double dip1 = sqrt(Sq(xi)*r2+Sq(1-xi)*x2-2*xi*(1-xi)*sprx);
    double K = (x2+sprx)/(x2*rpx2);
    res = 2 * exp(2*lnx) * K * Sr_interp_1D(x) * Sr_interp_1D(dip1);
      /*
    if(alpha_s_running==MIXED){
      res*=alpha_s_pos(min(xi*x,xi*sqrt(rpx2)));
    }
    else if(alpha_s_running==MIXEDBD){
      res*=alpha_s_pos(x);
    }
      */
  }
  else if(flag>1.5){
    // I2
    double dip1=sqrt(Sq(xi)*r2+Sq(1-xi)*x2-2*xi*(1-xi)*sprx);
    double dip2=sqrt(rpx2);
    double K=(x2+sprx)/(x2*rpx2);
    res=exp(2*lnx)*K*Sr_interp_1D(dip1);
    if(alpha_s_running==PARENT){
      res*=alpha_s_pos(r);
    }
    else if(alpha_s_running==DAUGHTER || alpha_s_running==SMALLEST){
      res*=alpha_s_pos(dip2);
    }
      
  }else if(flag>0.5){
    // J
    double dip1=sqrt(r2+Sq(1-xi)*x2+2*(1-xi)*sprx);
    double dip2=xi*x;
    double dip3=sqrt(rpx2);
    double K=(x2+sprx)/(x2*rpx2);
    if(alpha_s_running==PARENT){
      res=alpha_s_pos(r)*exp(2*lnx)*K*(Sr_interp_1D(dip1)-Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==DAUGHTER){
      res=alpha_s_pos(x)*exp(2*lnx)*K*(Sr_interp_1D(dip1)-Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==SMALLEST){
      res=exp(2*lnx)*K*(alpha_s_pos(x)*Sr_interp_1D(dip1)-alpha_s_pos(min(r, min(xi*x, xi*dip3)))*Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==MIXED){
      res=exp(2*lnx)*K*alpha_s_pos(min(xi*x,xi*dip3))*(Sr_interp_1D(dip1)-Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==MIXEDBD){
      res=exp(2*lnx)*K*alpha_s_pos(x)*(Sr_interp_1D(dip1)-Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else{
      res=exp(2*lnx)*K*(Sr_interp_1D(dip1)-Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
      
  }else if(flag>-0.5){
    // J-Jv(xi=1)
    double K;
    if(x<amplitude::minr/10 || sqrt(rpx2)<amplitude::minr/10){
      K=0;
    }else{
      K=r2/(x2*rpx2+1e-20);
    }
    res=exp(2*lnx)*K*(Sr_interp_1D(x)*Sr_interp_1D(sqrt(rpx2))-Sr_interp_1D(r));
    if(alpha_s_running==DAUGHTER){
      res*=alpha_s_pos(x);
    }
  }else if(flag>-1.5){
    // Jv, now include Jv2 stuffs already
    double dip1=sqrt(r2+Sq(1-xi)*x2-2*(1-xi)*sprx);
    double dip2=x;
    double dip3=sqrt(Sq(xi)*x2+r2+2*xi*sprx);
    double K=1/x2;
    if(alpha_s_running==PARENT){
      // take 1
      res=exp(2*lnx)*K*(alpha_s_pos(dip1)*Sr_interp_1D(dip1) - alpha_s_pos(r)*Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
      // take 2
      //res=exp(2*lnx)*K*alpha_s_pos(r)*(Sr_interp_1D(dip1) - Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==DAUGHTER){
      res=exp(2*lnx)*K*(alpha_s_pos(dip1)*Sr_interp_1D(dip1) - alpha_s_pos(x)*Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==SMALLEST){
      // take 1
      res=exp(2*lnx)*K*(alpha_s_pos(dip1)*Sr_interp_1D(dip1) - alpha_s_pos(min(r, min(xi*x, dip3)))*Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
      // take 2
      //res=exp(2*lnx)*K*alpha_s_pos(min(r, min(xi*x, dip3)))*(Sr_interp_1D(dip1) - Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==MIXED){
      res=exp(2*lnx)*K*alpha_s_pos(min(xi*x,xi*x+r))*(Sr_interp_1D(dip1) - Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else if(alpha_s_running==MIXEDBD){
      res=exp(2*lnx)*K*alpha_s_pos(x)*(Sr_interp_1D(dip1) - Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
    else{
      res = exp(2*lnx)*K*(Sr_interp_1D(dip1) - Sr_interp_1D(dip2)*Sr_interp_1D(dip3));
    }
      
  }else{
      // Jv2
      double dip1=sqrt(r2+Sq(1-xi)*x2-2*(1-xi)*sprx);
      double dip2=x;
      double dip3=sqrt(Sq(xi)*x2+r2+2*xi*sprx);
      double K=1/x2;
      if(alpha_s_running==PARENT){
        // take 1
        res=-exp(2*lnx)*K*alpha_s_pos(r)*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
        // take 2
        //res=-exp(2*lnx)*K*alpha_s_pos(r)*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
      }
      else if(alpha_s_running==DAUGHTER){
        res=-exp(2*lnx)*K*alpha_s_pos(x)*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
      }
      else if(alpha_s_running==SMALLEST){
        // take 1
        res=-exp(2*lnx)*K*alpha_s_pos(min(r, min(xi*x, dip3)))*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
        // take 2
        //res=-exp(2*lnx)*K*alpha_s_pos(min(r, min(xi*x, dip3)))*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
      }
      else if(alpha_s_running==MIXED){
        res=-exp(2*lnx)*K*alpha_s_pos(min(xi*x, dip3))*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
      }
      else{
          res=-exp(2*lnx)*K*Sr_interp_1D(dip2)*Sr_interp_1D(dip3);
      }
        
    }
  //if(alpha_s_running==DAUGHTER){
  //  res*=alpha_s_pos(x);
  //}
  if(gsl_finite(res)==1){
    return res;
  }else{
    return 0;
  }
}


double integrand_phi(double phi, void *userdata){
  using namespace amplitude;
  ((double *)userdata)[2]=phi;
  double result, error;
  gsl_function F;
  F.function=&integrand_x;
  F.params=userdata;
  gsl_integration_qag(&F,minlnr,maxlnr,epsabs_gsl,epsrel_gsl,gsl_maxpoints,
                      GSL_INTEG_GAUSS15,func_tmp::w_x,&result,&error);
  return result;
}


double func(double r, double xi, double flag){
  double userdata[4];
  userdata[0]=r;
  userdata[1]=xi;
  userdata[3]=flag;
  double result, error;
  gsl_function F;
  F.function=&integrand_phi;
  F.params=userdata;
  gsl_integration_qag(&F,0,M_PI,epsabs_gsl,epsrel_gsl,gsl_maxpoints,
                      GSL_INTEG_GAUSS15,func_tmp::w_phi,&result,&error);
  return (2*result)/Sq(2*M_PI); // 2: int over pi instead of 2*pi
}


double I2(double r, double xi){
  return func(r,xi,2);
}


double J(double r, double xi){
  return 2*func(r,xi,1);
}


double J1(double r, double xi){
  return func(r,xi,3);
}
    
    
double H2(double r, double xi){
  return func(r,xi,4);
}
    
    
double H3(double r, double xi){
  return func(r,xi,5);
}
    
    
double H5(double r, double xi){
  return func(r,xi,6);
}
    
    
double K3(double r, double xi){
  return func(r,xi,7);
}


double Jv(double r, double xi){
  return 2*func(r,xi,-1);
}


double Jv2(double r, double xi){
  return 2*func(r,xi,-2);
}


double JJv_xi1(double r){
  return func(r,1,0);
}
