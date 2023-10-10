namespace pdf{
  LHAPDF::PDF* lhpdf;
}


double xf(double x, double Q2){
    using namespace pdf;
    using namespace params;
    
    if(x<0 || x>1){
        cerr << "Bad x= " << x << endl;
        exit(1);
    }
    
    // 0: gluon
    // 1: d, -1: dbar
    // 2: u, -2: ubar
    // 3: s, -3: sbar
    int particle_flag;
    if(incoming.compare("g") == 0) particle_flag = 0;
    else if(incoming.compare("d") == 0) particle_flag = 1;
    else if(incoming.compare("u") == 0) particle_flag = 2;
    else if(incoming.compare("s") == 0) particle_flag = 3;
    else if(incoming.compare("dbar") == 0) particle_flag = -1;
    else if(incoming.compare("ubar") == 0) particle_flag = -2;
    else if(incoming.compare("sbar") == 0) particle_flag = -3;
    
    return lhpdf->xfxQ2(particle_flag,x,Q2);
}


/*
double xq(double x, double Q2){
  using namespace pdf;
  if(x<0 || x>1){
    cerr << "Bad x= " << x << endl;
    exit(1);
  }
  return lhpdf->xfxQ2(2,x,Q2);    // 2 refers to up quark
}


// Gluon PDF from LHAPDF
double xgl(double x, double Q2){
  using namespace pdf;
  if(x<0 || x>1){
    cerr << "Bad x= " << x << endl;
    exit(1);
  }
  return lhpdf->xfxQ2(0,x,Q2);    // 0 refers to gluon
}
*/

/*
double alpha_s_pos(double r){
  return (2*M_PI)/(beta0*log(1/(r*LambdaQCD)+alpha_s_mu_0));
}
 */


// Coupling constant runs differently in Henri's dipole. Switch to this one for position space rc.
double alpha_s_pos(double r){
  return (4*M_PI)/(0.2*beta0*log(pow(2.5, 10) + pow(4*1.21/Sq(LambdaQCD*r), 5)));
}



double alpha_s_mom(double Q2){
  return (4*M_PI)/(beta0*log(Q2/Sq(LambdaQCD)));
}




double str_to_double(string str){
  stringstream stream(str);
  double tmp;
  stream >> tmp;
  return tmp;
}


int str_to_int(string str){
  stringstream stream(str);
  int tmp;
  stream >> tmp;
  return tmp;
}


void gsl_error_handler(const char * reason,
                       const char * file,
                       int line,
                       int error){
  if(error==GSL_EMAXITER || error==GSL_ETOL || error==GSL_EUNDRFLW || error==GSL_EROUND){
    return;
  }else{
    std::cerr << file << ":" << line << ": Error " << error << ": " << reason << std::endl;
    exit(1);
  }
}


/*
double WSnormalize(){
    using namespace params;
    double integral, error;
    intdei(WSnormalize1D,0,epsrel_intde,&integral,&error);
    return 1/(2 * integral);
}


double WSnormalize1D(double z){
    using namespace params;
    double userdata[1];
    userdata[0] = z;
    gsl_function F;
    F.function = &WSnormintegrand;
    F.params = userdata;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(gsl_maxpoints);
    double result, error;
    gsl_integration_qagiu(&F,0,0,epsrel_gsl,gsl_maxpoints,w,&result,&error);
    gsl_integration_workspace_free(w);
    return result;
}


double WSnormintegrand(double b, void *userdata){
    using namespace params;
    double z = ((double *)userdata)[0];
    double denom = 1 + exp((pow(Sq(b) + Sq(z), 0.5) - RA) / WSd);
    return (2 * M_PI * b) / denom;
}
*/


