
 double Sr_0(double r){
     using namespace params;
     if(col.compare("pA") == 0){
         return exp(-0.125*sigma0*Anucleus*TA*pow(Sq(r)*Qs02,gamm)*log(1/(r*LambdaQCD)+ec*M_E));
     }
     return exp(-0.25*pow(Sq(r)*Qs02,gamm)*log(1/(r*LambdaQCD)+ec*M_E));
}


/*
double Sr_0(double r){
    using namespace amplitude;
 return Sr_interp_2D(r, Y0threshold);
}
 */


namespace amplitude{
  double minr, maxr, r_mult, x0, Y0threshold;
  int rpoints, ypoints;
  double *yvals, *rvals;
  vector<vector <double> > Sr;  // Sr[Y][r]=S(Y,r)
}


/*
void read_bkgrid(vector<double>& yvals_tmp){
  using namespace amplitude;

  AmplitudeLib N(datadir);

  //These are hard-coded by looking through dipole-resumbk-hera-parent-0.00.dip. 
  //A better method is needed. (Not sure if one could search through a .dip file as if it were .dat)
  //The mins are because the nlodisfit code screams if we call for r>=25.2552 or y>=20.8.)

  minr=1.000000000000000e-06;
  r_mult=1.187876286287484e+00;
  rpoints=98;
  x0=0.00000;

  double multFactor = 0.4; 
  vector<double> tmpvec;
  for(int i = 0; i < 52; i++){
    double this_y = multFactor * i;
    yvals_tmp.push_back(this_y);
    double this_r = minr;
    for(int j = 0; j <= rpoints; j++){  
      tmpvec.push_back(1 - N.DipoleAmplitude(this_r, this_y));
      this_r *= r_mult;
    }
    Sr.push_back(tmpvec);
    tmpvec.clear();
  }

}
*/



void read_bkgrid(vector<double>& yvals_tmp){
    using namespace amplitude;
    using namespace params;
    
    string bksol;
    if(col.compare("pA") == 0){
        if(b < 10){
            bksol = bksolpA + to_string(b).substr(0,1);
        }
        else{
            bksol = bksolpA + to_string(b).substr(0,2);
        }
    }
    else{
        bksol = bksolpp;
    }
   
  
        
  ifstream datafile(bksol.c_str());
  if(!datafile.is_open()){
    cerr << "Error opening the BK solution file" << endl;
    exit(1);
  }
  int confid=0;
  while(!datafile.eof() and confid<4){
    string line;
    getline(datafile,line);
    if(line.substr(0,3)=="###"){
      switch(confid){
      case 0:
        minr=str_to_double(line.substr(3,line.length()-3));
        break;
      case 1:
        r_mult=str_to_double(line.substr(3,line.length()-3));
        break;
      case 2:
        rpoints=str_to_int(line.substr(3,line.length()-3));
        break;
      case 3:
              
        x0=str_to_double(line.substr(3,line.length()-3));
              
        //x0=1.;
        //x0=str_to_double(line.substr(3,line.length()-3));
              
        //x0 = std::exp(-x0); // x0 when read from older HM datafiles is actually evolution rapidity at the initial condition, so the corresponding x is exp(-x0); comment this out when reading from DB or the new HM datafile.
              
        break;
      }
      confid++; 
    }
  }
  vector<double> tmpvec;
  while(!datafile.eof()){
    string line;
    getline(datafile,line);
    if(line.substr(0,3)=="###"){
      if(tmpvec.size()>0) Sr.push_back(tmpvec);
      yvals_tmp.push_back(str_to_double(line.substr(3,line.length()-3)));
      tmpvec.clear();
      continue;
    }
    tmpvec.push_back(1-str_to_double(line));
  }
  Sr.push_back(tmpvec);
  datafile.close();
}



void init_bksol(){
  using namespace amplitude;

  vector<double> yvals_tmp;

  read_bkgrid(yvals_tmp);
  ypoints=yvals_tmp.size();

  yvals=new double[ypoints];
  for(int i=0; i<ypoints; i++){
    yvals[i]=yvals_tmp[i];

  }
  yvals_tmp.clear();
  

  rvals=new double[rpoints];
  for(int i=0; i<rpoints; i++){
    rvals[i]=minr*gsl_pow_int(r_mult,i);
  }
  maxr=rvals[rpoints-1];
  minlnr=log(0.5*amplitude::minr);
  maxlnr=log(2*amplitude::maxr);
}
