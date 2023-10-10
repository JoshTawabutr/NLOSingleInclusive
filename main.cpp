#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include "gsl/gsl_math.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>

//#include </Users/tawabyx/Documents/GitHub/nlodisfit/src/amplitudelib.hpp>

#include <LHAPDF/LHAPDF.h>

using namespace std;

static inline double Sq(double x){return x*x;}

enum running_types {FIXED, MOM, PARENT, DAUGHTER, SMALLEST, MIXED, MIXEDBD};

namespace params{
  // First Pb BK solution
  //const string bksolpA = "./Pb/Pb_smallest_b_";
  //const string bksolpp = "./Pb/proton_resumbk_fit1_smallest";

  // KCBK fit1 solution
  //const string bksolpA = "./kcbk_fit1_sigma0_18.81mb/Pb_smallest_b_";
  //const string bksolpp = "./kcbk_fit1_sigma0_18.81mb/p_smallest";

  // KCBK first  bal+sd solution
  //const string bksolpp = "./dipole-kcbk-hera-bal+sd-4.61.dip";
 
  // KCBK fit3 solution
  const string bksolpp = "./KCBK_fit_3/proton.dat";
  const string bksolpA = "./KCBK_fit_3/Pb_b_";
 
  // KCBK bal+sd solution
  //const string bksolpp = "./KCBK_fit_5/proton.dat";
  //const string bksolpA = "./KCBK_fit_5/Pb_b_";

  // KCBK parent solution
  //const string bksolpp = "./KCBK_fit_1/proton.dat";
  //const string bksolpA = "./KCBK_fit_1/Pb_b_";

  // ResumBK bal+sd solution
  //const string bksolpp = "./ResumBK_fit_5/proton.dat";
  //const string bksolpA = "./ResumBK_fit_5/Pb_b_";

  // ResumBK parent solution
  //const string bksolpp = "./ResumBK_fit_1/proton";
  //const string bksolpA = "./ResumBK_fit_1/Pb_b_";

  // ResumBK first parent solution
  //const string bksolpp = "./dipole-resumbk-hera-parent-4.61.dip";

  // TBK parent solution
  //const string bksolpp = "./TBK_fit_1/proton.dat";
  //const string bksolpA = "./TBK_fit_1/Pb_b_";

  //const string bksol="/Users/tawabyx/Documents/GitHub/nlodisfit/data/dipole_resumbk_sdrc.dip";
  //const string bksol = "/Users/tawabyx/Documents/code_NLO_sinc/Au/pA_impactb0.000000";
  //const string bksol="bksol/Sgrid_31.dat";
  //AmplitudeLib N(bksol);

  // initial condition
  //const double Qs02=0.2;
  //const double Qs02=0.0964;  //First Pb run and ResumBK parent solution
  //const double Qs02=0.0833;      // KCBK fit1 and parent solution
  const double Qs02=0.0680;      // KCBK fit3
  //const double Qs02=0.0905;	 // KCBK bal+sd solution
  //const double Qs02=0.0950;      // ResumBK bal+sd solution
  //const double Qs02=0.0917;     // TBK parent (fit1) solution
  const double ec=1.;
  // running coupling
  //const running_types alpha_s_running=FIXED;
  const double alpha_s_fixed=0.2*M_PI/3.; // (alpha_bar=0.2)
  const double alpha_s_freeze=0.7;
  // choice of terms
  bool with_Nc = false;
  bool with_CF = false;
  bool with_gl = false;
  bool with_gq = false;
  bool with_gg = false;
  /*const bool with_Nc = true;
  const bool with_CF = true;*/
  const bool with_xi1 = false;  // True iff one wants to use the subtracted scheme (still not exactly CXY, see the paper)
  /*const bool with_gl = false;    // Include q->g terms
  const bool with_gq = false;   // Include g->q terms
  const bool with_gg = false;   // Include g->g terms */
  // kinematics
  //const double SQRTS=500;    // Forward RHIC (GeV)
  //const double SQRTS=5020;     // LHC (GeV)
  //const double SQRTS=5000;     // LHCb pp (GeV)
  const double SQRTS=8160;     // New LHCb (GeV)
  const double yh=3;
  // PDFs
  const string pdfname="MSTW2008nlo90cl";
  double mu2;
  // constants
  double minlnr, maxlnr;
  const int Nc=3, Nf=3;

  // const double CF=4./3.;
  const double CF=3./2.;     // To be consistent with the large-Nc limit taken elsewhere

  const double beta0=(11.*Nc-2.*Nf)/3.;
  const double c0=2*exp(-M_EULER);
  const double LambdaQCD=0.241;
  const double alpha_s_mu_0=exp((2*M_PI)/(beta0*alpha_s_freeze));

  //const int Anucleus = 197;   // Henri's Au data
  const int Anucleus = 208;   // Heikki's Pb data

  const double RA = (1.12 * pow(Anucleus, 1/3)) + (0.86 * pow(Anucleus, -1/3));
  const double WSd = 0.54;
  const double sigma_inel = 179.7733;
  const double gamm = 1.21;
  //const double sigma0 = 50.2628683108;   // Henri's Au data
  //const double sigma0 = 96.6153078;        // ResumBK and old KCBK data
  const double sigma0 = 94.4580282;          // KCBK fit3
  //const double sigma0 = 44.5837784;      // KCBK bal+sd
  //const double sigma0 = 50.0283412;        // KCBK parent
  //const double sigma0 = 39.3446708;        // ResumBK parent
  //const double sigma0 = 40.320583;      // ResumBK bal+sd
  //const double sigma0 = 31.7941922;        // TBK parent (fit1)
  double TA;

  const string TAfile = "./TAvalues_Pb_Heikki.dat";     // Pb

  // intde parameters
  const double epsrel_intde=1e-4;
  // gsl parameters
  const int gsl_maxpoints=1000;
  const double epsabs_gsl=1e-8, epsrel_gsl=1e-6;

  // input parameters
    string col;
    double b;
    double p;
    string incoming;
    string outgoing;
    running_types alpha_s_running;
}

using namespace params;

#include "utils.cpp"

// #include "bksol.cpp"
#include "bksol_nlodisfit.cpp"

#include "Sr_interp_2D.cpp"

#include "Sr_interp_1D.cpp"

#include "func.cpp"

#include "xi_interp.cpp"

#include "xi_int.cpp"

#include "init_interp.cpp"

#include "intde1.cpp"

#include "sigma_NLO.cpp"

#include "sigma_LO.cpp"

#include "init.cpp"


//extern "C"
//{
    // ih: input hadron, 1=pion, 2=kaon, 3=proton, 4=charged hadrons
    // ic: charge, 0=0, 1=+, -1=-
    // parton:
        //0    1    2    3    4    5    6    7    8     9    10
        //g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
    // io: order, 0=LO, 1=NLO
    // NOTE: cbar and bbar not in DSS, returns c or b instead
    // result is the fragmentation function (X factor divided out)
    
    /*void fdss_ (int &ih, int &ic, int &io, double &x, double& q2, double& u,
        double &ub, double &d, double &db, double &s, double &sb, double &c,
        double &b, double &g);*/
/*    void fdss_(int& hadron, int& charge, int& order, double& z,
              double& scalesqr, int& parton, double& result);

    extern struct{
        double fini;
    } fragini_;
}*/

int main(int argc, char* argv[]){
    
    using namespace params;
    
    double zmin = str_to_double(argv[1]);
    double zmax = str_to_double(argv[2]);
    double zstep = str_to_double(argv[3]);
    col = argv[4];
    b = str_to_double(argv[5]);
    incoming = argv[6];
    outgoing = argv[7];
    string rc = argv[8];
    p = str_to_double(argv[9]);
    double muratio = str_to_double(argv[10]);
    mu2 = Sq(muratio*p);
    if(rc.compare("fixed") == 0) alpha_s_running = FIXED;
    else if(rc.compare("mom") == 0) alpha_s_running = MOM;
    else if(rc.compare("parent") == 0) alpha_s_running = PARENT;
    else if(rc.compare("daughter") == 0) alpha_s_running = DAUGHTER;
    else if(rc.compare("mixed") == 0) alpha_s_running = MIXED;
    else if(rc.compare("mixedbd") == 0) alpha_s_running = MIXEDBD;
    else if(rc.compare("smallest") == 0) alpha_s_running = SMALLEST;
    else cout << "Invalid running coupling choice." << endl;
    
    ifstream datafile(TAfile.c_str());
    if(!datafile.is_open()){
        cerr << "Error opening the TA file" << endl;
        exit(1);
    }
    while(!datafile.eof()){
        string line;
        getline(datafile, line);
        if(line.size() == 0){
            break;
        }
        int loc = 0;
        while(line.compare(loc,1,",") != 0){
            loc++;
        }
        double this_b = str_to_double(line.substr(0, loc));
        if(this_b == b){
            TA = str_to_double(line.substr(loc+1, line.size()-loc-1));
            break;
        }
    }
    datafile.close();
    
    if(incoming.compare("g") == 0){
        if(outgoing.compare("g") == 0) with_gg = true;
        else with_gq = true;
    }
    else{
        if(outgoing.compare("g") == 0) with_gl = true;
        else{
            with_Nc = true;
            with_CF = true;
        }
    }

    init();

    double z = zmax;
    double k = p/zmax;
    //double k = 20.1;
    while(z >= zmin - 0.00001){
    //while(k <= 50.00001){

      //mu2 = Sq(4*k);

      k = p/z;
      double xp=(k/SQRTS)*exp(yh);
      double xg=(k/SQRTS)*exp(-yh);
      init_interp(xp,xg,k);
      cout << z << "," << k << "," << sigma_LO_k(k,xp) << "," << sigma_NLO_k(k,xp) << endl;
      //cout << k << "," << sigma_LO_k(k,xp) << "," << sigma_NLO_k(k,xp) << endl;
      clear_interp();
      z-=zstep;
      //k += 0.1;
   }
        
    return 0;
}

