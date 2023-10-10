void init(){
  LHAPDF::setVerbosity(0);
  pdf::lhpdf=LHAPDF::mkPDF(pdfname,0);

  init_bksol();
  init_Sr_interp_2D();

  gsl_set_error_handler(&gsl_error_handler);
}
