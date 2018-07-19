#include <HGCal/Reco/interface/PulseFitter.h>


#include <limits>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

//seems to be mandatory since we need static function
double _time[11],_energy[11];
double _maxTime=200.; 
double _trise=35.;
double _noise=5.;
double _ampl_norm = 1.85; // amplitude normalization factor for tau = 18, n = 3
double _tau = 20.;
int _n_ord = 3;

//new shape motivated here: H. Spieler - Semiconductor Detector Systems, pg. 179 for a CR2-RC
//implemented by A. Lobanov, suggested by F. Pitters, 06 April 2018
double pulseShape_fcn(double t, double tmax, double amp, double amp0 = 0., double tau = 19., int n_ord = 3){     //tau = 18 instead of 19 [ns] mattermost discussion on 13 April 2018

    const double ampl_norm = 1.85; // amplitude normalization factor for tau = 19, n = 3

    if( t>tmax-_trise )
	return (amp*ampl_norm * (1 - ((t-(tmax-_trise))/tau)/(n_ord+1)) * std::pow((t-(tmax-_trise))/tau, n_ord) * std::exp(-(t-(tmax-_trise))/tau)) + amp0;
    else return 0;

}

double pulseShape_chi2(const double *x)
{
  double sum = 0.0;

  for(size_t i=0; i<9; i++){
    if( _energy[i]<-170 || _time[i]>_maxTime ) continue;
    double zero = _energy[i]-pulseShape_fcn( _time[i],
	       x[0],x[1] );
    sum += zero * zero / _noise / _noise;
  }
  return sum;
}


PulseFitter::PulseFitter( int printLevel, double maxTime , double trise , double ampl_norm , double tau , int n_ord ) : m_printLevel(printLevel)
{
  _maxTime=maxTime;
  _trise=trise;
  _ampl_norm = ampl_norm;
  _tau = tau;
  _n_ord = n_ord;  
}

void PulseFitter::run(std::vector<Float16_t> &time, std::vector<Float16_t> &energy, PulseFitterResult &fit, Float16_t noise)
{
  if( time.size()!=energy.size() ){
    std::cout << "ERROR : we should have the same vector size in PulseFitter::run(std::vector<double> time, std::vector<double> energy, PulseFitterResult fit) -> return without fitting" << std::endl;
    return;
  }
  if( time.size()!=11 ){
    std::cout << "ERROR : we should have less than 13 time sample in PulseFitter::run(std::vector<double> time, std::vector<double> energy, PulseFitterResult fit) -> return without fitting" << std::endl;
    return;
  }

  float emax0(0),tmax0(0);
  for( uint16_t i=0; i<9; i++ ){
      _time[i] = time[i];
      _energy[i] = energy[i];
      if(_energy[i]>emax0) {emax0=energy[i];tmax0=time[i];}
  }

  if( noise>0 )
    _noise=noise;

  ROOT::Math::Minimizer* m = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
  m->SetMaxFunctionCalls(m_fitterParameter.nMaxIterations);
  m->SetMaxIterations(m_fitterParameter.nMaxIterations);
  m->SetTolerance(0.1);
  m->SetPrintLevel(m_printLevel);

  ROOT::Math::Functor f(&pulseShape_chi2, 2);

  m->SetFunction(f);

  m->Clear(); // just a precaution

  m->SetVariable(0, "tmax", tmax0, 0.5);
  m->SetVariableLimits(0,m_fitterParameter.tmaxRangeDown,m_fitterParameter.tmaxRangeUp);
  m->SetVariable(1, "amp", emax0, 1);

  m->SetVariableLimits(1,0,10000);

  m->Minimize();

  const double *xm = m->X();
  const double *errors = m->Errors();
  fit.tmax=xm[0];
  fit.trise=_trise;
  fit.amplitude=xm[1];
  fit.errortmax=errors[0];
  fit.erroramplitude=errors[1];
  fit.chi2=m->MinValue();
  fit.status=(fabs(xm[0]-m_fitterParameter.tmaxRangeDown)>std::numeric_limits<double>::epsilon()&&
	      fabs(xm[0]-m_fitterParameter.tmaxRangeUp)>std::numeric_limits<double>::epsilon()) ? m->Status() : 6;

  fit.ncalls=m->NCalls();

  delete m;
}

