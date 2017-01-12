#ifndef HGCAL_RECO_PARTICLETRACK_H
#define HGCAL_RECO_PARTICLETRACK_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h" //e.g. to get the cell's dimensions

#include "FWCore/Utilities/interface/Exception.h"


enum ConsiderationMethod {
  CONSIDERALL,
  CONSIDERSEVEN,
  CONSIDERNINETEEN,
  CONSIDERCLUSTERSALL,
  CONSIDERCLUSTERSSEVEN,
  CONSIDERCLUSTERSNINETEEN
};


enum WeightingMethod {
  DEFAULTWEIGHTING,
  SQUAREDWEIGHTING,
  LINEARWEIGHTING,
  LOGWEIGHTING_30_10,
  LOGWEIGHTING_30_15,
  LOGWEIGHTING_40_10,
  LOGWEIGHTING_40_15,
  LOGWEIGHTING_50_10,
  LOGWEIGHTING_50_15,
  LOGWEIGHTING_60_10,
  LOGWEIGHTING_60_15,
  LOGWEIGHTING_70_10,
  LOGWEIGHTING_70_15
};

enum TrackFittingMethod {
  DEFAULTFITTING,
  LINEFITANALYTICAL,
  LINEFITTGRAPHERRORS,
  POL2TGRAPHERRORS,
  POL3TGRAPHERRORS
};

enum FitPointWeightingMethod {
  NONE, 
  LINEAR,
  SQUARED,
  LOGARITHMIC,
  EXPONENTIAL
};

struct HitData {
  double x; double y; 
  double I; //intensity that is input to the weight calculation
  double E; //actual energy of the hit
  int ID;   //the ID corresponds to the cell ID, it is necessary for the pedestal subtraction
};

void parseAlignmentFile(std::map<int, double> &alignmentParameters, std::string path);

//
// class declarations
//

//class that performs an analytical straight line fit
class LineFitter {
  private:
    std::vector<double> _x;
    std::vector<double> _y;
    std::vector<double> _sigma_y;
    double _S, _S_x, _S_xx, _S_xy, _S_y;
  public:
    LineFitter(std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y);
    void addPoint(double x, double y, double sigma_y);
    void fit();
    bool converged();
    double getM();
    double getMError();
    double getB();
    double getBError();
    double getMBCovariance();

    double eval(double x);
    double evalError(double x);
};



class SensorHitMap {
  private:
    std::pair<double, double> centralHitPoint;
    double centralHitZ;
    std::pair<double, double> centralHitPointError;
    double layerLabZ;
    double layerZ_X0;
    double d_alpha, d_beta, d_gamma, d_x0, d_y0, d_z0; //alignment parameters
    int sensorSize;
    double CM_threshold;
    double ADC_per_MIP;
    std::map<int, HitData*> Hits; //those are all the hits in the layer
    std::map<int, std::vector<int>> clusterIndexes;
    std::vector<HitData*> HitsForPositioning;
    HitData* mostSignificantHit;

    //helpers to obtain the x-y coordinate
    HGCalTBCellVertices TheCell;
    std::pair<double, double> CellCenterXY;
    int CM_cells_count;
    double CM_sum;

    double totalWeight; //equivalent to the denominator in the according weighting method
    double totalEnergy;
    std::map<int, double> totalClusterEnergy;

    bool filterByCellType(int ID);
    void considerNClosest(int N_considered);
    void considerClusters(int N_considered);
    void poweredWeighting(int exponent);
    void logWeighting(double log_a, double log_b);

  public:
    SensorHitMap();
    ~SensorHitMap();
    void setLabZ(double z_cm, double z_X0);
    void setAlignmentParameters(double d_alpha, double d_beta, double d_gamma, double d_x0, double d_y0, double d_z0);
    void setADCPerMIP(double ADC_per_MIP);
    void setSensorSize(int s);
    void setPedestalThreshold(double t);
    //reduces the information from the Rechit towards what is necessary for the impact point calculation
    void addHit(HGCalTBRecHit Rechit);
    void registerClusterHit(HGCalTBDetId hit, int N_considered);
    std::pair<int, double> subtractCM();  //returns the sum of Common mode noise and the number of cells that enter the calculation
    void calculateCenterPosition(ConsiderationMethod considerationMethod, WeightingMethod weightingMethod);
    double getTotalEnergy();
    double getTotalClusterEnergy(int N_considered);
    double getTotalWeight();
    double getLabZ();
    double getZ_X0();
    double getIntrinsicHitZPosition();
    std::pair<double, double> getHitPosition(); //returns central hit in layer's own frame
    std::pair<double, double> getLabHitPosition();  //returns central hit in lab frame
    std::pair<double, double> getHitPositionError(); //calculated via RMS
    std::pair<double, double> getCenterOfClosestCell(std::pair<double, double> X_ref);

    //debug
    void printHits();
};


class ParticleTrack{
  public:
    ParticleTrack();
    ~ParticleTrack();
    void addFitPoint(SensorHitMap* sensor);
    void weightFitPoints(FitPointWeightingMethod method);
    void fitTrack(TrackFittingMethod method);
    std::pair<double, double> calculatePositionXY(double z);
    std::pair<double, double> calculatePositionErrorXY(double z);
    double getSumOfEnergies();
  private:
    int N_points; //cross check, should be identical to x.size() etc.
    //general information, the fit points
    std::vector<double> x;  
    std::vector<double> x_err;  
    std::vector<double> y;  
    std::vector<double> y_err;  
    std::vector<double> z;  
    std::vector<double> z_err;  
    std::vector<double> Energies;
    TrackFittingMethod lastAppliedMethod;

    //different fit functions
    void analyticalStraightLineFit();
    std::pair<double, double> positionFromAnalyticalStraightLine(double z);
    std::pair<double, double> positionFromAnalyticalStraightLineErrors(double z);
    LineFitter* linefit_x;
    LineFitter* linefit_y;

    void polFitTGraphErrors(int degree);  
    std::pair<double, double> positionFromPolFitTGraphErrors(int degree, double z);
    std::pair<double, double> positionErrorFromPolFitTGraphErrors(int degree, double z);
    TF1* ROOTpol_x;
    TF1* ROOTpol_y;
    TFitResultPtr fit_result_x;
    TFitResultPtr fit_result_y;
};

double weightToFitPointWeight(double w, double sum_w, FitPointWeightingMethod m);
#endif
