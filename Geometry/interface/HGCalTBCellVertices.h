#ifndef HGCAL_GEOMETRY_HGCALTBCELLVERTICES_H
#define HGCAL_GEOMETRY_HGCALTBCELLVERTICES_H

#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "math.h"
#include <vector>

class HGCalTBCellVertices
{
public:
	/** Constructor from cell ix & iv, valid sensorSizes are 128 and 256 */
	HGCalTBCellVertices();
	std::vector<std::pair<double, double>> GetCellCoordinates(int layer, int sensor_ix, int sensor_iv, int ix, int iv, int sensorsize);
	std::pair<double, double> GetCellCentreCoordinates(int layer, int sensor_ix, int sensor_iv, int ix, int iv, int sensorsize);
//  void CellType(int ix, int iv, bool ValidFlag);// 1 for full hex, 2 for half hex and 3 for the pentagons(to be implemented later)
private:
	double a = 1; // Size in terms of 1 unit of x/y co-ordinates of a cell side
	double x_a = sqrt(3) / 2; // cosine pi/6
	double y_a = 1 / 2.; // sine pi/6
	double vy_a = 3. / 2;

	std::vector<double> x_co_FullHex, y_co_FullHex; // stores the initial x,y coordinates of a hexagonal cell
	std::vector<std::pair<double, double>> Cell_co;
// Translation in x,v co-ordinates in terms of cartesian x,y.
	double  x0 = 2 * x_a * a; //Translation in Cartesian x for 1 unit of ix
	double vx0 = x_a * a; // Cartesian x component of translation for 1 unit of iv
	double vy0 = vy_a * a; // Cartesian y component of translation for 1 unit of iv

	double Xmax(int iv, double y);// returns the max x value for a cell to be in the given sensor

};

//std::ostream& operator<<(std::ostream&,const HGCalTBCellVertices& Cell_Ix_Iv);

#endif

