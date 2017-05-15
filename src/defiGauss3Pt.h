/* defiGauss3Pt
This class stores the gauss point data for three point integration.
*/

#ifndef defiGauss3Pt_h
#define defiGauss3Pt_h


//using namespace std;

class defiGauss3Pt
{
private:
	int d_numgp;		//number of Gauss integration points, d_numgp = 3;
	double d_point[3];	//location of Gauss point, i.e. d_point[0] = -0.774596669241483 etc.
	double d_weight[3];	//weight, i e. d_weight[0] = 0.555555555555556 etc.

public:
	//Constructor
	defiGauss3Pt();
	//Destructor
	~defiGauss3Pt();

	// Functions
	int getNumPts() const;				//return d_numgp;
	double getGaussPt(int gp) const;	//return d_point[gp];
	double getWeight(int gp) const;		//return d_weight[gp];
};


//Functions
//-----------------------------------------------------------------------------
defiGauss3Pt::~defiGauss3Pt() {}
defiGauss3Pt::defiGauss3Pt()
{
	d_numgp = 3;
	d_point[0] = -0.774596669241483;
	d_point[1] = 0.0;
	d_point[2] = 0.774596669241483;
	d_weight[0] = 0.555555555555556;
	d_weight[1] = 0.888888888888889;
	d_weight[2] = 0.555555555555556;
}

int defiGauss3Pt::getNumPts() const
{
	return d_numgp;               //Return d_numgp;
}

double defiGauss3Pt::getGaussPt(int gp) const
{
	return d_point[gp];           //Return d_point[gp]
}

double defiGauss3Pt::getWeight(int gp) const
{
	return d_weight[gp];          //Return d_weight[gp]
}


#endif