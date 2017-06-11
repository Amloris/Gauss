/* calcStress2D
The calcStress2D class stores stress data in 2D. It has sigma_xx,
sigma_yy, sigma_xy, and sigma_zz components.
*/

#ifndef calcStress2D_h
#define calcStress2D_h

class calcStress2D {
private:
	// stress data
	double m_sigxx, m_sigyy, m_sigzz, m_sigxy;	 //stress components, note that four components are stored

public:
	// constructors
	calcStress2D();
	calcStress2D(double xx, double yy, double zz, double xy);
	// destructor
	~calcStress2D();
	//functions
	void setSigXX(double xx);	  //Set m_sigxx = xx;
	void setSigYY(double yy);	  //Set m_sigyy = yy;
	void setSigZZ(double zz);	  //Set m_sigzz = zz;
	void setSigXY(double xy);	  //set m_sigxy = xy;
	void addSigXX(double xx);	  //Add xx to m_sigxx;
	void addSigYY(double yy);	  //Add yy to m_sigyy;
	void addSigZZ(double zz);	  //Add zz to m_sigzz;
	void addSigXY(double xy);	  //Add xy to m_sigxy;
	double getSigXX() const;	  //Return m_sigxx;
	double getSigYY() const;	  //Return m_sigyy;
	double getSigZZ() const;	  //Return m_sigzz;
	double getSigXY() const;	  //Return m_sigxy;
	void zero();		          //Initialize all components to zero
};


//Functions
//-----------------------------------------------------------------------------
calcStress2D::calcStress2D()  { /*std::cout << "Creating calcStress2D Class Object" << endl;*/ }
calcStress2D::~calcStress2D() { /*std::cout << "Deleting calcStress2D Class Object" << endl;*/ }
calcStress2D::calcStress2D(double xx, double yy, double zz, double xy)
{
	m_sigxx = xx;
	m_sigyy = yy;
	m_sigzz = zz;
	m_sigxy = xy;
}

void calcStress2D::setSigXX(double xx)
{	
	m_sigxx = xx;                 //Set m_sigxx = xx
}

void calcStress2D::setSigYY(double yy)
{
	m_sigyy = yy;                 //Set m_sigyy = yy
}

void calcStress2D::setSigZZ(double zz)
{
	m_sigzz = zz;                 //Set m_sigzz = zz
}

void calcStress2D::setSigXY(double xy)	  
{ 
	m_sigxy = xy;                 //set m_sigxy = xy
}

void calcStress2D::addSigXX(double xx)	  
{
	m_sigxx += xx;                //Add xx to m_sigxx
}

void calcStress2D::addSigYY(double yy)	  
{
	m_sigyy += yy;                //Add yy to m_sigyy
}

void calcStress2D::addSigZZ(double zz)	  
{
	m_sigzz += zz;                //Add zz to m_sigzz
}

void calcStress2D::addSigXY(double xy)	  
{
	m_sigxy += xy;                //Add xy to m_sigxy
}

double calcStress2D::getSigXX() const	  
{
	return m_sigxx;               //Return m_sigxx
}

double calcStress2D::getSigYY() const	  
{
	return m_sigyy;               //Return m_sigyy
}

double calcStress2D::getSigZZ() const	  
{
	return m_sigzz;               //Return m_sigzz
}

double calcStress2D::getSigXY() const	  
{
	return m_sigxy;               //Return m_sigxy
}

void calcStress2D::zero()		  
{
	m_sigxx = 0.0;                //Initialize all components to zero
	m_sigyy = 0.0;
	m_sigzz = 0.0;
	m_sigxy = 0.0;
}

#endif