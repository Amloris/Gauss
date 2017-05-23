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

#endif