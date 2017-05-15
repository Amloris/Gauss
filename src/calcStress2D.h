/* calcStress2D
The calcStress2D class stores stress data in 2D. It has sigma_xx,
sigma_yy, sigma_xy, and sigma_zz components.
*/

#ifndef calcStress2D_h
#define calcStress2D_h

class calcStress2D {
private:
	// stress data
	double m_sigxx, m_sigyy, m_sigzz, m_sigxy;	//stress components
												//note that four components are stored

public:
	// constructors
	calcStress2D();
	calcStress2D(double xx, double yy, double zz, double xy);
	// destructor
	~calcStress2D();
	//functions
	void setSigXX(double xx);	//	set m_sigxx = xx;
	void setSigYY(double yy);	//	set m_sigyy = yy;
	void setSigZZ(double zz);	//	set m_sigzz = zz;
	void setSigXY(double xy);	//	set m_sigxy = xy;
	void addSigXX(double xx);	//	add xx to m_sigxx;
	void addSigYY(double yy);	//	add yy to m_sigyy;
	void addSigZZ(double zz);	//	add zz to m_sigzz;
	void addSigXY(double xy);	//	add xy to m_sigxy;
	double getSigXX() const;	//	return m_sigxx;
	double getSigYY() const;	//	return m_sigyy;
	double getSigZZ() const;	//	return m_sigzz;
	double getSigXY() const;	//	return m_sigxy;
	void zero();		//initialize all components to zero
};

#endif