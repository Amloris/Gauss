/* defiDof
The defiDof class implements a degree of freedom.  Each dof has an equation
number with associated value.  Dof's with either no boundary condition or
with a natural boundary condition are active.  If the dof has an essential
boundary condition, the dof is not active.  During matrix assembly, if the
boundary condition is essential, the appropriate element stiffness
coefficient is multiplied by the boundary condition value and subtracted
from the global load vector.  If it is a natural boundary condition the
value is added to the load vector.

After the solution, the results (displacements) are stored in the variable
d_value.
*/

#ifndef defiDof_h
#define defiDof_h

using namespace std;

enum DOFType { UX = 0, UY };

class defiDof
{
private:
	double d_value;	              //The value of this degree of freedom (either UX or UY). For essential BC (dof is inactive), it is the prescribed displacement inputted by the user. For active dof, it stores the displacement value after the system equations are solved.
	int d_eqn;		              //The equation number of this dof. d_eqn has meaning only for active dof.
	bool d_active;	              //d_active = true for active dof, = false for inactive dof
					
	defiDof(defiDof &dof);        //Never used copy constructor

public:
	// Constructors
	defiDof();
	//Destructor
	~defiDof();

	// Functions
	void setNotActive();	      //Set d_active = false;
	bool isActive() const;	      //Return d_active;
	void setEqn(int count);	      //For active dof, set d_eqn = count; for inactive dof, set d_eqn = -1 (or any negative integer)
	int getEqn() const;		      //Return d_eqn;
	void setValue(double value);  //Set d_value = value;
	double getValue() const;	  //Return d_value;
	void print() const;	          //Print any useful info
};


//Functions
//-----------------------------------------------------------------------------
//defiDof::defiDof()  { std::cout << "Creating defiDof Class Object" << endl; }
defiDof::~defiDof() { std::cout << "Deleting defiDof Class Object" << endl; }
defiDof::defiDof()
{
	d_active = true;              //Assume the dof is active and unconstrained
	d_eqn = -1;                   //Assign a eqn number placeholder
	d_value = 0.0;                //Set computed displacement to zero
}

void defiDof::setNotActive()
{
	d_active = false;             //Deactivate the dof (dof of point is an essentialBC)
}

bool defiDof::isActive() const
{
	return d_active;              //Returns the status of the dof
}

void defiDof::setEqn(int count)
{
	if (d_active) 
	{ 
		d_eqn = count;            //Assigns eqn number if active
	}   
	else 
	{ 
		d_eqn = -1;               //Sets eqn number to negative if inactive
	}
}

int defiDof::getEqn() const
{
	return d_eqn;                 //Return the equation number of the dof
}

void defiDof::setValue(double value)
{
	d_value = value;              //Sets the dispacment value
}

double defiDof::getValue() const
{
	return d_value;               //Returns the displacement value
}

#endif