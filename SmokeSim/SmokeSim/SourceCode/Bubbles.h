#ifndef Bubbles_H_
#define Bubbles_H_

#pragma warning(disable: 4244 4267 4996)

#include <vector>
#include "vec.h"

class Bubbles
{
public:
	// default constructor
	Bubbles( void );

	// copy constructor
	Bubbles( const Bubbles& orig );

	// destructor
	~Bubbles( void );
	
	// assignment operator
	virtual Bubbles& operator=( const Bubbles& orig );

	// get bubble positions
	//std::vector<vec3>& m_positions( void );

private:
	std::vector<vec3> m_positions;
};

#endif