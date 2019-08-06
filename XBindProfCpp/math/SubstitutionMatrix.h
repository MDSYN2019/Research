#ifndef _SubstitutionMatrix_
#define _SubstitutionMatrix_

class SubstitutionMatrix
{
protected:
	double matrix[20][20];

public:
	SubstitutionMatrix();
	SubstitutionMatrix& operator=(const SubstitutionMatrix& other);
	~SubstitutionMatrix();
	double getValue(int i, int j);

};

#endif