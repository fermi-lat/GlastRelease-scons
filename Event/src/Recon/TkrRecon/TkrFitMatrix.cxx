//#define TIMETRIAL 1
#undef TIMETRIAL

#ifdef TIMETRIAL
#include "TkrFitMatrix.h"
#include "TkrFitPar.h"
#else /* Standard Gleam */
#include "Event/Recon/TkrRecon/TkrFitMatrix.h"
#include "Event/Recon/TkrRecon/TkrFitPar.h"
#endif /* TIMETRIAL */
using namespace Event;
//#include <iostream>

TkrFitMatrix::TkrFitMatrix() 
{

	m_11 = 0.0;
	m_12 = 0.0;
	m_13 = 0.0;
	m_14 = 0.0;

	m_21 = 0.0;
	m_22 = 0.0;
	m_23 = 0.0;
	m_24 = 0.0;

	m_31 = 0.0;
	m_32 = 0.0;
	m_33 = 0.0;
	m_34 = 0.0;


	m_41 = 0.0;
	m_42 = 0.0;
	m_43 = 0.0;
	m_44 = 0.0;


	m_zero = 0.0;  
}
// A matrix
TkrFitMatrix::TkrFitMatrix(
													 double e_11, double e_12, double e_13, double e_14,
													 double e_21, double e_22, double e_23, double e_24,
													 double e_31, double e_32, double e_33, double e_34,
													 double e_41, double e_42, double e_43, double e_44
													 )
{

	m_11 = e_11;
	m_12 = e_12;
	m_13 = e_13;
	m_14 = e_14;

	m_21 = e_21;
	m_22 = e_22;
	m_23 = e_23;
	m_24 = e_24;

	m_31 = e_31;
	m_32 = e_32;
	m_33 = e_33;
	m_34 = e_34;


	m_41 = e_41;
	m_42 = e_42;
	m_43 = e_43;
	m_44 = e_44;


	m_zero = 0.0;  

}

// Init from a HepMatrix
TkrFitMatrix::TkrFitMatrix(const HepMatrix &A) 
{

	m_11 = A(1,1);
	m_12 = A(1,2);
	m_13 = A(1,3);
	m_14 = A(1,4);

	m_21 = A(2,1);
	m_22 = A(2,2);
	m_23 = A(2,3);
	m_24 = A(2,4);

	m_31 = A(3,1);
	m_32 = A(3,2);
	m_33 = A(3,3);
	m_34 = A(3,4);

	m_41 = A(4,1);
	m_42 = A(4,2);
	m_43 = A(4,3);
	m_44 = A(4,4);


	m_zero = 0.0;  

}

// Init from a TkrFitMatrix
TkrFitMatrix::TkrFitMatrix(const TkrFitMatrix &A) 
{

	m_11 = A.m_11;
	m_12 = A.m_12;
	m_13 = A.m_13;
	m_14 = A.m_14;

	m_21 = A.m_21;
	m_22 = A.m_22;
	m_23 = A.m_23;
	m_24 = A.m_24;

	m_31 = A.m_31;
	m_32 = A.m_32;
	m_33 = A.m_33;
	m_34 = A.m_34;

	m_41 = A.m_41;
	m_42 = A.m_42;
	m_43 = A.m_43;
	m_44 = A.m_44;


	m_zero = 0.0;  

}


// Special constructor to produce propagation matrix F
TkrFitMatrix::TkrFitMatrix(double step) 
{

	m_11 = 1.0;
	m_12 = step;
	m_13 = 0.0;
	m_14 = 0.0;

	m_21 = 0.0;
	m_22 = 1.0;
	m_23 = 0.0;
	m_24 = 0.0;

	m_31 = 0.0;
	m_32 = 0.0;
	m_33 = 1.0;
	m_34 = step;

	m_41 = 0.0;
	m_42 = 0.0;
	m_43 = 0.0;
	m_44 = 1.0;


	m_zero = 0.0;
}

// Special constructor to  produce matrix H
TkrFitMatrix::TkrFitMatrix(int) 
{

	m_11 = 1.0;
	m_12 = 0.0;
	m_13 = 0.0;
	m_14 = 0.0;

	m_21 = 0.0;
	m_22 = 0.0;
	m_23 = 0.0;
	m_24 = 0.0;

	m_31 = 0.0;
	m_32 = 0.0;
	m_33 = 1.0;
	m_34 = 0.0;


	m_41 = 0.0;
	m_42 = 0.0;
	m_43 = 0.0;
	m_44 = 0.0;


	m_zero = 0.0;

}


//  function for calculating the transpose of TkrFitMatrix
TkrFitMatrix TkrFitMatrix::T() const 
{

	return TkrFitMatrix(
		m_11,m_21,m_31,m_41,
		m_12,m_22,m_32,m_42,
		m_13,m_23,m_33,m_43,
		m_14,m_24,m_34,m_44);

}

//*******************************************************
// TESTING FASTER(?) Invert method
//*******************************************************

void TkrFitMatrix::invert(int &ierr) 
{


	ierr = 0;

	// Find all NECESSARY 2x2 dets:  (18 of them)

	double Det2_12_01 = m_21*m_32 - m_22*m_31;
	double Det2_12_02 = m_21*m_33 - m_23*m_31;
	double Det2_12_03 = m_21*m_34 - m_24*m_31;			//
	double Det2_12_13 = m_22*m_34 - m_24*m_32;			//
	double Det2_12_23 = m_23*m_34 - m_24*m_33;			//
	double Det2_12_12 = m_22*m_33 - m_23*m_32;
	double Det2_13_01 = m_21*m_42 - m_22*m_41;
	double Det2_13_02 = m_21*m_43 - m_23*m_41;
	double Det2_13_03 = m_21*m_44 - m_24*m_41;
	double Det2_13_12 = m_22*m_43 - m_23*m_42;
	double Det2_13_13 = m_22*m_44 - m_24*m_42;
	double Det2_13_23 = m_23*m_44 - m_24*m_43;			//
	double Det2_23_01 = m_31*m_42 - m_32*m_41;
	double Det2_23_02 = m_31*m_43 - m_33*m_41;
	double Det2_23_03 = m_31*m_44 - m_34*m_41;
	double Det2_23_12 = m_32*m_43 - m_33*m_42;
	double Det2_23_13 = m_32*m_44 - m_34*m_42;
	double Det2_23_23 = m_33*m_44 - m_34*m_43;

	// Find all NECESSARY 3x3 dets:   (16 of them)

	double Det3_012_012 = m_11*Det2_12_12 - m_12*Det2_12_02 
		+ m_13*Det2_12_01;
	double Det3_012_013 = m_11*Det2_12_13 - m_12*Det2_12_03 
		+ m_14*Det2_12_01;			//
	double Det3_012_023 = m_11*Det2_12_23 - m_13*Det2_12_03 
		+ m_14*Det2_12_02;			//
	double Det3_012_123 = m_12*Det2_12_23 - m_13*Det2_12_13 
		+ m_14*Det2_12_12;			//
	double Det3_013_012 = m_11*Det2_13_12 - m_12*Det2_13_02 
		+ m_13*Det2_13_01;
	double Det3_013_013 = m_11*Det2_13_13 - m_12*Det2_13_03
		+ m_14*Det2_13_01;
	double Det3_013_023 = m_11*Det2_13_23 - m_13*Det2_13_03
		+ m_14*Det2_13_02;			//
	double Det3_013_123 = m_12*Det2_13_23 - m_13*Det2_13_13
		+ m_14*Det2_13_12;			//
	double Det3_023_012 = m_11*Det2_23_12 - m_12*Det2_23_02 
		+ m_13*Det2_23_01;
	double Det3_023_013 = m_11*Det2_23_13 - m_12*Det2_23_03
		+ m_14*Det2_23_01;
	double Det3_023_023 = m_11*Det2_23_23 - m_13*Det2_23_03
		+ m_14*Det2_23_02;
	double Det3_023_123 = m_12*Det2_23_23 - m_13*Det2_23_13
		+ m_14*Det2_23_12;			//
	double Det3_123_012 = m_21*Det2_23_12 - m_22*Det2_23_02 
		+ m_23*Det2_23_01;
	double Det3_123_013 = m_21*Det2_23_13 - m_22*Det2_23_03 
		+ m_24*Det2_23_01;
	double Det3_123_023 = m_21*Det2_23_23 - m_23*Det2_23_03 
		+ m_24*Det2_23_02;
	double Det3_123_123 = m_22*Det2_23_23 - m_23*Det2_23_13 
		+ m_24*Det2_23_12;

	// Find the 4x4 det:

	double det =    m_11*Det3_123_123 
		- m_12*Det3_123_023 
		+ m_13*Det3_123_013 
		- m_14*Det3_123_012;

	if ( det == 0 ) {  
		ierr = 1;

	} else {

		double oneOverDet = 1.0/det;
		double mn1OverDet = - oneOverDet;

		m_11 =  Det3_123_123 * oneOverDet;
		m_12 =  Det3_023_123 * mn1OverDet;
		m_13 =  Det3_013_123 * oneOverDet;
		m_14 =  Det3_012_123 * mn1OverDet;

		m_21 =  Det3_123_023 * mn1OverDet;
		m_22 =  Det3_023_023 * oneOverDet;
		m_23 =  Det3_013_023 * mn1OverDet;
		m_24 =  Det3_012_023 * oneOverDet;

		m_31 =  Det3_123_013 * oneOverDet;
		m_32 =  Det3_023_013 * mn1OverDet;
		m_33 =  Det3_013_013 * oneOverDet;
		m_34 =  Det3_012_013 * mn1OverDet;

		m_41 =  Det3_123_012 * mn1OverDet;
		m_42 =  Det3_023_012 * oneOverDet;
		m_43 =  Det3_013_012 * mn1OverDet;
		m_44 =  Det3_012_012 * oneOverDet;

	}

}

//*******************************************************



//  function for calculating the inverse of TkrFitMatrix
TkrFitMatrix TkrFitMatrix::inverse(int &ierr) const 
{
	TkrFitMatrix temp(*this);
	temp.invert(ierr);
	return temp;

}

//operator overload for element access TkrFitMatrix(row,column)
//!NOTE index starts with (1,1)
//!WARNING DOES NOT GIVE ERROR FEEDBACK
const double & TkrFitMatrix::operator() 
(const int &row, const int &column) const
{

	if (row == 1) {
		if(column == 1) return m_11;
		else if(column == 2) return m_12;
		else if(column == 3) return m_13;
		else if(column == 4) return m_14;
	} else if (row == 2) {
		if(column == 1) return m_21;
		else if(column == 2) return m_22;
		else if(column == 3) return m_23;
		else if(column == 4) return m_24;
	} else if (row == 3) {
		if(column == 1) return m_31;
		else if(column == 2) return m_32;
		else if(column == 3) return m_33;
		else if(column == 4) return m_34;
	} else if (row == 4) {
		if(column == 1) return m_41;
		else if(column == 2) return m_42;
		else if(column == 3) return m_43;
		else if(column == 4) return m_44;
	} else return m_zero;

	//hack to get no error should not be here
	return m_zero;

}

//!WARNING DOES NOT GIVE ERROR FEEDBACK    
double & TkrFitMatrix::operator() (const int &row, const int &column)
{

	if (row == 1) {
		if(column == 1) return m_11;
		else if(column == 2) return m_12;
		else if(column == 3) return m_13;
		else if(column == 4) return m_14;
	} else if (row == 2) {
		if(column == 1) return m_21;
		else if(column == 2) return m_22;
		else if(column == 3) return m_23;
		else if(column == 4) return m_24;
	} else if (row == 3) {
		if(column == 1) return m_31;
		else if(column == 2) return m_32;
		else if(column == 3) return m_33;
		else if(column == 4) return m_34;
	} else if (row == 4) {
		if(column == 1) return m_41;
		else if(column == 2) return m_42;
		else if(column == 3) return m_43;
		else if(column == 4) return m_44;
	} else return m_zero;

	//hack to get no error should not be here
	return m_zero;
}


//operator overload for TkrFitMatrix+TkrFitMatrix
const TkrFitMatrix TkrFitMatrix::operator +(const TkrFitMatrix& A) const
{

	double e_11 = m_11 + A.m_11;
	double e_12 = m_12 + A.m_12;
	double e_13 = m_13 + A.m_13;
	double e_14 = m_14 + A.m_14;

	double e_21 = m_21 + A.m_21;
	double e_22 = m_22 + A.m_22;
	double e_23 = m_23 + A.m_23;
	double e_24 = m_24 + A.m_24;

	double e_31 = m_31 + A.m_31;
	double e_32 = m_32 + A.m_32;
	double e_33 = m_33 + A.m_33;
	double e_34 = m_34 + A.m_34;

	double e_41 = m_41 + A.m_41;
	double e_42 = m_42 + A.m_42;
	double e_43 = m_43 + A.m_43;
	double e_44 = m_44 + A.m_44;

	return TkrFitMatrix(
		e_11,e_12,e_13,e_14,
		e_21,e_22,e_23,e_24,
		e_31,e_32,e_33,e_34,
		e_41,e_42,e_43,e_44);
}

//operator overload for TkrFitMatrix-TkrFitMatrix
const TkrFitMatrix TkrFitMatrix::operator -(const TkrFitMatrix& A) const
{

	double e_11 = m_11 - A.m_11;
	double e_12 = m_12 - A.m_12;
	double e_13 = m_13 - A.m_13;
	double e_14 = m_14 - A.m_14;

	double e_21 = m_21 - A.m_21;
	double e_22 = m_22 - A.m_22;
	double e_23 = m_23 - A.m_23;
	double e_24 = m_24 - A.m_24;

	double e_31 = m_31 - A.m_31;
	double e_32 = m_32 - A.m_32;
	double e_33 = m_33 - A.m_33;
	double e_34 = m_34 - A.m_34;

	double e_41 = m_41 - A.m_41;
	double e_42 = m_42 - A.m_42;
	double e_43 = m_43 - A.m_43;
	double e_44 = m_44 - A.m_44;

	return TkrFitMatrix(
		e_11,e_12,e_13,e_14,
		e_21,e_22,e_23,e_24,
		e_31,e_32,e_33,e_34,
		e_41,e_42,e_43,e_44);
}



//operator overload for TkrFitMatrix*TkrFitMatrix
//!NOTE: Assumed LH multiply
const TkrFitMatrix TkrFitMatrix::operator *(const TkrFitMatrix& A) const
{

	double e_11 = m_11*A.m_11 + m_12*A.m_21 + m_13*A.m_31 + m_14*A.m_41;
	double e_12 = m_11*A.m_12 + m_12*A.m_22 + m_13*A.m_32 + m_14*A.m_42;
	double e_13 = m_11*A.m_13 + m_12*A.m_23 + m_13*A.m_33 + m_14*A.m_43;
	double e_14 = m_11*A.m_14 + m_12*A.m_24 + m_13*A.m_34 + m_14*A.m_44;

	double e_21 = m_21*A.m_11 + m_22*A.m_21 + m_23*A.m_31 + m_24*A.m_41;
	double e_22 = m_21*A.m_12 + m_22*A.m_22 + m_23*A.m_32 + m_24*A.m_42;
	double e_23 = m_21*A.m_13 + m_22*A.m_23 + m_23*A.m_33 + m_24*A.m_43;
	double e_24 = m_21*A.m_14 + m_22*A.m_24 + m_23*A.m_34 + m_24*A.m_44;

	double e_31 = m_31*A.m_11 + m_32*A.m_21 + m_33*A.m_31 + m_34*A.m_41;
	double e_32 = m_31*A.m_12 + m_32*A.m_22 + m_33*A.m_32 + m_34*A.m_42;
	double e_33 = m_31*A.m_13 + m_32*A.m_23 + m_33*A.m_33 + m_34*A.m_43;
	double e_34 = m_31*A.m_14 + m_32*A.m_24 + m_33*A.m_34 + m_34*A.m_44;

	double e_41 = m_41*A.m_11 + m_42*A.m_21 + m_43*A.m_31 + m_44*A.m_41;
	double e_42 = m_41*A.m_12 + m_42*A.m_22 + m_43*A.m_32 + m_44*A.m_42;
	double e_43 = m_41*A.m_13 + m_42*A.m_23 + m_43*A.m_33 + m_44*A.m_43;
	double e_44 = m_41*A.m_14 + m_42*A.m_24 + m_43*A.m_34 + m_44*A.m_44;


	return TkrFitMatrix(
		e_11,e_12,e_13,e_14,
		e_21,e_22,e_23,e_24,
		e_31,e_32,e_33,e_34,
		e_41,e_42,e_43,e_44);
}


const TkrFitPar TkrFitMatrix::operator *(const TkrFitPar &A) const
{
	double xpos = m_11*A.m_xpos + m_12*A.m_xslope + m_13*A.m_ypos + m_14*A.m_yslope;

	double xslope = m_21*A.m_xpos + m_22*A.m_xslope + m_23*A.m_ypos + m_24*A.m_yslope;

	double ypos = m_31*A.m_xpos + m_32*A.m_xslope + m_33*A.m_ypos + m_34*A.m_yslope;

	double yslope = m_41*A.m_xpos + m_42*A.m_xslope + m_43*A.m_ypos + m_44*A.m_yslope;

	return TkrFitPar(xpos,xslope,ypos,yslope);
}


const TkrFitMatrix TkrFitMatrix::operator *(const double &x) const
{
	double e_11 = x*m_11;
	double e_12 = x*m_12;
	double e_13 = x*m_13;
	double e_14 = x*m_14;

	double e_21 = x*m_21;
	double e_22 = x*m_22;
	double e_23 = x*m_23;
	double e_24 = x*m_24;

	double e_31 = x*m_31;
	double e_32 = x*m_32;
	double e_33 = x*m_33;
	double e_34 = x*m_34;

	double e_41 = x*m_41;
	double e_42 = x*m_42;
	double e_43 = x*m_43;
	double e_44 = x*m_44;


	return TkrFitMatrix(
		e_11,e_12,e_13,e_14,
		e_21,e_22,e_23,e_24,
		e_31,e_32,e_33,e_34,
		e_41,e_42,e_43,e_44);
}
int TkrFitMatrix::compare(const TkrFitMatrix& A) const
{

	if( A.m_11 == m_11 || ( abs(A.m_11 - m_11) < 1e-010) )
	{
		if( A.m_12 == m_12 || ( abs(A.m_12 - m_12) < 1e-010) )
		{
			if( A.m_13 == m_13 || ( abs(A.m_13 - m_13) < 1e-010) )
			{
				if( A.m_14 == m_14 || ( abs(A.m_14 - m_14) < 1e-010) )
				{
					if( A.m_21 == m_21 || ( abs(A.m_21 - m_21) < 1e-010) )
					{
						if( A.m_22 == m_22 || ( abs(A.m_22 - m_22) < 1e-010) )
						{
							if( A.m_23 == m_23 || ( abs(A.m_23 - m_23) < 1e-010) )
							{
								if( A.m_24 == m_24 || ( abs(A.m_24 - m_24) < 1e-010) )
								{
									if( A.m_31 == m_31 || ( abs(A.m_31 - m_31) < 1e-010) )
									{
										if( A.m_32 == m_32 || ( abs(A.m_32 - m_32) < 1e-010) )
										{
											if( A.m_33 == m_33 || ( abs(A.m_33 - m_33) < 1e-010) )
											{
												if( A.m_34 == m_34 || ( abs(A.m_34 - m_34) < 1e-010) )
												{
													if( A.m_41 == m_41 || ( abs(A.m_41 - m_41) < 1e-010) )
													{
														if( A.m_42 == m_42 || ( abs(A.m_42 - m_42) < 1e-010) )
														{
															if( A.m_43 == m_43 || ( abs(A.m_43 - m_43) < 1e-010) )
															{
																if( A.m_44 == m_44 || ( abs(A.m_44 - m_44) < 1e-010) )
																{
																	return 0;
																} else {
																	//std::cout << abs(A.m_44 - m_44) << std::endl;
																	return 44;
																}
															} else {
																//std::cout << abs(A.m_43 - m_43) << std::endl;
																return 43;
															}
														} else {
															//std::cout << abs(A.m_42 - m_42) << std::endl;
															return 42;
														}
													} else {
														//std::cout << abs(A.m_41 - m_41) << std::endl;
														return 41;
													}
												} else {
													//std::cout << abs(A.m_34 - m_34) << std::endl;
													return 34;
												}
											} else {
												//std::cout << abs(A.m_33 - m_33) << std::endl;
												return 33;
											}
										} else {
											//std::cout << abs(A.m_32 - m_32) << std::endl;
											return 32;
										}
									} else {
										//std::cout << abs(A.m_31 - m_31) << std::endl;
										return 31;
									}
								} else {
									//std::cout << abs(A.m_24 - m_24) << std::endl;
									return 24;
								}
							} else {
								//std::cout << abs(A.m_23 - m_23) << std::endl;
								return 23;
							}
						} else {
							//std::cout << abs(A.m_22 - m_22) << std::endl;
							return 22;
						}
					} else {
						//std::cout << abs(A.m_21 - m_21) << std::endl;
						return 21;
					}
				} else {
					//std::cout << abs(A.m_14 - m_14) << std::endl;
					return 14;
				}
			} else {
				//std::cout << abs(A.m_13 - m_13) << std::endl;
				return 13;
			}
		} else {
			//std::cout << abs(A.m_12 - m_12) << std::endl;
			return 12;
		}
	} else {
		//std::cout << abs(A.m_11 - m_11) << std::endl;
		return 11;
	}
}