#include "Event/Recon/TkrRecon/TkrFitMatrix.h"
#include "Event/Recon/TkrRecon/TkrFitPar.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectList.h"
#include <vector>
#include <iostream>
#include <string>



int main() 
{
	// Make the generic objects to use in testing should have det()==0
	Event::TkrFitMatrix m_matrix0 = Event::TkrFitMatrix(10,1,2,3,4,20,5,6,7,8,11,9,3,2,1,16);
	Event::TkrFitPar m_par0 = Event::TkrFitPar(1,2,3,4);
	Event::TkrFitMatrix m_matrix1 = Event::TkrFitMatrix(30,1,1,2,12,1,10,9,8,7,6,5,4,3,2,1);
	Event::TkrFitPar m_par1 = Event::TkrFitPar(4,3,2,1);
	Event::TkrFitMatrix m_identity = Event::TkrFitMatrix(
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1);
	int error_code = 0;

	std::cout << "Testing TkrFitMatrix:" << std::endl;
	// invert matrix 0
	Event::TkrFitMatrix m_matrix0inv = m_matrix0.inverse(error_code);
	std::cout << "Inverted matrix m_matrix0 with error code: " << error_code << std::endl;
	// Invert Matrix 1
	Event::TkrFitMatrix m_matrix1inv = m_matrix1.inverse(error_code);
	std::cout << "Inverted matrix m_matrix1 with error code: " << error_code << std::endl;
	// Test inversion of matrix 0
	Event::TkrFitMatrix m_invtest0 = m_matrix0inv*m_matrix0;
	std::cout << "Inverse test 0 returned: " << m_invtest0.compare(m_identity) << std::endl;
	/*
	std::cout << "Inverse test 0:" << std::endl
		<< m_invtest0(1,1) << " " << m_invtest0(1,2) << " " << m_invtest0(1,3) << " "<< m_invtest0(1,4) << std::endl
		<< m_invtest0(2,1) << " " << m_invtest0(2,2) << " " << m_invtest0(2,3) << " "<< m_invtest0(2,4) << std::endl
		<< m_invtest0(3,1) << " " << m_invtest0(3,2) << " " << m_invtest0(3,3) << " "<< m_invtest0(3,4) << std::endl
		<< m_invtest0(4,1) << " " << m_invtest0(4,2) << " " << m_invtest0(4,3) << " "<< m_invtest0(4,4) << std::endl;
	*/

	// Test inversion of matrix 1
	Event::TkrFitMatrix m_invtest1 = m_matrix1inv*m_matrix1;
	/*
	std::cout << "Inverse test 1:" << std::endl
		<< m_invtest1(1,1) << " " << m_invtest1(1,2) << " " << m_invtest1(1,3) << " "<< m_invtest1(1,4) << std::endl
		<< m_invtest1(2,1) << " " << m_invtest1(2,2) << " " << m_invtest1(2,3) << " "<< m_invtest1(2,4) << std::endl
		<< m_invtest1(3,1) << " " << m_invtest1(3,2) << " " << m_invtest1(3,3) << " "<< m_invtest1(3,4) << std::endl
		<< m_invtest1(4,1) << " " << m_invtest1(4,2) << " " << m_invtest1(4,3) << " "<< m_invtest1(4,4) << std::endl;
  */
	std::cout << "Inverse test 1 returned: " << m_invtest1.compare(m_identity) << std::endl;

	// Test TkrFitMatrix*TkrFitPar
	std::cout << "TkrFitMatrix*TkrFitPar test 0 returned:" << ((m_identity*m_par0)==m_par0) << std::endl;
  std::cout << "TkrFitMatrix*TkrFitPar test 1 returned:" << ((m_identity*m_par1)==m_par1) << std::endl;

	return 0;
}
