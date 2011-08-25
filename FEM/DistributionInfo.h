/*
 * DistributionInfo.h
 *
 *  Created on: Sep 28, 2010
 *      Author: TF
 */

#ifndef DISTRIBUTIONINFO_H_
#define DISTRIBUTIONINFO_H_
#include <fstream>


// FEM
#include "FEMEnums.h"

using namespace std;

class DistributionInfo
{
   public:
      DistributionInfo();
      virtual ~DistributionInfo();

      /**
       * Sets the value for the distribution type
       * @param dis_type value for primary variable, possible values are documented in enum PrimaryVariable
       */
      void setProcessDistributionType (FiniteElement::DistributionType dis_type);

      /**
       * Get the distribution type of the process.
       * @return the distribution type of the process
       */
      FiniteElement::DistributionType getProcessDistributionType () const;

   private:
      /**
       * the distribution type of the process, see enum DistributionType for valid values
       */
      FiniteElement::DistributionType _dis_type;
};
//--------------------------------------------------------------
/*!
 \class LinearFunctionData
 \brief Define a linear function for IC, BC and ST
  
  WW 24.08.2011
*/
class LinearFunctionData
{
   public: 
     LinearFunctionData(ifstream &ins, const size_t num_var = -1);
     ~LinearFunctionData();

     double getValue(const int dom_i, const double x, const double y, const double z) const;
     double getValue(const double x, const double y, const double z) const;
	 size_t *getSubDomIndex() const {return subdom_index;} 

   private:
      
     size_t ndata;
     size_t *subdom_index; 
     // Coefficents for linear distribution function
     // f = a0+b0*x+c0*y+d0*z
     double *a0, *b0, *c0, *d0;                  

};
#endif                                            /* DISTRIBUTIONINFO_H_ */
