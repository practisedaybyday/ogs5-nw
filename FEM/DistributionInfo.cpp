/*
 * DistributionInfo.cpp
 *
 *  Created on: Sep 28, 2010
 *      Author: fischeth
 */

#include <vector>
#include <sstream>
#include <cstring>

#include "DistributionInfo.h"

DistributionInfo::DistributionInfo() :
_dis_type (FiniteElement::INVALID_DIS_TYPE)
{}

DistributionInfo::~DistributionInfo()
{}

void DistributionInfo::setProcessDistributionType (FiniteElement::DistributionType dis_type)

{
   _dis_type = dis_type;
}


FiniteElement::DistributionType DistributionInfo::getProcessDistributionType () const
{
   return _dis_type;
}

//---------------------------------------------------------------------------------
/*!
 \brief Constrcutor of class LinearFunctionData

  \param ifstream &ins: file
  \param num_var: number of data
  \param sub_domain: if sub_domain

  WW 24.08.2011
*/
LinearFunctionData::LinearFunctionData(ifstream &ins, const int num_var)
{
   string str_buff;
   vector<string> tokens;
   stringstream buff;
   size_t ibuf (0);
 
   bool is_sub_domain = false;
   if(num_var >0) 
   {
     is_sub_domain = true;
     ndata = num_var;
   } 
   else
     ndata = 1;



   char *pch;
   char seps[] = "+\n";
   char seps1[] = "*";
   double f_buff;
   a0 = new double[ndata];
   b0 = new double[ndata];
   c0 = new double[ndata];
   d0 = new double[ndata];

   if(is_sub_domain)  
     subdom_index = new size_t[ndata];
   else
     subdom_index = NULL;


   for (size_t i = 0; i < ndata; i++)
   {
      a0[i] = b0[i] = c0[i] = d0[i] = 0.0;

      if(is_sub_domain) 
	  {
         ins >> ibuf >> str_buff >> ws;
         subdom_index[i] = ibuf;
	  }
	  else
         ins >> str_buff >> ws;
 
      pch = strtok(const_cast<char*> (str_buff.c_str()), seps);
      buff << pch;
      buff >> a0[i];
      buff.clear();
      while (pch != NULL)
      {
          pch = strtok(NULL, seps);
          if (pch == NULL)
             break;
          string token = pch;
             tokens.push_back(token);
      }
      for (size_t k=0; k < tokens.size(); k++)
      {
          pch = strtok( const_cast<char*> (tokens[k].c_str()), seps1);
          buff << pch;
          buff >> f_buff;
          buff.clear();
          pch = strtok(NULL, seps1);
          switch (pch[0])
          {
             case 'x':
               b0[i] = f_buff;
               break;
             case 'y':
               c0[i] = f_buff;
               break;
             case 'z':
               d0[i] = f_buff;
               break;
          }
      }
      tokens.clear();
   }
}
/*!
 \brief Destrcutor of class LinearFunctionData

  WW 24.08.2011
*/
LinearFunctionData::~LinearFunctionData()
{
   delete [] a0;
   delete [] b0;
   delete [] c0;
   delete [] d0;
   if(subdom_index)
     delete [] subdom_index;
   
   subdom_index = NULL;
   a0 = b0 = c0 = d0 = NULL;

}

/*!
 \brief Get Value of class LinearFunctionData

  WW 24.08.2011
*/
double LinearFunctionData::getValue(const int dom_i, const double x, const double y, const double z) const
{
  for (size_t i = 0; i < ndata; i++)
  {
     if(dom_i == subdom_index[i])
	 {
         return a0[i] + b0[i] * x + c0[i] * y + d0[i] * z;
	 } 
  }
  return 0.; 
}
/*!
 \brief Get Value of class LinearFunctionData

  WW 24.08.2011
*/
double LinearFunctionData::getValue(const double x, const double y, const double z) const
{
   return a0[0] + b0[0] * x + c0[0] * y + d0[0] * z;
}

