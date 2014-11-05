/*!
 * \file variable_adjoint_turbulent.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.3 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/variable_structure.hpp"

CAdjTurbVariable::CAdjTurbVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  
	dmuT_dUTvar = NULL;
	dRTstar_dUTvar = NULL;
	dFT_dUTvar = NULL;
	EddyViscSens = NULL;
  
}

CAdjTurbVariable::CAdjTurbVariable(double val_psinu_inf, unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
  
  /*--- Array initialization ---*/
  
	dmuT_dUTvar = NULL;
	dRTstar_dUTvar = NULL;
	dFT_dUTvar = NULL;
	EddyViscSens = NULL;
  
	/*--- Initialization of variables ---*/
  
	for (unsigned short iVar = 0; iVar < nVar; iVar++) {
		Solution[iVar] = val_psinu_inf;
		Solution_Old[iVar] = val_psinu_inf;
	}
  
	Residual_Old = new double [nVar];
  
	if (config->GetSpatialOrder() == SECOND_ORDER_LIMITER) Limiter = new double [nVar];
  
}

CAdjTurbVariable::~CAdjTurbVariable(void) {
  
	if (dmuT_dUTvar   != NULL) delete [] dmuT_dUTvar;
	if (EddyViscSens  != NULL) delete [] EddyViscSens;
  
}
