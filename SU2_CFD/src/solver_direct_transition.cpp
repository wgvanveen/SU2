/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CTransLMSolver::CTransLMSolver(void) : CTurbSolver() {}

CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolver() {
	unsigned short iVar, iDim, nLineLets;
	unsigned long iPoint, index;
	double Viscosity_Inf, tu_Inf, dull_val, rey, mach;
  ifstream restart_file;
  char *cstr;
  string text_line;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
	bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  cout << "Entered constructor for CTransLMSolver -AA\n";
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constans in the solver structure ---*/
  nDim = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
<<<<<<< HEAD

  node = new CVariable*[geometry->GetnPoint()];

  /*--- Dimension of the problem --> 2 Transport equations (intermittency,Reth) ---*/
  nVar = 2;

  if (iMesh == MESH_0) {

    /*--- Define some auxillary vectors related to the residual ---*/
    Residual     = new double[nVar]; Residual_RMS = new double[nVar];
    Residual_i   = new double[nVar]; Residual_j   = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max    = new unsigned long[nVar];
=======
  
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem --> 2 Transport equations (intermittency,Reth) ---*/
	nVar = 2;
	
	if (iMesh == MESH_0) {
		
		/*--- Define some auxillary vectors related to the residual ---*/
		Residual     = new double[nVar]; Residual_RMS = new double[nVar];
		Residual_i   = new double[nVar]; Residual_j   = new double[nVar];
    Residual_Max = new double[nVar];

    /*--- Define some structures for locating max residuals ---*/
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
    Point_Max_Coord = new double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }
>>>>>>> master

    /*--- Define some auxiliar vector related with the solution ---*/
    Solution   = new double[nVar];
    Solution_i = new double[nVar]; Solution_j = new double[nVar];

    /*--- Define some auxiliar vector related with the geometry ---*/
    Vector_i = new double[nDim]; Vector_j = new double[nDim];

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {

			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar];
			Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar];
				Jacobian_j[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
      
<<<<<<< HEAD

      if (config->GetKind_Linear_Solver_Prec() == LINELET) {
=======
      if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
          (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
>>>>>>> master
        nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
        if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
      }

    }

    /*--- Computation of gradients by least squares ---*/
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new double [nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
      cvector = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        cvector[iVar] = new double [nDim];
    }
  }

  /*--- Read farfield conditions from config ---*/
  Density_Inf       = config->GetDensity_FreeStreamND();
  Viscosity_Inf     = config->GetViscosity_FreeStreamND();
  Intermittency_Inf = config->GetIntermittency_FreeStream();
	tu_Inf            = config->GetTurbulenceIntensity_FreeStream()*100;

  /*-- Initialize REth from correlation --*/
  if (tu_Inf <= 1.3) {
    REth_Inf = (1173.51-589.428*tu_Inf+0.2196/(tu_Inf*tu_Inf));
  } else {
    REth_Inf = 331.5*pow(tu_Inf-0.5658,-0.671);
  }
  rey = config->GetReynolds();
  mach = config->GetMach();

  //  REth_Inf *= mach/rey;
  cout << "tu_Inf = " << tu_Inf << endl;
  cout << "REth_Inf = " << REth_Inf << ", rey: "<< rey << " -AA" << endl;

  /*--- Define some auxiliar vector related with the flow solution ---*/
  FlowPrimVar_i = new double [nDim+7]; FlowPrimVar_j = new double [nDim+7];

  /*--- Restart the solution from file information ---*/
  if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CTransLMVariable(Density_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);
  }

  else {
		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();

        /*--- Modify file name for an unsteady restart ---*/
        if (dual_time) {
            int Unst_RestartIter;
            if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
                Unst_RestartIter = int(config->GetUnst_RestartIter())-1;
            else
                Unst_RestartIter = int(config->GetUnst_RestartIter())-2;
            filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
        }

        /*--- Open the restart file, throw an error if this fails. ---*/
		restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no turbulent restart file!!" << endl;
      exit(1);
    }

<<<<<<< HEAD
		/*--- In case this is a parallel simulation, we need to perform the
         Global2Local index transformation first. ---*/
		long *Global2Local;
		Global2Local = new long[geometry->GetGlobal_nPointDomain()];
		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
			Global2Local[iPoint] = -1;
		}
		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
=======
	}
	else {
    cout << "No LM restart yet!!" << endl; // TODO, Aniket
    int j;
    cin >> j;
		string mesh_filename = config->GetSolution_FlowFileName();
		cstr = new char [mesh_filename.size()+1];
		strcpy (cstr, mesh_filename.c_str());
		restart_file.open(cstr, ios::in);
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			exit(EXIT_FAILURE);
>>>>>>> master
		}

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;

        /*--- The first line is the header ---*/
        getline (restart_file, text_line);

		while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);

			/*--- Retrieve local index. If this node from the restart file lives
             on a different processor, the value of iPoint_Local will be -1.
             Otherwise, the local index for this node on the current processor
             will be returned and used to instantiate the vars. ---*/
			iPoint_Local = Global2Local[iPoint_Global];

			if (iPoint_Local >= 0) {

				double Density;
				if (compressible) {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Density >> dull_val>> dull_val>> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Density >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];

				}
				if (incompressible) {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				}

                if (freesurface) {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >>  dull_val >> dull_val >> Solution[0] >> Solution[1];
				}

				/*--- Instantiate the solution at this node, note that the eddy viscosity should be recomputed ---*/
    			node[iPoint_Local] = new CTransLMVariable(Density, Solution[0]/Density, Solution[1]/Density, nDim, nVar, config);

    			/*--- Compute gamma_eff, which will be needed in first iteration of turbulence model ---*/
    			node[iPoint_Local]->SetGammaEff(1.0); // TODO: Need to recompute gamma_eff for restart
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
         at any halo/periodic nodes. The initial solution can be arbitrary,
         because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CTransLMVariable(Density_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);
		}

		/*--- Close the restart file ---*/
    restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
  }

    /*--- MPI solution ---*/
    Set_MPI_Solution(geometry, config);

}

CTransLMSolver::~CTransLMSolver(void){
  unsigned short iVar, iDim;

  delete [] Residual; delete [] Residual_Max;
  delete [] Residual_i; delete [] Residual_j;
  delete [] Solution;
  delete [] Solution_i; delete [] Solution_j;
  delete [] Vector_i; delete [] Vector_j;

  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i; delete [] Jacobian_j;

  for (iDim = 0; iDim < this->nDim; iDim++)
    delete [] Smatrix[iDim];
  delete [] Smatrix;

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] cvector[iVar];
  delete [] cvector;
}

void CTransLMSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
	double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_gamma_eff = NULL, *Buffer_Send_gamma_eff = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
            nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_U = new double [nBufferR_Vector];
            Buffer_Send_U = new double[nBufferS_Vector];

            Buffer_Receive_gamma_eff = new double [nBufferR_Scalar];
            Buffer_Send_gamma_eff = new double[nBufferS_Scalar];

            /*--- Copy the solution that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                Buffer_Send_gamma_eff[iVertex] = node[iPoint]->GetGammaEff();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
            }

#ifndef NO_MPI

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_gamma_eff, nBufferS_Scalar, MPI::DOUBLE, send_to, 1,
                                     Buffer_Receive_gamma_eff, nBufferR_Scalar, MPI::DOUBLE, receive_from, 1);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                Buffer_Receive_gamma_eff[iVertex] = node[iPoint]->GetGammaEff();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_U;
            delete [] Buffer_Send_gamma_eff;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

                /*--- Copy conservative variables. ---*/
                node[iPoint]->SetGammaEff(Buffer_Receive_gamma_eff[iVertex]);
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_gamma_eff;
            delete [] Buffer_Receive_U;

        }

	}

}

void CTransLMSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_U = new double [nBufferR_Vector];
            Buffer_Send_U = new double[nBufferS_Vector];

            /*--- Copy the solution old that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
            }

#ifndef NO_MPI

            //      /*--- Send/Receive using non-blocking communications ---*/
            //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
            //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
            //      send_request.Wait(status);
            //      recv_request.Wait(status);

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_U;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

                /*--- Copy transformed conserved variables back into buffer. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_U;

        }

	}
}

void CTransLMSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
    *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_Limit = new double [nBufferR_Vector];
            Buffer_Send_Limit = new double[nBufferS_Vector];

            /*--- Copy the solution old that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
            }

#ifndef NO_MPI

            //      /*--- Send/Receive using non-blocking communications ---*/
            //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
            //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
            //      send_request.Wait(status);
            //      recv_request.Wait(status);

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_Limit;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

                /*--- Retrieve the supplied periodic information. ---*/
                angles = config->GetPeriodicRotation(iPeriodic_Index);

                /*--- Store angles separately for clarity. ---*/
                theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
                cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
                sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

                /*--- Compute the rotation matrix. Note that the implicit
                 ordering is rotation about the x-axis, y-axis,
                 then z-axis. Note that this is the transpose of the matrix
                 used during the preprocessing stage. ---*/
                rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

                /*--- Copy conserved variables before performing transformation. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];

                /*--- Rotate the momentum components. ---*/
                if (nDim == 2) {
                    Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
                    Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
                }
                else {
                    Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
                    rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
                    Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
                    rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
                    Solution[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
                    rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
                }

                /*--- Copy transformed conserved variables back into buffer. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetLimiter(iVar, Solution[iVar]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_Limit;

        }

	}
}

void CTransLMSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
    *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

    double **Gradient = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
        Gradient[iVar] = new double[nDim];

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_Gradient = new double [nBufferR_Vector];
            Buffer_Send_Gradient = new double[nBufferS_Vector];

            /*--- Copy the solution old that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
            }

#ifndef NO_MPI

            //      /*--- Send/Receive using non-blocking communications ---*/
            //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
            //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
            //      send_request.Wait(status);
            //      recv_request.Wait(status);

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_Gradient;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

                /*--- Retrieve the supplied periodic information. ---*/
                angles = config->GetPeriodicRotation(iPeriodic_Index);

                /*--- Store angles separately for clarity. ---*/
                theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
                cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
                sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

                /*--- Compute the rotation matrix. Note that the implicit
                 ordering is rotation about the x-axis, y-axis,
                 then z-axis. Note that this is the transpose of the matrix
                 used during the preprocessing stage. ---*/
                rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

                /*--- Copy conserved variables before performing transformation. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];

                /*--- Need to rotate the gradients for all conserved variables. ---*/
                for (iVar = 0; iVar < nVar; iVar++) {
                    if (nDim == 2) {
                        Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
                        Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
                    }
                    else {
                        Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
                        Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
                        Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
                    }
                }

                /*--- Store the received information ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_Gradient;

        }

	}

    for (iVar = 0; iVar < nVar; iVar++)
        delete [] Gradient[iVar];
    delete [] Gradient;

}

void CTransLMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		LinSysRes.SetBlock_Zero(iPoint);
  Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CTransLMSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

}

void CTransLMSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
	unsigned short iVar;
	unsigned long iPoint, total_index, IterLinSol;
	double Delta, Delta_flow, Vol;
  
  
  /*--- Set maximum residual to zero ---*/
  
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
	/*--- Build implicit system ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta_flow = Vol / (solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Delta = Delta_flow;
    Jacobian.AddVal2Diag(iPoint,Delta);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      
      /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
      
      LinSysRes[total_index] = -LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]*Vol);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
	
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
	/*--- Update solution (system written in terms of increments) ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
<<<<<<< HEAD
    if (iPoint==99) {
      double *sol = node[iPoint]->GetSolution();

    }
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->AddClippedSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar], 0.0, 1E10);
    }
    if (iPoint==99) {
      double *sol = node[iPoint]->GetSolution();

    }
  }
  
    /*--- MPI solution ---*/
    Set_MPI_Solution(geometry, config);
    
    /*--- Compute the root mean square residual ---*/
    SetResidual_RMS(geometry, config);

    iPoint=99;
    double *sol = node[iPoint]->GetSolution();

=======
    for (iVar = 0; iVar < nVar; iVar++)
      node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
>>>>>>> master
}

void CTransLMSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {

  double *Trans_i, *Trans_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool second_order  = ((config->GetSpatialOrder() == SECOND_ORDER) || (config->GetSpatialOrder() == SECOND_ORDER_LIMITER));
  bool high_order_diss = true;
  bool grid_movement = config->GetGrid_Movement();
  bool limiter = (config->GetKind_SlopeLimit() != NONE);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Conservative variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Transulent variables w/o reconstruction ---*/
    
    Trans_i = node[iPoint]->GetSolution();
    Trans_j = node[jPoint]->GetSolution();
    numerics->SetTransVar(Trans_i,Trans_j);
    
    /*--- Grid Movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    if (second_order) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Mean flow primitive variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = solver_container[FLOW_SOL]->node[iPoint]->GetLimiter_Primitive();
        Limiter_j = solver_container[FLOW_SOL]->node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          FlowPrimVar_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);
      
      /*--- Turbulent variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        Solution_i[iVar] = Trans_i[iVar] + Project_Grad_i;
        Solution_j[iVar] = Trans_j[iVar] + Project_Grad_j;
      }
      
      numerics->SetTransVar(Solution_i,Solution_j);
      
    }
    
    /*--- Add and subtract residual ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);
    
    /*--- Implicit part ---*/
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
  }

}


void CTransLMSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																				CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive());
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(),
                                 solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive());
    
    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), solver_container[TURB_SOL]->node[jPoint]->GetSolution());
    numerics->SetTurbVarGradient(solver_container[TURB_SOL]->node[iPoint]->GetGradient(), solver_container[TURB_SOL]->node[jPoint]->GetGradient());

    /*--- Transulent variables w/o reconstruction, and its gradients ---*/
    numerics->SetTransVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTransVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Menter's first blending function (only SST)---*/
    if (config->GetKind_Trans_Model() == SST)
      numerics->SetF1blending(node[iPoint]->GetF1blending(),node[jPoint]->GetF1blending());
    
    /*--- Compute residual, and Jacobians ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
	
}

void CTransLMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
  unsigned long iPoint, iVertex;
  unsigned short iMarker;
  double gamma_eff = 1.;
  bool boundary;
  ofstream sagt_debug;

<<<<<<< HEAD
    //cout << "Setting Trans residual -AA " << endl;
    //cout << "\nBeginAA" << endl;
=======
  //cout << "Setting Trans residual -AA " << endl;
  cout << "\nBeginAA" << endl;
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    cout << "\niPoint: " << iPoint << endl;
		
		/*--- Conservative variables w/o reconstruction ---*/
		numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);
		
		/*--- Gradient of the primitive and conservative variables ---*/
		numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
		
		/*--- Laminar and eddy viscosity ---*/
		numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
    numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);
		
		/*--- Turbulent variables w/o reconstruction, and its gradient ---*/
		numerics->SetTransVar(node[iPoint]->GetSolution(), NULL);
		// numerics->SetTransVarGradient(node[iPoint]->GetGradient(), NULL);  // Is this needed??
		
		/*--- Set volume ---*/
		numerics->SetVolume(geometry->node[iPoint]->GetVolume());
		
		/*--- Set distance to the surface ---*/
		numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
		
		/*--- Compute the source term ---*/
		numerics->ComputeResidual_TransLM(Residual, Jacobian_i, NULL, config, gamma_sep);
		
    /*-- Store gamma_sep in variable class --*/
    node[iPoint]->SetGammaSep(gamma_sep);

		/*--- Subtract residual and the Jacobian ---*/
		LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
>>>>>>> master

   //static int start_counter=0;
   //if (start_counter < 100) {
   //  cout << "Skipping source term" << endl;
   //  start_counter += 1;
   //  return;
   //} 
  

  // DEBUG
  sagt_debug.open("sagt_debug.plt");
  sagt_debug << "TITLE = \"SAGT (Langtry+Menter) Transition model debug file \" " << endl;
  sagt_debug << "VARIABLES = \"iPoint\" \"itmc\" \"Re_th_bar\" \"re_theta_t\" \"flen\" \"re_theta_c\" "; 
  sagt_debug << "\"val_resid[0]\" \"val_resid[1]\" \"strain\" \"vorticity\" \"tu\" \"lambda\" \"f_lambda\" ";
  sagt_debug << "\"time_scale\" \"f_theta\" \"du_ds\" \"dist_i\" \"Volume\" \"delta\" \"var1\" ";
  sagt_debug << "\"dv_dx\" \"du_dy\" \"re_v\" \"f_onset\" \"f_onset1\" \"f_onset2\" \"f_onset3\" \"gamma_eff\" " << endl;
  sagt_debug << "ZONE DATAPACKING=POINT" << endl;

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

	  /*--- Conservative variables w/o reconstruction ---*/
	  numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

	  /*--- Gradient of the primitive and conservative variables ---*/
	  numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

	  /*--- Laminar and eddy viscosity ---*/
	  numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
	  numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);

	  /*--- Turbulent variables ---*/
	  numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), NULL);

	  /*--- Transition Variables ---*/
	  numerics->SetTransVar(node[iPoint]->GetSolution(), NULL);

	  /*--- Set volume ---*/
	  numerics->SetVolume(geometry->node[iPoint]->GetVolume());

	  /*--- Set distance to the surface ---*/
	  numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

	  /*--- Set distance to the surface ---*/
	  boundary = geometry->node[iPoint]->GetBoundary();

	  /*--- Compute the source term ---*/
    sagt_debug << iPoint << " "; 
	  numerics->ComputeResidual_TransLM(Residual, Jacobian_i, gamma_eff, config, boundary, sagt_debug);

    /*--- Store gamma_eff (eq. 19) in CVariable so it can be accessed by the turbulent solver --- */
    node[iPoint]->SetGammaEff(gamma_eff);
    
    // DEBUG: AA, check the derivatives against finite differencing
    // if (iPoint==9801 && config->GetExtIter()==5) {

    //   ofstream jac_check; 
    //   jac_check.open("jac_check.txt");

    //   double TransVar_i[2];
    //   TransVar_i[0] = 0.5; TransVar_i[1] = 200;
	  //   numerics->SetTransVar(TransVar_i, NULL);
	  //   numerics->ComputeResidual_TransLM(Residual, Jacobian_i, gamma_eff, config, boundary, sagt_debug);

    //   jac_check << endl << endl;
    //   jac_check << "0.0 " << "\t";
    //   jac_check << Jacobian_i[0][0] << "\t " << Jacobian_i[0][1] << "\t";
    //   jac_check << Jacobian_i[1][0] << "\t " << Jacobian_i[1][1] << endl;

    //   for (int j=-4; j<10; j++) {
    //     double delta_u = pow(10,-1.0*j);
    //     double TransVar_i_plus[2], Residual_plus[2], Jacobian_fd[2][2];
    //     //double *TransVar_i;
    //     //TransVar_i = node[iPoint]->GetSolution();

    //     // Perturb TransVar_i[0]
    //     TransVar_i_plus[0] = TransVar_i[0] + delta_u;
    //     TransVar_i_plus[1] = TransVar_i[1];
	  //     numerics->SetTransVar(TransVar_i_plus, NULL);
    //     numerics->ComputeResidual_TransLM(Residual_plus, Jacobian_i, gamma_eff, config, boundary, sagt_debug);
    //     Jacobian_fd[0][0] = (Residual_plus[0]-Residual[0])/delta_u;
    //     Jacobian_fd[1][0] = (Residual_plus[1]-Residual[1])/delta_u;

    //     // Perturb TransVar_i[1]
    //     TransVar_i_plus[0] = TransVar_i[0];
    //     TransVar_i_plus[1] = TransVar_i[1]+delta_u;
	  //     numerics->SetTransVar(TransVar_i_plus, NULL);
    //     numerics->ComputeResidual_TransLM(Residual_plus, Jacobian_i, gamma_eff, config, boundary, sagt_debug);
    //     Jacobian_fd[0][1] = (Residual_plus[0]-Residual[0])/delta_u;
    //     Jacobian_fd[1][1] = (Residual_plus[1]-Residual[1])/delta_u;

    //     jac_check << delta_u << " " ;
    //     jac_check << Jacobian_fd[0][0] << "\t " << Jacobian_fd[0][1] << "\t";
    //     jac_check << Jacobian_fd[1][0] << "\t " << Jacobian_fd[1][1] << endl;
    //   }
    //   exit(1); 
    //} 

	  /*--- Subtract residual and the jacobian ---*/
	  LinSysRes.SubtractBlock(iPoint, Residual);
	  Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }
  
   sagt_debug.close();
  
}

void CTransLMSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																			 CConfig *config, unsigned short iMesh) {
}

void CTransLMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  unsigned short iVar;
  int total_index;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /* --- Impose boundary values (Dirichlet) ---*/
      Solution[0] = 0.0;
      Solution[1] = 0.0;
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- includes 1 in the diagonal ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }

}

void CTransLMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set transition variable at the wall, and at infinity ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0] = Intermittency_Inf*Density_Inf;
      //Solution_j[1] = REth_Inf*Density_Inf;
      Solution_j[1] = Solution_i[1]; 

      conv_numerics->SetTransVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute residuals and jacobians ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Add residuals and jacobians ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

    }
  }
  
  delete [] Normal;

  }

  void CTransLMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double *V_inlet, *V_domain, *Normal;
  
  Normal = new double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
    
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the transition variable states. Use free-stream
       values for the turbulent state at the inflow. ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0]= Intermittency_Inf*Density_Inf;
      Solution_j[1]= REth_Inf*Density_Inf;
      
      conv_numerics->SetTransVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        
      }
    }
    
    /*--- Free locally allocated memory ---*/
    delete[] Normal;
    //cin.get();
  }

  void CTransLMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker) {

  /*--- Local variables and initialization. ---*/
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Allocate the value at the outlet ---*/
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the transition variables. Here we use a Neumann BC such
       that the transition variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TransVar_internal,
       Solution_j --> TransVar_outlet ---*/
      Solution_i = node[iPoint]->GetSolution();
      Solution_j = node[iPoint]->GetSolution();
      conv_numerics->SetTransVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_outlet);
      visc_numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(),
                                 solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive());
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), solver_container[TURB_SOL]->node[iPoint]->GetSolution());
      visc_numerics->SetTurbVarGradient(solver_container[TURB_SOL]->node[iPoint]->GetGradient(), solver_container[TURB_SOL]->node[iPoint]->GetGradient());
      visc_numerics->SetTransVar(Solution_i, Solution_j);
      visc_numerics->SetTransVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
    
  }

  void CTransLMSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker) {

  }
