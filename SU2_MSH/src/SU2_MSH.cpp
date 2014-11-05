/*!
 * \file SU2_MSH.cpp
 * \brief Main file of Mesh Adaptation Code (SU2_MSH).
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

#include "../include/SU2_MSH.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	
	/*--- Variable definitions ---*/
  
	char file_name[MAX_STRING_SIZE];
  unsigned short nZone = 1;
  
#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
#endif
	
	/*--- Definition of the config problem ---*/
  
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_MSH, ZONE_0, nZone, 0, VERB_HIGH);
	else { strcpy (file_name, "default.cfg"); config = new CConfig(file_name, SU2_MSH, ZONE_0, nZone, 0, VERB_HIGH); }
	
	/*--- Definition of the Class for the geometry ---*/
  
	CGeometry *geometry; geometry = new CGeometry;
	geometry = new CPhysicalGeometry(config, ZONE_0, nZone);
  
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute elements surrounding points, points surrounding points, and elements surronding elements ---*/
  
	cout << "Setting local point and element connectivity." <<endl;
	geometry->SetPoint_Connectivity(); geometry->SetElement_Connectivity();
	
	/*--- Check the orientation before computing geometrical quantities ---*/
  
	cout << "Check numerical grid orientation." <<endl;
	geometry->SetBoundVolume(); geometry->Check_IntElem_Orientation(config); geometry->Check_BoundElem_Orientation(config);
	
	/*--- Create the edge structure ---*/
  
	cout << "Identify faces, edges and vertices." <<endl;
	geometry->SetFaces(); geometry->SetEdges(); geometry->SetVertex(config); geometry->SetCG();
	
	/*--- Create the control volume structures ---*/
  
	cout << "Set control volume structure." << endl;
	geometry->SetControlVolume(config, ALLOCATE); geometry->SetBoundControlVolume(config, ALLOCATE);

	
	if ((config->GetKind_Adaptation() != NONE) && (config->GetKind_Adaptation() != PERIODIC)) {
		
		cout << endl <<"--------------------- Start numerical grid adaptation -------------------" << endl;
		
		/*-- Definition of the Class for grid adaptation ---*/
    
		CGridAdaptation *grid_adaptation;
		grid_adaptation = new CGridAdaptation(geometry, config);
		
		/*--- Read the flow solution and/or the adjoint solution
		 and choose the elements to adapt ---*/
    
		if ((config->GetKind_Adaptation() != FULL)
				&& (config->GetKind_Adaptation() != WAKE) && (config->GetKind_Adaptation() != TWOPHASE)
				&& (config->GetKind_Adaptation() != SMOOTHING) && (config->GetKind_Adaptation() != SUPERSONIC_SHOCK))
			grid_adaptation->GetFlowSolution(geometry, config);
		
		switch (config->GetKind_Adaptation()) {
			case NONE:
				break;
			case SMOOTHING:
				config->SetSmoothNumGrid(true);
				grid_adaptation->SetNo_Refinement(geometry, 1);
				break;
			case FULL:
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case WAKE:
				grid_adaptation->SetWake_Refinement(geometry, 1);
				break;
			case TWOPHASE:
				grid_adaptation->SetTwoPhase_Refinement(geometry, 1);
				break;
			case SUPERSONIC_SHOCK:
				grid_adaptation->SetSupShock_Refinement(geometry, config);
				break;
			case FULL_FLOW:
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case FULL_ADJOINT:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case FULL_LINEAR:
				grid_adaptation->GetLinSolution(geometry, config);
				grid_adaptation->SetComplete_Refinement(geometry, 1);
				break;
			case GRAD_FLOW:
				grid_adaptation->SetIndicator_Flow(geometry, config, 1);
				break;
			case GRAD_ADJOINT:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->SetIndicator_Adj(geometry, config, 1);
				break;
			case GRAD_FLOW_ADJ:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->SetIndicator_FlowAdj(geometry, config);
				break;
			case COMPUTABLE:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->GetFlowResidual(geometry, config);
				grid_adaptation->SetIndicator_Computable(geometry, config);
				break;
			case REMAINING:
				cout << "Adaptation method not implemented."<< endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
				exit(1);
				break;
			case ROBUST:
				grid_adaptation->GetFlowResidual(geometry, config);
				grid_adaptation->GetAdjResidual(geometry, config);
				grid_adaptation->SetIndicator_Robust(geometry, config);
				break;
			case COMPUTABLE_ROBUST:
				grid_adaptation->GetAdjSolution(geometry, config);
				grid_adaptation->GetLinResidual(geometry, config);
				grid_adaptation->SetIndicator_Computable_Robust(geometry, config);
				break;
			default :
				cout << "The adaptation is not defined" << endl;
		}
		
		/*--- Perform an homothetic adaptation of the grid ---*/
    
		CPhysicalGeometry *geo_adapt; geo_adapt = new CPhysicalGeometry;
		
		cout << "Homothetic grid adaptation" << endl;
		if (geometry->GetnDim() == 2) grid_adaptation->SetHomothetic_Adaptation2D(geometry, geo_adapt, config);
		if (geometry->GetnDim() == 3) grid_adaptation->SetHomothetic_Adaptation3D(geometry, geo_adapt, config);
    
		/*--- Smooth the numerical grid coordinates ---*/
    
		if (config->GetSmoothNumGrid()) {
			cout << "Preprocessing for doing the implicit smoothing." << endl;
			geo_adapt->SetPoint_Connectivity(); geo_adapt->SetElement_Connectivity();
			geo_adapt->SetBoundVolume(); geo_adapt->Check_IntElem_Orientation(config); geo_adapt->Check_BoundElem_Orientation(config);
			geo_adapt->SetEdges(); geo_adapt->SetVertex(config);
			cout << "Implicit smoothing of the numerical grid coordinates." << endl;
			geo_adapt->SetCoord_Smoothing(5, 1.5, config);
		}
		
		/*--- Original and adapted grid ---*/
    strcpy (file_name, "original_grid.plt");
    geometry->SetTecPlot(file_name, true);
    strcpy (file_name, "original_surface.plt");
    geometry->SetBoundTecPlot(file_name, true, config);
    
		/*--- Write the adapted grid sensor ---*/
    
    strcpy (file_name, "adapted_grid.plt");
    geo_adapt->SetTecPlot(file_name, true);
    strcpy (file_name, "adapted_surface.plt");
    geo_adapt->SetBoundTecPlot(file_name, true, config);
		
		/*--- Write the new adapted grid, including the modified boundaries surfaces ---*/
    
		geo_adapt->SetMeshFile(config, config->GetMesh_Out_FileName());
    
    
		/*--- Write the restart file ---*/
    
		if ((config->GetKind_Adaptation() != SMOOTHING) && (config->GetKind_Adaptation() != FULL) &&
				(config->GetKind_Adaptation() != WAKE) && (config->GetKind_Adaptation() != TWOPHASE) &&
				(config->GetKind_Adaptation() != SUPERSONIC_SHOCK))
			grid_adaptation->SetRestart_FlowSolution(config, geo_adapt, config->GetRestart_FlowFileName());
		
		if ((config->GetKind_Adaptation() == GRAD_FLOW_ADJ) || (config->GetKind_Adaptation() == GRAD_ADJOINT)
				|| (config->GetKind_Adaptation() == FULL_ADJOINT) || (config->GetKind_Adaptation() == ROBUST)
				|| (config->GetKind_Adaptation() == COMPUTABLE) || (config->GetKind_Adaptation() == COMPUTABLE_ROBUST) ||
				(config->GetKind_Adaptation() == REMAINING))
			grid_adaptation->SetRestart_AdjSolution(config, geo_adapt, config->GetRestart_AdjFileName());
		
		if ((config->GetKind_Adaptation() == FULL_LINEAR) || (config->GetKind_Adaptation() == COMPUTABLE_ROBUST)) {
			grid_adaptation->SetRestart_LinSolution(config, geo_adapt, config->GetRestart_LinFileName());
		}
	}
	else {
    
    if (config->GetKind_Adaptation() == PERIODIC) {
      
      cout << endl <<"-------------------- Setting the periodic boundaries --------------------" << endl;
      
      /*--- Set periodic boundary conditions ---*/
      
      geometry->SetPeriodicBoundary(config);
      
      /*--- Original grid for debugging purposes ---*/
      
      strcpy (file_name, "periodic_original.plt"); geometry->SetTecPlot(file_name);
      
      /*--- Create a new grid with the right periodic boundary ---*/
      
      CGeometry *periodic; periodic = new CPeriodicGeometry(geometry, config);
      periodic->SetPeriodicBoundary(geometry, config);
      periodic->SetMeshFile(geometry, config, config->GetMesh_Out_FileName());
      
      /*--- Output of the grid for debuging purposes ---*/
      
      strcpy (file_name, "periodic_halo.plt"); periodic->SetTecPlot(file_name);
      
    }
    
    if (config->GetKind_Adaptation() == NONE) {
      strcpy (file_name, "original_grid.plt");
      geometry->SetTecPlot(file_name);
      geometry->SetMeshFile (config, config->GetMesh_Out_FileName());
    }
    
	}
  
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	
	/*--- End solver ---*/
  
	cout << endl <<"------------------------- Exit Success (SU2_MSH) ------------------------" << endl << endl;
  
	return EXIT_SUCCESS;
}

