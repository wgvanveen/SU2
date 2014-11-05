/*!
 * \file SU2_SOL.cpp
 * \brief Main file for the solution export/conversion code (SU2_SOL).
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

#include "../include/SU2_SOL.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	/*--- Variable definitions ---*/
	unsigned short iZone, nZone;
	ofstream ConvHist_file;
	char file_name[MAX_STRING_SIZE];
	int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
	/*--- MPI initialization, and buffer setting ---*/
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
  
	/*--- Pointer to different structures that will be used throughout the entire code ---*/
	COutput *output = NULL;
	CGeometry **geometry = NULL;
	CSolver **solver = NULL;
	CConfig **config = NULL;
	
	/*--- Definition of the containers per zones ---*/
	solver = new CSolver*[MAX_ZONES];
	config = new CConfig*[MAX_ZONES];
	geometry = new CGeometry *[MAX_ZONES];
	
	/*--- Only one zone is allowed ---*/
	nZone = 1;
	
	for (iZone = 0; iZone < nZone; iZone++) {
		
		/*--- Definition of the configuration class per zones ---*/
		if (argc == 2) config[iZone] = new CConfig(argv[1], SU2_SOL, iZone, nZone, 0, VERB_HIGH);
		else { strcpy (file_name, "default.cfg"); config[iZone] = new CConfig(file_name, SU2_SOL,
                                                                          iZone, nZone, 0, VERB_HIGH); }
		
#ifdef HAVE_MPI
		/*--- Change the name of the input-output files for a parallel computation ---*/
		config[iZone]->SetFileNameDomain(rank+1);
#endif
    
		/*--- Definition of the geometry class and open the mesh file ---*/
		geometry[iZone] = new CPhysicalGeometry(config[iZone], iZone+1, nZone);
    
    /*--- Create the vertex structure (required for MPI) ---*/
    if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
    geometry[iZone]->SetVertex(config[iZone]);
    
  }
  
#ifdef HAVE_MPI
  /*--- Synchronization point after the geometrical definition subroutine ---*/
	MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solution Postprocessing -----------------------" << endl;
  
#ifdef HAVE_MPI
  /*--- Synchronization point after the solution subroutine ---*/
	MPI_Barrier(MPI_COMM_WORLD);
#endif
  
	/*--- Definition of the output class (one for all the zones) ---*/
	output = new COutput();
  
  /*---  Check whether this is an unsteady simulation, and call the
   solution merging routines accordingly.---*/
  
  if (config[ZONE_0]->GetWrt_Unsteady()) {
    
    /*--- Unsteady simulation: merge all unsteady time steps. First,
     find the frequency and total number of files to write. ---*/
    
    double Physical_dt, Physical_t;
    unsigned long iExtIter = 0;
    bool StopCalc = false;
    bool SolutionInstantiated = false;
    
    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config[ZONE_0]->GetWrt_Unsteady() && config[ZONE_0]->GetRestart())
      iExtIter = config[ZONE_0]->GetUnst_RestartIter();
    
    while (iExtIter < config[ZONE_0]->GetnExtIter()) {
      
      /*--- Check several conditions in order to merge the correct time step files. ---*/
      Physical_dt = config[ZONE_0]->GetDelta_UnstTime();
      Physical_t  = (iExtIter+1)*Physical_dt;
      if (Physical_t >=  config[ZONE_0]->GetTotal_UnstTime())
        StopCalc = true;
        
      if ((iExtIter+1 == config[ZONE_0]->GetnExtIter()) ||
          ((iExtIter % config[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
           !((config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
             (config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
          (StopCalc) ||
          (((config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
            (config[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
           ((iExtIter == 0) || (iExtIter % config[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

        /*--- Set the current iteration number in the config class. ---*/
        config[ZONE_0]->SetExtIter(iExtIter);
        
        /*--- Read in the restart file for this time step ---*/
        for (iZone = 0; iZone < nZone; iZone++) {
          
          /*--- Either instantiate the solution class or load a restart file. ---*/
          if (SolutionInstantiated == false && (iExtIter == 0 ||
              (config[ZONE_0]->GetRestart() && (iExtIter == config[ZONE_0]->GetUnst_RestartIter() ||
                                                iExtIter % config[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                                                iExtIter+1 == config[ZONE_0]->GetnExtIter())))) {
            solver[iZone] = new CBaselineSolver(geometry[iZone], config[iZone], MESH_0);
            SolutionInstantiated = true;
          }
          else
            solver[iZone]->LoadRestart(geometry, &solver, config[iZone], int(MESH_0));
        }

            if (rank == MASTER_NODE)
          cout << "Writing the volume solution for time step " << iExtIter << "." << endl;
        output->SetBaselineResult_Files(solver, geometry, config, iExtIter, nZone);
      }
      
      iExtIter++;
      if (StopCalc) break;
    }
    
  } else if (config[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL) {

	  /*--- Time-spectral simulation: merge files for each time instance (each zone). ---*/
	  unsigned short nTimeSpectral = config[ZONE_0]->GetnTimeInstances();
	  unsigned short iTimeSpectral;
	  for (iTimeSpectral = 0; iTimeSpectral < nTimeSpectral; iTimeSpectral++) {

		  /*--- Set the current instance number in the config class to "ExtIter." ---*/
		  config[ZONE_0]->SetExtIter(iTimeSpectral);

		  /*--- Read in the restart file for this time step ---*/
		  /*--- N.B. In SU2_SOL, nZone != nTimeInstances ---*/
		  for (iZone = 0; iZone < nZone; iZone++) {

			  /*--- Either instantiate the solution class or load a restart file. ---*/
			  if (iTimeSpectral == 0)
				  solver[iZone] = new CBaselineSolver(geometry[iZone], config[iZone], MESH_0);
			  else
				  solver[iZone]->LoadRestart(geometry, &solver, config[iZone], int(MESH_0));
		  }

		  /*--- Print progress in solution writing to the screen. ---*/
		  if (rank == MASTER_NODE) {
			  cout << "Writing the volume solution for time instance " << iTimeSpectral << "." << endl;
		  }

		  output->SetBaselineResult_Files(solver, geometry, config, iTimeSpectral, nZone);
	  }
  } else {

	  /*--- Steady simulation: merge the single solution file. ---*/

	  for (iZone = 0; iZone < nZone; iZone++) {
		  /*--- Definition of the solution class ---*/
		  solver[iZone] = new CBaselineSolver(geometry[iZone], config[iZone], MESH_0);
	  }

	  output->SetBaselineResult_Files(solver, geometry, config, 0, nZone);

  }
  

#ifdef HAVE_MPI
  /*--- Finalize MPI parallelization ---*/
  MPI_Finalize();
#endif
  
  /*--- End solver ---*/
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_SOL) ------------------------" << endl << endl;
  
  return EXIT_SUCCESS;
}
