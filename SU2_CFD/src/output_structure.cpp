/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
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

#include "../include/output_structure.hpp"


COutput::COutput(void) {
  
  /*--- Initialize point and connectivity counters to zero. ---*/
  nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
  nGlobal_Elem      = 0;
  nSurf_Elem        = 0;
  nGlobal_Tria      = 0;
  nGlobal_Quad      = 0;
  nGlobal_Tetr      = 0;
  nGlobal_Hexa      = 0;
  nGlobal_Wedg      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;
  
  /*--- Initialize CGNS write flag ---*/
  wrote_base_file = false;
  
  /*--- Initialize CGNS write flag ---*/
  wrote_CGNS_base = false;
  
  /*--- Initialize Tecplot write flag ---*/
  wrote_Tecplot_base = false;
  
  /*--- Initialize Paraview write flag ---*/
  wrote_Paraview_base = false;
  
}

COutput::~COutput(void) { }

void COutput::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry,
                                 CSolver *FlowSolver, unsigned long iExtIter,
                                 unsigned short val_iZone) {
  
  unsigned short iMarker;
  unsigned long iPoint, iVertex, Global_Index;
  double PressCoeff = 0.0, SkinFrictionCoeff, HeatFlux;
  double xCoord = 0.0, yCoord = 0.0, zCoord = 0.0, Mach, Pressure;
  char cstr[200];
  
  unsigned short solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();
  
#ifndef HAVE_MPI
  
  char buffer [50];
  ofstream SurfFlow_file;
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 10)   && (int(val_iZone) < 100))   sprintf (buffer, "_000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 100)  && (int(val_iZone) < 1000))  sprintf (buffer, "_00%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
    if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
  }
  else
    sprintf (buffer, ".csv");
  
  strcat (cstr, buffer);
  SurfFlow_file.precision(15);
  SurfFlow_file.open(cstr, ios::out);
  
  SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
  if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
  SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";
  
  switch (solver) {
    case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
    case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
    case TNE2_EULER: SurfFlow_file << "\"Mach_Number\"" << endl; break;
    case TNE2_NAVIER_STOKES: SurfFlow_file << "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        xCoord = geometry->node[iPoint]->GetCoord(0);
        yCoord = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) zCoord = geometry->node[iPoint]->GetCoord(2);
        
        /*--- The output should be in inches ---*/
        
        if (config->GetSystemMeasurements() == US) {
          xCoord *= 12.0; yCoord *= 12.0;
          if (nDim == 3) zCoord *= 12.0;
        }
        
        Pressure = FlowSolver->node[iPoint]->GetPressure();
        PressCoeff = FlowSolver->GetCPressure(iMarker,iVertex);
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        switch (solver) {
          case EULER :
            Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
            SurfFlow_file << scientific << Mach << endl;
            break;
          case NAVIER_STOKES: case RANS:
            SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
            HeatFlux = FlowSolver->GetHeatFlux(iMarker,iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << ", " << HeatFlux << endl;
            break;
          case TNE2_EULER:
            Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
            SurfFlow_file << scientific << Mach << endl;
            break;
          case TNE2_NAVIER_STOKES:
            SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
            HeatFlux = FlowSolver->GetHeatFlux(iMarker,iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << ", " << HeatFlux << endl;
        }
      }
    }
  }
  
  SurfFlow_file.close();
  
#else
  
  int rank, iProcessor, nProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nVertex_Surface = 0, nLocalVertex_Surface = 0;
  unsigned long MaxLocalVertex_Surface = 0;
  
  /*--- Find the max number of surface vertices among all
   partitions and set up buffers. The master node will handle the
   writing of the CSV file after gathering all of the data. ---*/
  
  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
      }
  
  /*--- Communicate the number of local vertices on each partition
   to the master node ---*/
  
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- Send and Recv buffers ---*/
  
  double *Buffer_Send_Coord_x = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_Coord_x = NULL;
  
  double *Buffer_Send_Coord_y = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_Coord_y = NULL;
  
  double *Buffer_Send_Coord_z = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_Coord_z = NULL;
  
  double *Buffer_Send_Press = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_Press = NULL;
  
  double *Buffer_Send_CPress = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_CPress = NULL;
  
  double *Buffer_Send_Mach = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_Mach = NULL;
  
  double *Buffer_Send_SkinFriction = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_SkinFriction = NULL;
  
  double *Buffer_Send_HeatTransfer = new double [MaxLocalVertex_Surface];
  double *Buffer_Recv_HeatTransfer = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Recv_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Press   = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_CPress  = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Mach    = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_SkinFriction = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_HeatTransfer = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_GlobalIndex  = new unsigned long [nProcessor*MaxLocalVertex_Surface];
  }
  
  /*--- Loop over all vertices in this partition and load the
   data of the specified type into the buffer to be sent to
   the master node. ---*/
  
  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Buffer_Send_Press[nVertex_Surface] = FlowSolver->node[iPoint]->GetPressure();
          Buffer_Send_CPress[nVertex_Surface] = FlowSolver->GetCPressure(iMarker,iVertex);
          Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
          Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
          if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2); }
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
          }
          
          Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          
          if (solver == EULER)
            Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
          if ((solver == NAVIER_STOKES) || (solver == RANS))
            Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker,iVertex);
          if (solver == TNE2_EULER)
            Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
          if (solver == TNE2_NAVIER_STOKES) {
            Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker,iVertex);
            Buffer_Send_HeatTransfer[nVertex_Surface] = FlowSolver->GetHeatFlux(iMarker,iVertex);
          }
          nVertex_Surface++;
        }
      }
  
  /*--- Send the information to the master node ---*/
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) MPI_Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Press, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Press, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (solver == EULER) MPI_Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if ((solver == NAVIER_STOKES) || (solver == RANS)) MPI_Gather(Buffer_Send_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (solver == TNE2_EULER) MPI_Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (solver == TNE2_NAVIER_STOKES) {
    MPI_Gather(Buffer_Send_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_HeatTransfer, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_HeatTransfer, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  MPI_Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks the data and writes the surface CSV file ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = config->GetSurfFlowCoeff_FileName();
    ofstream SurfFlow_file;
    
    /*--- Remove the domain number from the surface csv filename ---*/
    if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());
    
    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
      if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
      if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
      
    } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
      if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
      if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
      if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
    }
    else
      sprintf (buffer, ".csv");
    
    strcat (cstr, buffer);
    SurfFlow_file.precision(15);
    SurfFlow_file.open(cstr, ios::out);
    
    SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
    if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
    SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";
    
    switch (solver) {
      case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
      case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\"" << endl; break;
      case TNE2_EULER: SurfFlow_file << "\"Mach_Number\"" << endl; break;
      case TNE2_NAVIER_STOKES: SurfFlow_file << "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
    }
    
    /*--- Loop through all of the collected data and write each node's values ---*/
    
    unsigned long Total_Index;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
        
        /*--- Current index position and global index ---*/
        Total_Index  = iProcessor*MaxLocalVertex_Surface+iVertex;
        Global_Index = Buffer_Recv_GlobalIndex[Total_Index];
        
        /*--- Retrieve the merged data for this node ---*/
        xCoord = Buffer_Recv_Coord_x[Total_Index];
        yCoord = Buffer_Recv_Coord_y[Total_Index];
        if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
        Pressure   = Buffer_Recv_Press[Total_Index];
        PressCoeff = Buffer_Recv_CPress[Total_Index];
        
        /*--- Write the first part of the data ---*/
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        
        /*--- Write the solver-dependent part of the data ---*/
        switch (solver) {
          case EULER :
            Mach = Buffer_Recv_Mach[Total_Index];
            SurfFlow_file << scientific << Mach << endl;
            break;
          case NAVIER_STOKES: case RANS:
            SkinFrictionCoeff = Buffer_Recv_SkinFriction[Total_Index];
            SurfFlow_file << scientific << SkinFrictionCoeff << endl;
            break;
          case TNE2_EULER:
            Mach = Buffer_Recv_Mach[Total_Index];
            SurfFlow_file << scientific << Mach << endl;
            break;
          case TNE2_NAVIER_STOKES:
            SkinFrictionCoeff = Buffer_Recv_SkinFriction[Total_Index];
            SurfFlow_file << scientific << SkinFrictionCoeff << endl;
            HeatFlux = Buffer_Recv_HeatTransfer[Total_Index];
            SurfFlow_file << scientific << HeatFlux << endl;
            break;
        }
      }
    }
    
    /*--- Close the CSV file ---*/
    SurfFlow_file.close();
    
    /*--- Release the recv buffers on the master node ---*/
    
    delete [] Buffer_Recv_Coord_x;
    delete [] Buffer_Recv_Coord_y;
    if (nDim == 3) delete [] Buffer_Recv_Coord_z;
    delete [] Buffer_Recv_Press;
    delete [] Buffer_Recv_CPress;
    delete [] Buffer_Recv_Mach;
    delete [] Buffer_Recv_SkinFriction;
    delete [] Buffer_Recv_HeatTransfer;
    delete [] Buffer_Recv_GlobalIndex;
    
  }
  
  /*--- Release the memory for the remaining buffers and exit ---*/
  
  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_Press;
  delete [] Buffer_Send_CPress;
  delete [] Buffer_Send_Mach;
  delete [] Buffer_Send_SkinFriction;
  delete [] Buffer_Send_HeatTransfer;
  delete [] Buffer_Send_GlobalIndex;
  
#endif
  
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone) {
  
#ifndef HAVE_MPI
  
  unsigned long iPoint, iVertex;
  double *Solution, xCoord, yCoord, zCoord, *IntBoundary_Jump;
  unsigned short iMarker;
  char cstr[200], buffer[50];
  ofstream SurfAdj_file;
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
    if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
    if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
    
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
  }
  else
    sprintf (buffer, ".csv");
  
  strcat(cstr, buffer);
  SurfAdj_file.precision(15);
  SurfAdj_file.open(cstr, ios::out);
  
  if (geometry->GetnDim() == 2) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          IntBoundary_Jump = AdjSolver->node[iPoint]->GetIntBoundary_Jump();
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            xCoord *= 12.0;
            yCoord *= 12.0;
          }
          
          SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord << endl;
        }
    }
  }
  
  if (geometry->GetnDim() == 3) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          zCoord = geometry->node[iPoint]->GetCoord(2);
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            xCoord *= 12.0;
            yCoord *= 12.0;
            zCoord *= 12.0;
          }
          
          SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord << endl;
        }
    }
  }
  
  SurfAdj_file.close();
  
#else
  int rank, iProcessor, nProcessor;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned short nDim = geometry->GetnDim(), iMarker;
  double *Solution, *Normal, *d, *Coord;
  unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
  MaxLocalVertex_Surface = 0, nBuffer_Scalar;
  unsigned long *Buffer_Receive_nVertex = NULL;
  ofstream SurfAdj_file;
  
  /*--- Write the surface .csv file ---*/
  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
      }
  
  if (rank == MASTER_NODE)
    Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  double *Buffer_Send_Coord_x = new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Coord_y= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Coord_z= new double[MaxLocalVertex_Surface];
  unsigned long *Buffer_Send_GlobalPoint= new unsigned long[MaxLocalVertex_Surface];
  double *Buffer_Send_Sensitivity= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_PsiRho= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Phi_x= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Phi_y= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_Phi_z= new double[MaxLocalVertex_Surface];
  double *Buffer_Send_PsiE= new double[MaxLocalVertex_Surface];
  
  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Solution = AdjSolver->node[iPoint]->GetSolution();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          d = AdjSolver->node[iPoint]->GetForceProj_Vector();
          Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
          Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
          Buffer_Send_Sensitivity[nVertex_Surface] =  AdjSolver->GetCSensitivity(iMarker,iVertex);
          Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
          Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
          Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
          if (nDim == 2) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
          if (nDim == 3) {
            Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
            Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
            Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
          }
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
          }
          
          nVertex_Surface++;
        }
      }
  
  double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
  *Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
  *Buffer_Receive_PsiE = NULL;
  unsigned long *Buffer_Receive_GlobalPoint = NULL;
  
  if (rank == MASTER_NODE) {
    Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Sensitivity = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiRho = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_x = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_y = new double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Phi_z = new double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiE = new double [nProcessor*MaxLocalVertex_Surface];
  }
  
  nBuffer_Scalar = MaxLocalVertex_Surface;
  
  /*--- Send the information to the Master node ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) MPI_Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) MPI_Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node is the one who writes the surface files ---*/
  if (rank == MASTER_NODE) {
    unsigned long iVertex, GlobalPoint, position;
    char cstr[200], buffer[50];
    ofstream SurfAdj_file;
    string filename = config->GetSurfAdjCoeff_FileName();
    
    /*--- Remove the domain number from the surface csv filename ---*/
    if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());
    
    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    
    if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
      if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
      if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
      if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));
      
    } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
      if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
      if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
    }
    else
      sprintf (buffer, ".csv");
    
    strcat (cstr, buffer);
    SurfAdj_file.open(cstr, ios::out);
    SurfAdj_file.precision(15);
    
    /*--- Write the 2D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 2) {
      
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
      
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];
          
          SurfAdj_file << scientific << GlobalPoint <<
          ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
          ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
          ", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
          ", "<< Buffer_Receive_Coord_y[position]  << endl;
        }
    }
    
    /*--- Write the 3D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 3) {
      
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
      
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];
          
          SurfAdj_file << scientific << GlobalPoint <<
          ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
          ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
          ", " << Buffer_Receive_PsiE[position] <<", "<< Buffer_Receive_Coord_x[position] <<
          ", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position] << endl;
        }
    }
    
  }
  
  if (rank == MASTER_NODE) {
    delete [] Buffer_Receive_nVertex;
    delete [] Buffer_Receive_Coord_x;
    delete [] Buffer_Receive_Coord_y;
    if (nDim == 3) delete [] Buffer_Receive_Coord_z;
    delete [] Buffer_Receive_Sensitivity;
    delete [] Buffer_Receive_PsiRho;
    delete [] Buffer_Receive_Phi_x;
    delete [] Buffer_Receive_Phi_y;
    if (nDim == 3) delete [] Buffer_Receive_Phi_z;
    delete [] Buffer_Receive_PsiE;
    delete [] Buffer_Receive_GlobalPoint;
  }
  
  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_GlobalPoint;
  delete [] Buffer_Send_Sensitivity;
  delete [] Buffer_Send_PsiRho;
  delete [] Buffer_Send_Phi_x;
  delete [] Buffer_Send_Phi_y;
  delete [] Buffer_Send_Phi_z;
  delete [] Buffer_Send_PsiE;
  
  SurfAdj_file.close();
  
#endif
}

void COutput::SetSurfaceCSV_Linearized(CConfig *config, CGeometry *geometry, CSolver *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/
  
  if (!wrote_base_file) {
    
    /*--- Merge volumetric grid. ---*/
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
      cout <<"Merging volumetric triangle grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TRIANGLE    );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
      cout <<"Merging volumetric rectangle grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, RECTANGLE   );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
      cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
      cout <<"Merging volumetric hexahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Wedg != 0))
      cout <<"Merging volumetric wedge grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, WEDGE       );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
      cout <<"Merging volumetric pyramid grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, PYRAMID     );
    
    /*--- Merge surface grid. ---*/
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
      cout <<"Merging surface line grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, LINE);
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
      cout <<"Merging surface triangle grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, TRIANGLE);
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
      cout <<"Merging surface rectangle grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, RECTANGLE);
    
    /*--- Update total number of volume elements after merge. ---*/
    
    nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
    nGlobal_Hexa + nGlobal_Pyra + nGlobal_Wedg;
    
    /*--- Update total number of surface elements after merge. ---*/
    
    nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;
    
    /*--- Write the connectivity to the base binary output file, then
     clear the memory immediately for the rest of the computation. ---*/
    
    unsigned short FileFormat = config->GetOutput_FileFormat();
    if (rank == MASTER_NODE && FileFormat == CGNS_SOL) {
      SetCGNS_Connectivity(config, geometry, val_iZone);
      DeallocateConnectivity(config, geometry, false);
    }
    
  }
}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint;
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv, RecvFrom;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  nGlobal_Doma = geometry->GetnPointDomain();
  
  /*--- Allocate the coordinates data structure. ---*/
  
  Coords = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over the mesh to collect the coords of the local points. ---*/
  
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (!Local_Halo[iPoint]) {
      
      /*--- Retrieve the current coordinates at this node. ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][jPoint] = geometry->node[iPoint]->GetCoord(iDim);
        
        /*--- If US system, the output should be in inches ---*/
        
        if (config->GetSystemMeasurements() == US) {
          Coords[iDim][jPoint] *= 12.0;
        }
        
      }
      
      /*--- Increment a counter since we may be skipping over
       some halo nodes during this loop. ---*/
      
      jPoint++;
    }
  }
  
  delete [] Local_Halo;
  
#else
  
  /*--- MPI preprocessing ---*/
  int iProcessor, nProcessor, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex, iMarker;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0, periodicNodes = 0;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        RecvFrom = abs(SendRecv)-1;
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
  }
  Buffer_Send_nPoin[0] = nLocalPoint;
  
  /*--- Communicate the total number of nodes on this domain. ---*/
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
             Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  
  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  double *Buffer_Send_X = new double[MaxLocalPoint];
  double *Buffer_Recv_X = NULL;
  
  double *Buffer_Send_Y = new double[MaxLocalPoint];
  double *Buffer_Recv_Y = NULL;
  
  double *Buffer_Send_Z, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new double[MaxLocalPoint];
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_X = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new double[nGlobal_Poin];
    }
  }
  
  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  /*--- Loop over this partition to collect the coords of the local points. ---*/
  double *Coords_Local; jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halos and write only if requested ---*/
    if (!Local_Halo[iPoint] || Wrt_Halo) {
      
      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();
      
      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
      
      /*--- If US system, the output should be in inches ---*/
      
      if (config->GetSystemMeasurements() == US) {
        Buffer_Send_X[jPoint] *= 12.0;
        Buffer_Send_Y[jPoint] *= 12.0;
        if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
      }
      
      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
    }
  }
  
  /*--- Gather the coordinate data on the master node using MPI. ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    MPI_Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  MPI_Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks and sorts this variable by global index ---*/
  
  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
        Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
        Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
        if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
        jPoint++;
      }
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      jPoint = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary data buffers. ---*/
  
  delete [] Local_Halo;
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (nDim == 3) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (nDim == 3)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }
  
#endif
  
}

void COutput::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int *Conn_Elem;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
      nLocalElem = geometry->GetnElemQuad();
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      nLocalElem = geometry->GetnElemTetr();
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      nLocalElem = geometry->GetnElemHexa();
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case WEDGE:
      nLocalElem = geometry->GetnElemWedg();
      NODES_PER_ELEMENT = N_POINTS_WEDGE;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(EXIT_FAILURE); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex;
  int SendRecv, RecvFrom;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Check if this is a halo node. ---*/
      isHalo = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        if (Local_Halo[iPoint])
          isHalo = true;
      }
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the temporary array. Do not merge the
       halo cells (periodic BC). Note that we are adding one to
       the index value because CGNS/Tecplot use 1-based indexing. ---*/
      
      if (!isHalo) {
        nElem_Total++;
        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
          Conn_Elem[jNode] = (int)geometry->elem[iElem]->GetNode(iNode) + 1;
          
          /*--- Increment jNode as the counter. ---*/
          jNode++;
        }
      }
    }
  }
  
  delete [] Local_Halo;
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int iProcessor, jProcessor, nProcessor;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex, iMarker;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0, pElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem, notPeriodic, notHalo, addedPeriodic;
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  int *Buffer_Send_Elem = new int[nBuffer_Scalar];
  int *Buffer_Recv_Elem = NULL;
  
  int *Buffer_Send_Halo = new int[MaxLocalElem];
  int *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new int[nProcessor*nBuffer_Scalar];
    Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
    Conn_Elem = new int[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
            (SendRecv < 0)) {
          Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[nProcessor];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[nProcessor*maxAddedPeriodic];
  
  for(iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                MPI_COMM_WORLD);
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain ---*/
        notPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                       (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                       (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  
  jNode = 0; jElem = 0;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the send buffer. ---*/
      
      Buffer_Send_Halo[jElem] = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        
        /*--- Store the global index values directly. ---*/
        
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        Buffer_Send_Elem[jNode] = (int)geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/
        
        if (Local_Halo[iPoint]) {
          Buffer_Send_Halo[jElem] = true;
        }
        
        /*--- Increment jNode as the counter. We need this because iElem
         may include other elements that we skip over during this loop. ---*/
        
        jNode++;
      }
      jElem++;
    }
  }
  
  /*--- Gather the element connectivity information. ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_INT, Buffer_Recv_Elem, nBuffer_Scalar, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Halo, MaxLocalElem, MPI_INT, Buffer_Recv_Halo, MaxLocalElem, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[nProcessor*MaxLocalElem];
    for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
          
          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;
          
        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }
    
    /*--- Store the unique connectivity list for this element type. ---*/
    
    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
#endif
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case RECTANGLE:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case WEDGE:
        nGlobal_Wedg = nElem_Total;
        if (nGlobal_Wedg > 0) Conn_Wedg = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(EXIT_FAILURE); break;
    }
  }
  
}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned short iMarker;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int *Conn_Elem;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  nLocalElem = 0;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          nLocalElem++;
        }
      }
    }
  }
  
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(EXIT_FAILURE); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  unsigned long iVertex;
  int SendRecv, RecvFrom;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Check if this is a halo node. ---*/
          isHalo = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            if (Local_Halo[iPoint])
              isHalo = true;
          }
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the temporary array. Do not merge any
           halo cells (periodic BC). Note that we are adding one to
           the index value because CGNS/Tecplot use 1-based indexing. ---*/
          if (!isHalo) {
            nElem_Total++;
            for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
              Conn_Elem[jNode] = (int)geometry->bound[iMarker][iElem]->GetNode(iNode) + 1;
              
              /*--- Increment jNode as the counter. ---*/
              jNode++;
            }
          }
        }
      }
  
  delete [] Local_Halo;
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int iProcessor, jProcessor, nProcessor;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0, pElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem, notPeriodic, notHalo, addedPeriodic;
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  int *Buffer_Send_Elem = new int[nBuffer_Scalar];
  int *Buffer_Recv_Elem = NULL;
  
  int *Buffer_Send_Halo = new int[MaxLocalElem];
  int *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new int[nProcessor*nBuffer_Scalar];
    Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
    Conn_Elem = new int[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
            (SendRecv < 0)) {
          Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[nProcessor];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[nProcessor*maxAddedPeriodic];
  
  for(iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                MPI_COMM_WORLD);
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain ---*/
        notPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                       (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                       (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  jNode = 0; jElem = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for(iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if(geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the send buffer. ---*/
          
          Buffer_Send_Halo[jElem] = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            
            /*--- Store the global index values directly. ---*/
            
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            Buffer_Send_Elem[jNode] = (int)geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/
            
            if (Local_Halo[iPoint])
              Buffer_Send_Halo[jElem] = true;
            
            /*--- Increment jNode as the counter. We need this because iElem
             may include other elements that we skip over during this loop. ---*/
            
            jNode++;
          }
          jElem++;
        }
      }
  
  /*--- Gather the element connectivity information. ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_INT, Buffer_Recv_Elem, nBuffer_Scalar, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Halo, MaxLocalElem, MPI_INT, Buffer_Recv_Halo, MaxLocalElem, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[nProcessor*MaxLocalElem];
    for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
          
          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;
          
        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }
    
    /*--- Store the unique connectivity list for this element type. ---*/
    
    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
#endif
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case RECTANGLE:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(EXIT_FAILURE); break;
    }
  }
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, jVar = 0, iSpecies, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0, iVar_Eddy = 0, iVar_Sharp = 0;
  unsigned short iVar_GridVel = 0, iVar_PressCp = 0, iVar_Density = 0, iVar_Lam = 0, iVar_MachMean = 0,
  iVar_Tempv = 0, iVar_EF =0, iVar_Temp = 0, iVar_Mach = 0, iVar_Press = 0, iVar_TempLam = 0,
  iVar_ViscCoeffs = 0, iVar_Sens = 0, iVar_FEA = 0, iVar_Extra = 0;
  
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;
  
  double *Aux_Frict = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL, *Aux_Sens = NULL;
  
  bool grid_movement  = (config->GetGrid_Movement());
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool transition     = (config->GetKind_Trans_Model() == LM);
  bool flow           = (( config->GetKind_Solver() == EULER             ) ||
                         ( config->GetKind_Solver() == NAVIER_STOKES     ) ||
                         ( config->GetKind_Solver() == RANS              ) ||
                         ( config->GetKind_Solver() == ADJ_EULER         ) ||
                         ( config->GetKind_Solver() == ADJ_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == ADJ_RANS          )   );
  
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  double RefAreaCoeff = config->GetRefAreaCoeff();
  double Gamma = config->GetGamma();
  double RefVel2;
  
  /*--- Set the non-dimensionalization ---*/
  if (flow) {
    if (grid_movement) {
      Gas_Constant = config->GetGas_ConstantND();
      Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
    }
    RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
    RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
    factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  }
  
  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/
  switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case RANS :
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL;
      if (transition) ThirdIndex=TRANS_SOL;
      else ThirdIndex = NONE;
      break;
    case TNE2_EULER : case TNE2_NAVIER_STOKES:
      FirstIndex = TNE2_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case POISSON_EQUATION:
      FirstIndex = POISSON_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case WAVE_EQUATION:
      FirstIndex = WAVE_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case HEAT_EQUATION:
      FirstIndex = HEAT_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case LINEAR_ELASTICITY:
      FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES :
      FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_TNE2_EULER : case ADJ_TNE2_NAVIER_STOKES :
      FirstIndex = ADJTNE2_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_RANS :
      FirstIndex = ADJFLOW_SOL;
      if (config->GetFrozen_Visc()) SecondIndex = NONE;
      else SecondIndex = ADJTURB_SOL;
      ThirdIndex = NONE;
      break;
    case LIN_EULER : case LIN_NAVIER_STOKES : ThirdIndex = NONE;
      FirstIndex = LINFLOW_SOL; SecondIndex = NONE;
      break;
    default: SecondIndex = NONE; ThirdIndex = NONE;
      break;
  }
  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  
  nVar_Total = nVar_Consv;
  
  /*--- Add the limiters ---*/
  if (config->GetWrt_Limiters()) nVar_Total += nVar_Consv;
  
  /*--- Add the residuals ---*/
  if (config->GetWrt_Residuals()) nVar_Total += nVar_Consv;
  
  /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
  if (grid_movement) {
    iVar_GridVel = nVar_Total;
    if (geometry->GetnDim() == 2) nVar_Total += 2;
    else if (geometry->GetnDim() == 3) nVar_Total += 3;
  }
  
  if ((config->GetKind_Regime() == FREESURFACE)) {
    /*--- Density ---*/
    iVar_Density = nVar_Total;
    nVar_Total += 1;
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Pressure, Temperature, Cp, Mach ---*/
    iVar_PressCp = nVar_Total;
    nVar_Total += 3;
    iVar_MachMean = nVar_Total;
    nVar_Total += 1;
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Laminar Viscosity ---*/
    iVar_Lam = nVar_Total;
    nVar_Total += 1;
    /*--- Skin Friction, Heat Flux, & yPlus ---*/
    iVar_ViscCoeffs = nVar_Total;
    nVar_Total += 3;
  }
  
  if (Kind_Solver == RANS) {
    /*--- Eddy Viscosity ---*/
    iVar_Eddy = nVar_Total;
    nVar_Total += 1;
  }
  
  if (config->GetWrt_SharpEdges()) {
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      /*--- Sharp edges ---*/
      iVar_Sharp = nVar_Total;
      nVar_Total += 1;
    }
  }
  
  if ((Kind_Solver == TNE2_EULER)         ||
      (Kind_Solver == TNE2_NAVIER_STOKES)   ) {
    /*--- Mach ---*/
    iVar_Mach = nVar_Total;
    nVar_Total++;
    /*--- Pressure ---*/
    iVar_Press = nVar_Total;
    nVar_Total++;
    /*--- Temperature ---*/
    iVar_Temp = nVar_Total;
    nVar_Total++;
    /*--- Vib-El. Temperature ---*/
    iVar_Tempv = nVar_Total;
    nVar_Total++;
  }
  
  if (Kind_Solver == TNE2_NAVIER_STOKES) {
    /*--- Diffusivity, viscosity, & thermal conductivity ---*/
    iVar_TempLam = nVar_Total;
    nVar_Total += config->GetnSpecies()+3;
  }
  
  if (Kind_Solver == POISSON_EQUATION) {
    iVar_EF = geometry->GetnDim();
    nVar_Total += geometry->GetnDim();
  }
  
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == ADJ_TNE2_EULER         ) ||
      ( Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
    /*--- Surface sensitivity coefficient, and solution sensor ---*/
    iVar_Sens   = nVar_Total;
    nVar_Total += 2;
  }
  
  if (Kind_Solver == LINEAR_ELASTICITY) {
    /*--- Surface sensitivity coefficient, and solution sensor ---*/
    iVar_FEA   = nVar_Total;
    nVar_Total += 2;
  }
  
  if (config->GetExtraOutput()) {
    if (Kind_Solver == RANS) {
      iVar_Extra  = nVar_Total;
      nVar_Extra  = solver[TURB_SOL]->GetnOutputVariables();
      nVar_Total += nVar_Extra;
    }
    if ((Kind_Solver == TNE2_EULER)         ||
        (Kind_Solver == TNE2_NAVIER_STOKES)) {
      iVar_Extra  = nVar_Total;
      nVar_Extra  = solver[TNE2_SOL]->GetnVar();
      nVar_Total += nVar_Extra;
    }
  }
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  
  unsigned long nTotalPoints = 0;
  int SendRecv, RecvFrom;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- In case there is grid movement ---*/
  double *Grid_Vel;
  
  /*--- First, loop through the mesh in order to find and store the
   value of some coefficients at any surface nodes. They
   will be placed in an auxiliary vector and then communicated like
   all other volumetric variables. ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    Aux_Frict = new double [geometry->GetnPoint()];
    Aux_Heat  = new double [geometry->GetnPoint()];
    Aux_yPlus = new double [geometry->GetnPoint()];
    
    for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatFlux(iMarker,iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
        }
      }
  }
  
  if ((Kind_Solver == ADJ_EULER)         ||
      (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)            ) {
    
    Aux_Sens = new double [geometry->GetnPoint()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
  }
  
  if ((Kind_Solver == ADJ_TNE2_EULER)         ||
      (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)   ) {
    
    Aux_Sens = new double [geometry->GetnPoint()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJTNE2_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (includes original periodic nodes). ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halo nodes & only write if requested ---*/
    
    if (!Local_Halo[iPoint]) {
      
      /*--- Solution (first, second and third system of equations) ---*/
      jVar = 0;
      for (iVar = 0; iVar < nVar_First; iVar++) {
        Data[jVar][jPoint] = solver[FirstIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Second; iVar++) {
        Data[jVar][jPoint] = solver[SecondIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Third; iVar++) {
        Data[jVar][jPoint] = solver[ThirdIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      /*--- Limiters (first, second and third system of equations) ---*/
      if (config->GetWrt_Limiters()) {
        
        if (solver[FirstIndex]->node[iPoint]->GetLimiter_Primitive() != NULL) {
          for (iVar = 0; iVar < nVar_First; iVar++) {
            Data[jVar][jPoint] = solver[FirstIndex]->node[iPoint]->GetLimiter_Primitive(iVar);
            jVar++;
          }
        }
        else { for (iVar = 0; iVar < nVar_First; iVar++) { Data[jVar][jPoint] = 0.0; jVar++; } }
        
        if (solver[SecondIndex]->node[iPoint]->GetLimiter() != NULL) {
          for (iVar = 0; iVar < nVar_Second; iVar++) {
            Data[jVar][jPoint] = solver[SecondIndex]->node[iPoint]->GetLimiter(iVar);
            jVar++;
          }
        }
        else { for (iVar = 0; iVar < nVar_Second; iVar++) { Data[jVar][jPoint] = 0.0; jVar++; } }
        
        if (solver[ThirdIndex]->node[iPoint]->GetLimiter() != NULL) {
          for (iVar = 0; iVar < nVar_Third; iVar++) {
            Data[jVar][jPoint] = solver[ThirdIndex]->node[iPoint]->GetLimiter(iVar);
            jVar++;
          }
        }
        else { for (iVar = 0; iVar < nVar_Third; iVar++) { Data[jVar][jPoint] = 0.0; jVar++; } }
        
      }
      
      /*--- Residual (first, second and third system of equations) ---*/
      if (config->GetWrt_Residuals()) {
        for (iVar = 0; iVar < nVar_First; iVar++) {
          Data[jVar][jPoint] = solver[FirstIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Second; iVar++) {
          Data[jVar][jPoint] = solver[SecondIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Third; iVar++) {
          Data[jVar][jPoint] = solver[ThirdIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
      }
      
      /*--- For unsteady problems with grid movement, write the mesh velocities ---*/
      if (grid_movement) {
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
          Data[jVar][jPoint] = Grid_Vel[iDim];
          jVar++;
        }
      }
      
      /*--- Any extra data for output files ---*/
      switch (Kind_Solver) {
        case EULER:
          
          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
          if (compressible) {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
          }
          if (incompressible || freesurface) {
            if (freesurface) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressureInc(); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
          }
          if (config->GetWrt_SharpEdges()) {
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          }
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus ---*/
        case NAVIER_STOKES:
          if (compressible) {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          if (incompressible || freesurface) {
            if (freesurface) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressureInc(); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          }
          if (config->GetWrt_SharpEdges()) {
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          }
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus, eddy viscosity ---*/
        case RANS:
          if (compressible) {
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); jVar++;
          }
          if (incompressible || freesurface) {
            if (freesurface) {
              Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc(); jVar++;
            }
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressureInc(); jVar++;
            Data[jVar][jPoint] = 0.0; jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosityInc(); jVar++;
          }
          if (config->GetWrt_SharpEdges()) {
            Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          }
          break;
          /*--- Write poisson field. ---*/
        case POISSON_EQUATION:
          for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
            Data[jVar][jPoint] = -1.0*solver[POISSON_SOL]->node[iPoint]->GetGradient(0,iDim);
            jVar++;
          }
          break;
          
        case TNE2_EULER:
          /*--- Write Mach number ---*/
          Data[jVar][jPoint] = sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2())
          / solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
          jVar++;
          /*--- Write Pressure ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
          jVar++;
          /*--- Write Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
          jVar++;
          /*--- Write Vib.-El. Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
          jVar++;
          break;
          
        case TNE2_NAVIER_STOKES:
          /*--- Write Mach number ---*/
          Data[jVar][jPoint] = sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2())
          / solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
          jVar++;
          /*--- Write Pressure ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
          jVar++;
          /*--- Write Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
          jVar++;
          /*--- Write Vib.-El. Temperature ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
          jVar++;
          /*--- Write species diffusion coefficients ---*/
          for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
            Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff()[iSpecies];
            jVar++;
          }
          /*--- Write viscosity ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
          jVar++;
          /*--- Write thermal conductivity ---*/
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
          jVar++;
          Data[jVar][jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
          break;
          
        case ADJ_EULER:      case ADJ_NAVIER_STOKES:     case ADJ_RANS:
          
          Data[jVar][jPoint] = Aux_Sens[iPoint]; jVar++;
          if (config->GetKind_ConvNumScheme() == SPACE_CENTERED)
          { Data[jVar][jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(); jVar++; }
          if (config->GetKind_ConvNumScheme() == SPACE_UPWIND)
          { Data[jVar][jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0); jVar++; }
          break;
          
        case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
          Data[jVar][jPoint] = Aux_Sens[iPoint]; jVar++;
          if (config->GetKind_ConvNumScheme() == SPACE_CENTERED)
          { Data[jVar][jPoint] = solver[ADJTNE2_SOL]->node[iPoint]->GetSensor(); jVar++; }
          if (config->GetKind_ConvNumScheme() == SPACE_UPWIND)
          { Data[jVar][jPoint] = solver[ADJTNE2_SOL]->node[iPoint]->GetLimiter(0); jVar++; }
          break;
          
        case LINEAR_ELASTICITY:
          Data[jVar][jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
          jVar++;
          Data[jVar][jPoint] = solver[FEA_SOL]->node[iPoint]->GetFlow_Pressure();
          jVar++;
          break;
          
      }
    }
    
    if (config->GetExtraOutput()) {
      if (Kind_Solver == RANS) {
        for (unsigned short iVar = 0; iVar < nVar_Extra; iVar++) {
          Data[jVar][jPoint] =  solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+iVar];
          jVar++;
        }
      }
      if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
        for (unsigned short iVar = 0; iVar < nVar_Extra; iVar++) {
          Data[jVar][jPoint] =  solver[TNE2_SOL]->OutputVariables[iPoint*nVar_Extra+iVar];
          jVar++;
        }
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  int rank, iProcessor, nProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  /*--- Local variables needed for merging with MPI ---*/
  unsigned short CurrentIndex;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        RecvFrom = abs(SendRecv)-1;
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
    
  }
  Buffer_Send_nPoint[0] = nLocalPoint;
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  double *Buffer_Send_Var = new double[MaxLocalPoint];
  double *Buffer_Recv_Var = NULL;
  
  double *Buffer_Send_Res = new double[MaxLocalPoint];
  double *Buffer_Recv_Res = NULL;
  
  double *Buffer_Send_Vol = new double[MaxLocalPoint];
  double *Buffer_Recv_Vol = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Auxiliary vectors for surface coefficients ---*/
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    Aux_Frict = new double[geometry->GetnPoint()];
    Aux_Heat  = new double[geometry->GetnPoint()];
    Aux_yPlus = new double[geometry->GetnPoint()];
  }
  
  if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)  || (Kind_Solver == ADJ_TNE2_EULER)    ||
      (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)                            ) {
    Aux_Sens = new double[geometry->GetnPoint()];
  }
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Res = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Vol = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new double[nGlobal_Poin];
    }
  }
  
  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    
    /*--- Logic for which solution class to draw from. ---*/
    
    jVar = iVar;
    CurrentIndex = FirstIndex;
    if ((SecondIndex != NONE) && (iVar > nVar_First-1)) {
      jVar = iVar - nVar_First;
      CurrentIndex = SecondIndex;
    }
    if ((SecondIndex != NONE) && (ThirdIndex != NONE) && (iVar > (nVar_First + nVar_Second-1))) {
      jVar = iVar - nVar_First - nVar_Second;
      CurrentIndex = ThirdIndex;
    }
    
    /*--- Loop over this partition to collect the current variable ---*/
    
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        
        Buffer_Send_Var[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar);
        
        if (config->GetWrt_Limiters()) {
          Buffer_Send_Vol[jPoint] = solver[CurrentIndex]->node[iPoint]->GetLimiter_Primitive(jVar);
        }
        
        if (config->GetWrt_Residuals()) {
          Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes.GetBlock(iPoint, jVar);
        }
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        
        jPoint++;
        
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    if (config->GetWrt_Limiters()) {
      MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    }
    
    if (config->GetWrt_Residuals()) {
      MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    }
    
    if (iVar == 0) {
      MPI_Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          
          if (config->GetWrt_Limiters()) {
            Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
          }
          
          if (config->GetWrt_Residuals()) {
            unsigned short ExtraIndex;
            ExtraIndex = nVar_Consv;
            if (config->GetWrt_Limiters()) ExtraIndex = 2*nVar_Consv;
            Data[iVar+ExtraIndex][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          }
          
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
  }
  
  /*--- Additional communication routine for the grid velocity. Note that
   we are reusing the same temporary buffers from above for efficiency.
   Also, in the future more routines like this could be used to write
   an arbitrary number of additional variables to the file. ---*/
  
  if (grid_movement) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0; double *Grid_Vel;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the three grid velocity components. ---*/
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        Buffer_Send_Var[jPoint] = Grid_Vel[0];
        Buffer_Send_Res[jPoint] = Grid_Vel[1];
        if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (geometry->GetnDim() == 3) {
      MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_GridVel;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          if (geometry->GetnDim() == 3)
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the Density in Free-surface problems ---*/
  if (config->GetKind_Regime() == FREESURFACE) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Density;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
  }
  
  /*--- Communicate Pressure, Cp, and Mach ---*/
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the coefficient of pressure at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure();
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature();
          Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff;
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressureInc();
          Buffer_Send_Res[jPoint] = 0.0;
          Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff;
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_PressCp;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate Mach---*/
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
          solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/
          sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref()));
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_MachMean;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  /*--- Laminar Viscosity ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        if (compressible) {
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        }
        if (incompressible || freesurface) {
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Lam;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate skin friction, heat transfer, y+ ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the viscous coefficients at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatFlux(iMarker,iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
        }
      }
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
          Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
          Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_ViscCoeffs;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the Eddy Viscosity ---*/
  if (Kind_Solver == RANS) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        if (compressible) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        }
        if (incompressible || freesurface) {
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosityInc();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Eddy;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
  }
  
  /*--- Communicate the Sharp Edges ---*/
  if (config->GetWrt_SharpEdges()) {
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the pressure and mach variables. ---*/
          Buffer_Send_Var[jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance();
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Sharp;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            jPoint++;
          }
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
  }
  
  /*--- Communicate additional variables for the two-temperature solvers ---*/
  if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
    
    /*--- Mach number ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] =
        sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2())
        /solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Mach;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
    /*--- Pressure ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Press;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
    /*--- Temperature ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Temp;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index             = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
    /*--- Vib-el Temperature ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_Tempv;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index             = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  if (Kind_Solver == TNE2_NAVIER_STOKES) {
    /*--- Species diffusion coefficients ---*/
    // Loop over this partition to collect the current variable
    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the Mach number variables. ---*/
          Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff()[iSpecies];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      if (rank == MASTER_NODE) {
        jPoint = 0;
        iVar = iVar_TempLam+iSpecies;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            jPoint++;
          }
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    /*--- Laminar viscosity ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_TempLam+config->GetnSpecies();
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
    /*--- Thermal conductivity ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_TempLam+config->GetnSpecies()+1;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
    /*--- Vib-el Thermal conductivity ---*/
    // Loop over this partition to collect the current variable
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the Mach number variables. ---*/
        Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      iVar = iVar_TempLam+config->GetnSpecies()+2;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the surface sensitivity ---*/
  
  if ((Kind_Solver == ADJ_EULER)         ||
      (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)          ||
      (Kind_Solver == ADJ_TNE2_EULER)    ||
      (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)) {
    
    /*--- First, loop through the mesh in order to find and store the
     value of the surface sensitivity at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker,iVertex);
        }
      }
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
        
        Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
        if (config->GetKind_ConvNumScheme() == SPACE_CENTERED)
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(iPoint);
        if (config->GetKind_ConvNumScheme() == SPACE_UPWIND)
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0);
        
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Sens;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the Linear elasticity ---*/
  if ( Kind_Solver == LINEAR_ELASTICITY ) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        Buffer_Send_Var[jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
        Buffer_Send_Res[jPoint] = solver[FEA_SOL]->node[iPoint]->GetFlow_Pressure();
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_FEA;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  
  if (config->GetExtraOutput()) {
    
    for (jVar = 0; jVar < nVar_Extra; jVar++) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Get this variable into the temporary send buffer. ---*/
          if (Kind_Solver == RANS) {
            Buffer_Send_Var[jPoint] = solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
          }
          if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
          }
          jPoint++;
          
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Extra;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+jVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            jPoint++;
          }
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
  }
  
  
  /*--- Immediately release the temporary buffers. ---*/
  
  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_Res;
  delete [] Buffer_Send_Vol;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_Res;
    delete [] Buffer_Recv_Vol;
    delete [] Buffer_Recv_GlobalIndex;
  }
  
#endif
  
  /*--- Release memory needed for surface coefficients ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    delete [] Aux_Frict; delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == ADJ_TNE2_EULER         ) ||
      ( Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
    delete [] Aux_Sens;
  }
  
}

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;
  
  nVar_Total = config->fields.size() - 1;
  
  /*--- Merge the solution either in serial or parallel. ---*/
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv, RecvFrom;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (!Local_Halo[iPoint]) {
      /*--- Solution (first, and second system of equations) ---*/
      unsigned short jVar = 0;
      for (iVar = 0; iVar < nVar_Total; iVar++) {
        Data[jVar][jPoint] = solver->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  int rank, nProcessor, iProcessor;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  /*--- Local variables needed for merging with MPI ---*/
  
  unsigned long iVertex, iMarker;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        RecvFrom = abs(SendRecv)-1;
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
    
  }
  Buffer_Send_nPoint[0] = nLocalPoint;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  double *Buffer_Send_Var = new double[MaxLocalPoint];
  double *Buffer_Recv_Var = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new double[nGlobal_Poin];
    }
    
  }
  
  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos and write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        Buffer_Send_Var[jPoint] = solver->node[iPoint]->GetSolution(iVar);
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (iVar == 0) {
      MPI_Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  
  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }
  
#endif
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  ofstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
  if (config->GetAdjoint()) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetRestart_FlowFileName();
  }
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, int(val_iZone));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, int(iExtIter));
  }
  
  /*--- Open the restart file and write the solution. ---*/
  restart_file.open(filename.c_str(), ios::out);
  restart_file.precision(15);
  
  /*--- Write the header line based on the particular solver ----*/
  restart_file << "\"PointID\"";
  
  /*--- Mesh coordinates are always written to the restart first ---*/
  if (nDim == 2) {
    restart_file << "\t\"x\"\t\"y\"";
  } else {
    restart_file << "\t\"x\"\t\"y\"\t\"z\"";
  }
  
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    restart_file << "\t\"Conservative_" << iVar+1<<"\"";
  }
  if (config->GetWrt_Limiters()) {
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      restart_file << "\t\"Limiter_" << iVar+1<<"\"";
    }
  }
  if (config->GetWrt_Residuals()) {
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      restart_file << "\t\"Residual_" << iVar+1<<"\"";
    }
  }
  
  /*--- Mesh velocities for dynamic mesh cases ---*/
  if (grid_movement) {
    if (nDim == 2) {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"";
    } else {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"\t\"Grid_Velz\"";
    }
  }
  
  /*--- Solver specific output variables ---*/
  if (config->GetKind_Regime() == FREESURFACE) {
    restart_file << "\t\"Density\"";
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"Pressure_Coefficient\"\t\"Mach\"";
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient\"\t\"Heat_Flux\"\t\"Y_Plus\"";
  }
  
  if (Kind_Solver == RANS) {
    restart_file << "\t\"Eddy_Viscosity\"";
  }
  
  if (config->GetWrt_SharpEdges()) {
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      restart_file << "\t\"Sharp_Edge_Dist\"";
    }
  }
  
  if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) {
    restart_file << "\t\"Mach\"\t\"Pressure\"\t\"Temperature\"\t\"Temperature_ve\"";
  }
  
  if (Kind_Solver == TNE2_NAVIER_STOKES) {
    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
      restart_file << "\t\"DiffusionCoeff_" << iSpecies << "\"";
    restart_file << "\t\"Laminar_Viscosity\"\t\"ThermConductivity\"\t\"ThermConductivity_ve\"";
  }
  
  if (Kind_Solver == POISSON_EQUATION) {
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
      restart_file << "\t\"poissonField_" << iDim+1 << "\"";
  }
  
  if ((Kind_Solver == ADJ_EULER              ) ||
      (Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      (Kind_Solver == ADJ_RANS               ) ||
      (Kind_Solver == ADJ_TNE2_EULER         ) ||
      (Kind_Solver == ADJ_TNE2_NAVIER_STOKES )   ) {
    restart_file << "\t\"Surface_Sensitivity\"\t\"Solution_Sensor\"";
  }
  
  if (Kind_Solver == LINEAR_ELASTICITY) {
    restart_file << "\t\"Von_Mises_Stress\"\t\"Flow_Pressure\"";
  }
  
  if (config->GetExtraOutput()) {
    string *headings = NULL;
    //if (Kind_Solver == RANS){
    headings = solver[TURB_SOL]->OutputHeadingNames;
    //}
    
    for (iVar = 0; iVar < nVar_Extra; iVar++) {
      if (headings == NULL){
        restart_file << "\t\"ExtraOutput_" << iVar+1<<"\"";
      }else{
        restart_file << "\t\""<< headings[iVar] <<"\"";
      }
    }
  }
  
  restart_file << endl;
  
  /*--- Write the restart file ---*/
  
  for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
    
    /*--- Index of the point ---*/
    restart_file << iPoint << "\t";
    
    /*--- Write the grid coordinates first ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      restart_file << scientific << Coords[iDim][iPoint] << "\t";
    }
    
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      restart_file << scientific << Data[iVar][iPoint] << "\t";
    }
    restart_file << endl;
  }
  
  restart_file.close();
  
}

void COutput::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables and initialization ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for coordinate data ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      delete [] Coords[iDim];
    }
    delete [] Coords;
    
  }
}

void COutput::DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for connectivity data ---*/
    if (surf_sol) {
      if (nGlobal_Line > 0) delete [] Conn_Line;
      if (nGlobal_BoundTria > 0) delete [] Conn_BoundTria;
      if (nGlobal_BoundQuad > 0) delete [] Conn_BoundQuad;
    }
    else {
      if (nGlobal_Tria > 0) delete [] Conn_Tria;
      if (nGlobal_Quad > 0) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0) delete [] Conn_Hexa;
      if (nGlobal_Wedg > 0) delete [] Conn_Wedg;
      if (nGlobal_Pyra > 0) delete [] Conn_Pyra;
    }
    
  }
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for solution data ---*/
    for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) {
      delete [] Data[iVar];
    }
    delete [] Data;
    
  }
}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
  char cstr[200], buffer[50], turb_resid[1000];
  unsigned short iMarker, iMarker_Monitoring, iSpecies;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff;
  
  bool rotating_frame = config->GetRotating_Frame();
  bool aeroelastic    = config->GetAeroelastic_Simulation();
  bool equiv_area     = config->GetEquivArea();
  bool turbulent      = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS));
  bool transition     = (config->GetKind_Trans_Model() == LM);
  bool frozen_turb    = config->GetFrozen_Visc();
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool inv_design     = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());
  bool output_1d      = config->GetWrt_1D_Output();
  bool output_per_surface = false;
  if(config->GetnMarker_Monitoring() > 1) output_per_surface = true;
  
  bool isothermal = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC   ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)   )
      isothermal = true;
  
  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  strcpy (cstr, filename.data());
  
  if (config->GetWrt_Unsteady() && config->GetRestart()) {
    long iExtIter = config->GetUnst_RestartIter();
    if (int(iExtIter) < 10) sprintf (buffer, "_0000%d", int(iExtIter));
    if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d", int(iExtIter));
    if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d", int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d", int(iExtIter));
    if (int(iExtIter) >= 10000) sprintf (buffer, "_%d", int(iExtIter));
    strcat(cstr,buffer);
  }
  
  if (config->GetOutput_FileFormat() == TECPLOT) sprintf (buffer, ".dat");
  else if (config->GetOutput_FileFormat() == TECPLOT_BINARY)  sprintf (buffer, ".plt");
  else if ((config->GetOutput_FileFormat() == CGNS_SOL) ||
           (config->GetOutput_FileFormat() == PARAVIEW))  sprintf (buffer, ".csv");
  strcat(cstr,buffer);
  
  ConvHist_file->open(cstr, ios::out);
  ConvHist_file->precision(15);
  
  /*--- Begin of the header ---*/
  
  char begin[]= "\"Iteration\"";
  
  /*--- Header for the coefficients ---*/
  
  char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
  char heat_coeff[]= ",\"HeatFlux_Total\",\"HeatFlux_Maximum\"";
  char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
  char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
  char free_surface_coeff[]= ",\"CFreeSurface\"";
  char wave_coeff[]= ",\"CWave\"";
  char fea_coeff[]= ",\"CFEA\"";
  char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
  char oneD_outputs[]= ",\"Avg_TotalPress\",\"Avg_Mach\",\"Avg_Temperature\",\"MassFlowRate\",\"FluxAvg_Pressure\",\"FluxAvg_Density\",\"FluxAvg_Velocity\",\"FluxAvg_Enthalpy\"";
  char Cp_inverse_design[]= ",\"Cp_Diff\"";
  char Heat_inverse_design[]= ",\"HeatFlux_Diff\"";
  
  /* Find the markers being monitored and create a header for them */
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
    monitoring_coeff += ",\"CLift_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CDrag_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CSideForce_" + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CFz_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMz_"    + Monitoring_Tag + "\"";
    aeroelastic_coeff += ",\"plunge_" + Monitoring_Tag + "\"";
    aeroelastic_coeff += ",\"pitch_"  + Monitoring_Tag + "\"";
  }
  
  /*--- Header for the residuals ---*/
  
  char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
  char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";
  switch (config->GetKind_Turb_Model()) {
    case SA:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case ML:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case SST:	sprintf (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
  }
  char trans_resid[] = ",\"Res_Trans[0]\",\"Res_Trans[1]\"";
  char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";
  char levelset_resid[]= ",\"Res_LevelSet\"";
  char adj_levelset_resid[]= ",\"Res_AdjLevelSet\"";
  char wave_resid[]= ",\"Res_Wave[0]\",\"Res_Wave[1]\"";
  char fea_resid[]= ",\"Res_FEA\"";
  char heat_resid[]= ",\"Res_Heat\"";
  
  /*--- End of the header ---*/
  
  char end[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY)) {
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }
  
  /*--- Write the header, case depending ---*/
  switch (config->GetKind_Solver()) {
      
    case EULER : case NAVIER_STOKES: case RANS :
      ConvHist_file[0] << begin << flow_coeff;
      if (isothermal) ConvHist_file[0] << heat_coeff;
      if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
      if (inv_design) {
        ConvHist_file[0] << Cp_inverse_design;
        if (isothermal) ConvHist_file[0] << Heat_inverse_design;
      }
      if (rotating_frame) ConvHist_file[0] << rotating_frame_coeff;
      ConvHist_file[0] << flow_resid;
      if (turbulent) ConvHist_file[0] << turb_resid;
      if (transition) ConvHist_file[0] << trans_resid;
      if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
      if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
      if (output_1d) ConvHist_file[0] << oneD_outputs;
      ConvHist_file[0] << end;
      if (freesurface) {
        ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
        ConvHist_file[0] << flow_resid << levelset_resid << end;
      }
      break;
      
    case TNE2_EULER : case TNE2_NAVIER_STOKES:
      ConvHist_file[0] << begin << flow_coeff;
      if (config->GetKind_Solver() == TNE2_NAVIER_STOKES) {
        ConvHist_file[0] << heat_coeff;
        if (inv_design) ConvHist_file[0] << Heat_inverse_design;
      }
      for (iSpecies = 0; iSpecies < config->GetnSpecies()+5; iSpecies++)
        ConvHist_file[0] << ",\"Residual[" << iSpecies << "]\"";
      ConvHist_file[0] << end;
      break;
      
    case ADJ_EULER      : case ADJ_NAVIER_STOKES      : case ADJ_RANS:
    case ADJ_TNE2_EULER : case ADJ_TNE2_NAVIER_STOKES :
      ConvHist_file[0] << begin << adj_coeff << adj_flow_resid;
      if ((turbulent) && (!frozen_turb)) ConvHist_file[0] << adj_turb_resid;
      ConvHist_file[0] << end;
      if (freesurface) {
        ConvHist_file[0] << begin << adj_coeff << adj_flow_resid << adj_levelset_resid << end;
      }
      break;
      
    case WAVE_EQUATION:
      ConvHist_file[0] << begin << wave_coeff;
      ConvHist_file[0] << wave_resid << end;
      break;
      
    case HEAT_EQUATION:
      ConvHist_file[0] << begin << heat_coeff;
      ConvHist_file[0] << heat_resid << end;
      break;
      
    case LINEAR_ELASTICITY:
      ConvHist_file[0] << begin << fea_coeff;
      ConvHist_file[0] << fea_resid << end;
      break;
      
  }
  
  if (config->GetOutput_FileFormat() == TECPLOT || config->GetOutput_FileFormat() == TECPLOT_BINARY) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }
  
}


void COutput::SetConvergence_History(ofstream *ConvHist_file,
                                     CGeometry ***geometry,
                                     CSolver ****solver_container,
                                     CConfig **config,
                                     CIntegration ***integration,
                                     bool DualTime_Iteration,
                                     double timeused,
                                     unsigned short val_iZone) {
  
  bool output_1d  = config[val_iZone]->GetWrt_1D_Output();
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  
  int rank;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = MASTER_NODE;
#endif
  
  /*--- If 1-D outputs requested, calculated them. Requires info from all nodes,
   Get area-averaged and flux-averaged values at the specified surface ---*/
  
  if (output_1d) {
    switch (config[val_iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case FLUID_STRUCTURE_EULER:   case FLUID_STRUCTURE_NAVIER_STOKES:   case FLUID_STRUCTURE_RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
        OneDimensionalOutput(solver_container[val_iZone][FinestMesh][FLOW_SOL], geometry[val_iZone][FinestMesh], config[val_iZone]);
    }
  }
  
  /*--- Output using only the master node ---*/
  if (rank == MASTER_NODE) {
    
    unsigned long iIntIter = config[val_iZone]->GetIntIter();
    unsigned long iExtIter = config[val_iZone]->GetExtIter();
    
    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/
    char begin[1000], direct_coeff[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000], resid_aux[1000], levelset_resid[1000], adj_levelset_resid[1000], wave_coeff[1000],
    heat_coeff[1000], fea_coeff[1000], wave_resid[1000], heat_resid[1000], fea_resid[1000], end[1000];
    char oneD_outputs[1000];
    double dummy = 0.0, *Coord;
    unsigned short iVar, iMarker, iMarker_Monitoring;
    
    unsigned long LinSolvIter = 0, iPointMaxResid;
    double timeiter = timeused/double(iExtIter+1);
    
    unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
    unsigned short nSpecies = config[val_iZone]->GetnSpecies();
    
    bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
    bool freesurface = (config[val_iZone]->GetKind_Regime() == FREESURFACE);
    
    bool rotating_frame = config[val_iZone]->GetRotating_Frame();
    bool aeroelastic = config[val_iZone]->GetAeroelastic_Simulation();
    bool equiv_area = config[val_iZone]->GetEquivArea();
    bool inv_design = (config[val_iZone]->GetInvDesign_Cp() || config[val_iZone]->GetInvDesign_HeatFlux());
    bool isothermal = false;
    for (iMarker = 0; iMarker < config[val_iZone]->GetnMarker_All(); iMarker++)
      if ((config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
          (config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC) ||
          (config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC))
        isothermal = true;
    bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                      (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool transition = (config[val_iZone]->GetKind_Trans_Model() == LM);
    bool adjoint = config[val_iZone]->GetAdjoint();
    bool fluid_structure = ((config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) ||
                            (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool wave = (config[val_iZone]->GetKind_Solver() == WAVE_EQUATION);
    bool heat = (config[val_iZone]->GetKind_Solver() == HEAT_EQUATION);
    bool fea = (config[val_iZone]->GetKind_Solver() == LINEAR_ELASTICITY);
    bool TNE2 = ((config[val_iZone]->GetKind_Solver() == TNE2_EULER) || (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) ||
                 (config[val_iZone]->GetKind_Solver() == ADJ_TNE2_EULER) || (config[val_iZone]->GetKind_Solver() == ADJ_TNE2_NAVIER_STOKES));
    bool flow = (config[val_iZone]->GetKind_Regime() == EULER) || (config[val_iZone]->GetKind_Regime() == NAVIER_STOKES) ||
    (config[val_iZone]->GetKind_Regime() == RANS) || (config[val_iZone]->GetKind_Regime() == ADJ_EULER) ||
    (config[val_iZone]->GetKind_Regime() == ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Regime() == ADJ_RANS);
    
    bool output_per_surface = false;
    if(config[val_iZone]->GetnMarker_Monitoring() > 1) output_per_surface = true;
    
    /*--- Initialize variables to store information from all domains (direct solution) ---*/
    double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
    Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0, Total_CHeat = 0.0, Total_CpDiff = 0.0, Total_HeatFluxDiff = 0.0,
    Total_CFEA = 0.0, Total_Heat = 0.0, Total_MaxHeat = 0.0;
    
    double OneD_AvgStagPress = 0.0, OneD_AvgMach = 0.0, OneD_AvgTemp = 0.0, OneD_MassFlowRate = 0.0,
    OneD_FluxAvgPress = 0.0, OneD_FluxAvgDensity = 0.0, OneD_FluxAvgVelocity = 0.0, OneD_FluxAvgEntalpy = 0.0;
    
    /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
    double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
    double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;
    
    /*--- Residual arrays ---*/
    double *residual_flow         = NULL,
    *residual_turbulent    = NULL,
    *residual_transition   = NULL,
    *residual_TNE2         = NULL,
    *residual_levelset     = NULL;
    double *residual_adjflow      = NULL,
    *residual_adjturbulent = NULL,
    *residual_adjTNE2      = NULL,
    *residual_adjlevelset  = NULL;
    double *residual_wave         = NULL;
    double *residual_fea          = NULL;
    double *residual_heat         = NULL;
    
    /*--- Coefficients Monitored arrays ---*/
    double *aeroelastic_plunge = NULL,
    *aeroelastic_pitch  = NULL,
    *Surface_CLift      = NULL,
    *Surface_CDrag      = NULL,
    *Surface_CSideForce = NULL,
    *Surface_CFx        = NULL,
    *Surface_CFy        = NULL,
    *Surface_CFz        = NULL,
    *Surface_CMx        = NULL,
    *Surface_CMy        = NULL,
    *Surface_CMz        = NULL;
    
    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0,
    nVar_Trans = 0, nVar_TNE2 = 0, nVar_Wave = 0, nVar_Heat = 0, nVar_FEA = 0,
    nVar_AdjFlow = 0, nVar_AdjTNE2 = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0;
    
    /*--- Direct problem variables ---*/
    if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+1;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()){
        case SA:	nVar_Turb = 1; break;
        case ML:	nVar_Turb = 1; break;
        case SST: nVar_Turb = 2; break;
      }
    }
    if (transition) nVar_Trans = 2;
    if (TNE2) nVar_TNE2 = config[val_iZone]->GetnSpecies()+nDim+2;
    if (wave) nVar_Wave = 2;
    if (fea) nVar_FEA = nDim;
    if (heat) nVar_Heat = 1;
    if (freesurface) nVar_LevelSet = 1;
    
    /*--- Adjoint problem variables ---*/
    if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+1;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()){
        case SA:	nVar_AdjTurb = 1; break;
        case ML:	nVar_AdjTurb = 1; break;
        case SST: nVar_AdjTurb = 2; break;
      }
    }
    if (TNE2) nVar_AdjTNE2 = config[val_iZone]->GetnSpecies()+nDim+2;
    if (freesurface) nVar_AdjLevelSet = 1;
    
    /*--- Allocate memory for the residual ---*/
    residual_flow       = new double[nVar_Flow];
    residual_turbulent  = new double[nVar_Turb];
    residual_transition = new double[nVar_Trans];
    residual_TNE2       = new double[nVar_TNE2];
    residual_levelset   = new double[nVar_LevelSet];
    residual_wave       = new double[nVar_Wave];
    residual_fea        = new double[nVar_FEA];
    residual_heat       = new double[nVar_Heat];
    
    residual_adjflow      = new double[nVar_AdjFlow];
    residual_adjturbulent = new double[nVar_AdjTurb];
    residual_adjTNE2      = new double[nVar_AdjTNE2];
    residual_adjlevelset  = new double[nVar_AdjLevelSet];
    
    /*--- Allocate memory for the coefficients being monitored ---*/
    aeroelastic_plunge = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    aeroelastic_pitch  = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CLift      = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CDrag      = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSideForce = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz        = new double[config[ZONE_0]->GetnMarker_Monitoring()];
    
    /*--- Write information from nodes ---*/
    switch (config[val_iZone]->GetKind_Solver()) {
        
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case FLUID_STRUCTURE_EULER:   case FLUID_STRUCTURE_NAVIER_STOKES:   case FLUID_STRUCTURE_RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
        
        /*--- Flow solution coefficients ---*/
        Total_CLift       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
        Total_CDrag       = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
        Total_CSideForce  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
        Total_CEff        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEff();
        Total_CMx         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
        Total_CMy         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
        Total_CMz         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
        Total_CFx         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
        Total_CFy         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
        Total_CFz         = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();
        
        if (freesurface) {
          Total_CFreeSurface = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFreeSurface();
        }
        
        if (isothermal) {
          Total_Heat     = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_HeatFlux();
          Total_MaxHeat  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_MaxHeatFlux();
        }
        
        if (equiv_area) {
          Total_CEquivArea    = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
          Total_CNearFieldOF  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
          
          /*--- Note that there is a redefinition of the nearfield based functionals ---*/
          Total_CEquivArea    = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CEquivArea;
          Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0-config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
        }
        
        if (inv_design) {
          Total_CpDiff  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
          if (isothermal) {
            Total_HeatFluxDiff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_HeatFluxDiff();
          }
        }
        
        if (rotating_frame) {
          Total_CT      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
          Total_CQ      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
          Total_CMerit  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
        }
        
        if (aeroelastic) {
          /*--- Look over the markers being monitored and get the desired values ---*/
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            aeroelastic_plunge[iMarker_Monitoring] = config[val_iZone]->GetAeroelastic_plunge(iMarker_Monitoring);
            aeroelastic_pitch[iMarker_Monitoring]  = config[val_iZone]->GetAeroelastic_pitch(iMarker_Monitoring);
          }
        }
        
        if (output_per_surface) {
          /*--- Look over the markers being monitored and get the desired values ---*/
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Surface_CLift[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
            Surface_CDrag[iMarker_Monitoring]      = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
            Surface_CSideForce[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce(iMarker_Monitoring);
            Surface_CFx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring);
            Surface_CFy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring);
            Surface_CFz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring);
            Surface_CMx[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
            Surface_CMy[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
            Surface_CMz[iMarker_Monitoring]        = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
          }
        }
        
        if (fluid_structure) {
          Total_CFEA  = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetTotal_CFEA();
        }
        
        if (output_1d) {
          
          /*--- Get area-averaged and flux-averaged values at the specified surface ---*/
          
          OneD_AvgStagPress = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_TotalPress();
          OneD_AvgMach = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_Mach();
          OneD_AvgTemp = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_Temp();
          OneD_MassFlowRate = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_MassFlowRate();
          
          OneD_FluxAvgPress = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgPress();
          OneD_FluxAvgDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgDensity();
          OneD_FluxAvgVelocity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgVelocity();
          OneD_FluxAvgEntalpy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgEntalpy();
          
        }
        
        /*--- Flow Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Flow; iVar++)
          residual_flow[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);
        
        /*--- Turbulent residual ---*/
        
        if (turbulent) {
          for (iVar = 0; iVar < nVar_Turb; iVar++)
            residual_turbulent[iVar] = solver_container[val_iZone][FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Transition residual ---*/
        
        if (transition) {
          for (iVar = 0; iVar < nVar_Trans; iVar++)
            residual_transition[iVar] = solver_container[val_iZone][FinestMesh][TRANS_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Free Surface residual ---*/
        
        if (freesurface) {
          for (iVar = 0; iVar < nVar_LevelSet; iVar++)
            residual_levelset[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(nDim+1);
        }
        
        /*--- FEA residual ---*/
        if (fluid_structure) {
          for (iVar = 0; iVar < nVar_FEA; iVar++)
            residual_fea[iVar] = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Iterations of the linear solver ---*/
        
        LinSolvIter = (unsigned long) solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetIterLinSolver();
        
        /*--- Adjoint solver ---*/
        
        if (adjoint) {
          
          /*--- Adjoint solution coefficients ---*/
          
          Total_Sens_Geo   = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
          Total_Sens_Mach  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
          Total_Sens_AoA   = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
          Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
          Total_Sens_Temp  = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();
          
          /*--- Adjoint flow residuals ---*/
          
          for (iVar = 0; iVar < nVar_AdjFlow; iVar++) {
            residual_adjflow[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(iVar);
          }
          
          /*--- Adjoint turbulent residuals ---*/
          
          if (turbulent) {
            if (!config[val_iZone]->GetFrozen_Visc()) {
              for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
                residual_adjturbulent[iVar] = solver_container[val_iZone][FinestMesh][ADJTURB_SOL]->GetRes_RMS(iVar);
            }
          }
          
          /*--- Adjoint level set residuals ---*/
          
          if (freesurface) {
            for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++)
              residual_adjlevelset[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(nDim+1);
          }
          
        }
        
        break;
        
      case TNE2_EULER:     case TNE2_NAVIER_STOKES:
      case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
        
        /*--- Coefficients ---*/
        
        Total_CLift       = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CLift();
        Total_CDrag       = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CDrag();
        Total_CSideForce  = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CSideForce();
        Total_CEff        = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CEff();
        Total_CMx         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMx();
        Total_CMy         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMy();
        Total_CMz         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMz();
        Total_CFx         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFx();
        Total_CFy         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFy();
        Total_CFz         = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFz();
        
        if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) {
          Total_Heat           = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_HeatFlux();
          Total_MaxHeat        = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_MaxHeatFlux();
          if (inv_design) {
            Total_HeatFluxDiff = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_HeatFluxDiff();
          }
        }
        
        /*--- Residuals ---*/
        
        for (iVar = 0; iVar < nVar_TNE2; iVar++)
          residual_TNE2[iVar] = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetRes_RMS(iVar);
        
        /*--- Iterations of the linear solver ---*/
        LinSolvIter = (unsigned long) solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetIterLinSolver();
        
        /*--- Adjoint solver ---*/
        if (adjoint) {
          
          /*--- Adjoint solution coefficients ---*/
          
          Total_Sens_Geo   = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Geo();
          Total_Sens_Mach  = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Mach();
          Total_Sens_AoA   = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_AoA();
          Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Press();
          Total_Sens_Temp  = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Temp();
          
          /*--- Adjoint flow residuals ---*/
          for (iVar = 0; iVar < nVar_AdjTNE2; iVar++) {
            residual_adjTNE2[iVar] = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetRes_RMS(iVar);
          }
        }
        
        break;
        
      case WAVE_EQUATION:
        
        /*--- Wave coefficients  ---*/
        
        Total_CWave = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();
        
        /*--- Wave Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Wave; iVar++) {
          residual_wave[iVar] = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetRes_RMS(iVar);
        }
        
        break;
        
      case HEAT_EQUATION:
        
        /*--- Heat coefficients  ---*/
        
        Total_CHeat = solver_container[val_iZone][FinestMesh][HEAT_SOL]->GetTotal_CHeat();
        
        /*--- Wave Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Heat; iVar++) {
          residual_heat[iVar] = solver_container[val_iZone][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
        }
        
        break;
        
      case LINEAR_ELASTICITY:
        
        /*--- FEA coefficients ---*/
        
        Total_CFEA = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();
        
        /*--- Plasma Residuals ---*/
        
        for (iVar = 0; iVar < nVar_FEA; iVar++) {
          residual_fea[iVar] = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
        }
        
        break;
        
    }
    
    /*--- Header frecuency ---*/
    
    bool Unsteady = ((config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                     (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    
    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*20)) == 0));
    
    if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {
      
      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {
        
        /*--- Write the begining of the history file ---*/
        sprintf (begin, "%12d", int(iExtIter));
        
        /*--- Write the end of the history file ---*/
        sprintf (end, ", %12.10f, %12.10f\n", double(LinSolvIter), timeused/60.0);
        
        /*--- Write the solution and residual of the history file ---*/
        switch (config[val_iZone]->GetKind_Solver()) {
            
          case EULER : case NAVIER_STOKES: case RANS:
          case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS:
          case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
            
            /*--- Direct coefficients ---*/
            sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff);
            if (isothermal)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat);
            if (equiv_area)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
            if (inv_design) {
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CpDiff);
              Total_CpDiff  = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
              if (isothermal) {
                sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat, Total_CpDiff, Total_HeatFluxDiff);
              }
            }
            if (rotating_frame)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);
            
            if (freesurface) {
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                       Total_CFz, Total_CEff, Total_CFreeSurface);
            }
            if (fluid_structure)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
                       Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);
            
            if (aeroelastic) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to aeroelastic coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  sprintf(aeroelastic_coeff, ", %12.10f",aeroelastic_plunge[iMarker_Monitoring]);
                }
                else {
                  sprintf(surface_coeff, ", %12.10f",aeroelastic_plunge[iMarker_Monitoring]);
                  strcat(aeroelastic_coeff, surface_coeff);
                }
                sprintf(surface_coeff, ", %12.10f",aeroelastic_pitch[iMarker_Monitoring]);
                strcat(aeroelastic_coeff, surface_coeff);
              }
            }
            
            if (output_per_surface) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to monitoring coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  sprintf(monitoring_coeff, ", %12.10f",Surface_CLift[iMarker_Monitoring]);
                }
                else {
                  sprintf(surface_coeff, ", %12.10f",Surface_CLift[iMarker_Monitoring]);
                  strcat(monitoring_coeff, surface_coeff);
                }
                sprintf(surface_coeff, ", %12.10f",Surface_CDrag[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CSideForce[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CFx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CFy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CFz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
              }
            }
            
            
            /*--- Flow residual ---*/
            if (nDim == 2) {
              if (compressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
              if (incompressible || freesurface) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), dummy, dummy );
            }
            else {
              if (compressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]) );
              if (incompressible || freesurface) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
            }
            
            /*--- Turbulent residual ---*/
            if (turbulent){
              switch(nVar_Turb) {
                case 1: sprintf (turb_resid, ", %12.10f", log10 (residual_turbulent[0])); break;
                case 2: sprintf (turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); break;
              }
            }
            /*---- Averaged stagnation pressure at an exit ---- */
            if (output_1d) {
              sprintf( oneD_outputs, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", OneD_AvgStagPress, OneD_AvgMach, OneD_AvgTemp, OneD_MassFlowRate, OneD_FluxAvgPress, OneD_FluxAvgDensity, OneD_FluxAvgVelocity, OneD_FluxAvgEntalpy);
            }
            
            /*--- Transition residual ---*/
            if (transition){
              sprintf (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }
            
            /*--- Free surface residual ---*/
            if (freesurface) {
              sprintf (levelset_resid, ", %12.10f", log10 (residual_levelset[0]));
            }
            
            /*--- Fluid structure residual ---*/
            if (fluid_structure) {
              if (nDim == 2) sprintf (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (residual_fea[0]), log10 (residual_fea[1]));
              else sprintf (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), log10 (residual_fea[1]), log10 (residual_fea[2]));
            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
              
              /*--- Adjoint flow residuals ---*/
              if (nDim == 2) {
                if (compressible) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
                if (incompressible || freesurface) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, 0.0, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]) );
              }
              else {
                if (compressible) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]), log10 (residual_adjflow[4]) );
                if (incompressible || freesurface) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
              }
              
              /*--- Adjoint turbulent residuals ---*/
              if (turbulent)
                if (!config[val_iZone]->GetFrozen_Visc())
                  sprintf (adj_turb_resid, ", %12.10f", log10 (residual_adjturbulent[0]));
              
              /*--- Adjoint free surface residuals ---*/
              if (freesurface) sprintf (adj_levelset_resid, ", %12.10f", log10 (residual_adjlevelset[0]));
            }
            
            break;
            
          case TNE2_EULER :    case TNE2_NAVIER_STOKES:
          case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
            
            /*--- Direct coefficients ---*/
            if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) {
              if (!(inv_design))
                sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                         Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                         Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                         Total_CEff, Total_Heat, Total_MaxHeat);
              else
                sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                         Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                         Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                         Total_CEff, Total_Heat, Total_MaxHeat, Total_HeatFluxDiff);
            }
            else
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                       Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                       Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                       Total_CEff);
            
            /*--- Direct problem residual ---*/
            for (iVar = 0; iVar < nSpecies+nDim+2; iVar++) {
              sprintf (resid_aux, ", %12.10f", log10 (residual_TNE2[iVar]));
              if (iVar == 0) strcpy(flow_resid, resid_aux);
              else strcat(flow_resid, resid_aux);
            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
              
              /*--- Adjoint flow residuals ---*/
              for (iVar = 0; iVar < nSpecies+nDim+2; iVar++) {
                sprintf (resid_aux, ", %12.10f", log10 (residual_adjTNE2[iVar]));
                if (iVar == 0) strcpy(adj_flow_resid, resid_aux);
                else strcat(adj_flow_resid, resid_aux);
              }
            }
            
            break;
            
          case WAVE_EQUATION:
            
            sprintf (direct_coeff, ", %12.10f", Total_CWave);
            sprintf (wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_wave[0]), log10 (residual_wave[1]), dummy, dummy, dummy );
            
            break;
            
          case HEAT_EQUATION:
            
            sprintf (direct_coeff, ", %12.10f", Total_CHeat);
            sprintf (heat_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_heat[0]), dummy, dummy, dummy, dummy );
            
            break;
            
          case LINEAR_ELASTICITY:
            
            sprintf (direct_coeff, ", %12.10f", Total_CFEA);
            sprintf (fea_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), dummy, dummy, dummy, dummy );
            
            break;
            
        }
      }
      
      /*--- Write the screen header---*/
      if ((write_heads) && !(!DualTime_Iteration && Unsteady)) {
        
        if (!Unsteady) {
          switch (config[val_iZone]->GetKind_Solver()) {
            case EULER :                  case NAVIER_STOKES:
            case FLUID_STRUCTURE_EULER :  case FLUID_STRUCTURE_NAVIER_STOKES:
              cout << endl << "Min Delta Time: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
              ". Max Delta Time: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".";
              break;
              
            case TNE2_EULER: case TNE2_NAVIER_STOKES:
            case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
              cout << endl << "Min Delta Time: " << solver_container[val_iZone][MESH_0][TNE2_SOL]->GetMin_Delta_Time()<< ". Max Delta Time: " << solver_container[val_iZone][MESH_0][TNE2_SOL]->GetMax_Delta_Time() << ".";
              break;
          }
        }
        else {
          if (flow) {
            cout << endl << "Min DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
            ".Max DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
            ".Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
          }
          else {
            cout << endl << "Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
          }
        }
        
        switch (config[val_iZone]->GetKind_Solver()) {
          case EULER :                  case NAVIER_STOKES:
          case FLUID_STRUCTURE_EULER :  case FLUID_STRUCTURE_NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;
            cout <<"Max residual point " << iPointMaxResid << " is located at (" << Coord[0] << ", " << Coord[1];
            if (nDim == 3) cout << ", " << Coord[2];
            cout <<   ")." << endl;
            
            /*--- Print out the number of non-physical points and reconstructions ---*/
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << " There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
              cout << " There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            if (!fluid_structure) {
              if (incompressible) cout << "   Res[Press]" << "     Res[Velx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
              else if (freesurface) cout << "   Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "     CLevelSet" << endl;
              else if (rotating_frame && nDim == 3) cout << "     Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
              else if (aeroelastic) cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
              else if (equiv_area) cout << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
              else cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            }
            else if (fluid_structure) cout << "     Res[Rho]" << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            
            break;
            
          case TNE2_EULER :  case TNE2_NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetRes_Max(0)) << "." << endl;
            cout <<"Max residual point " << iPointMaxResid << " is located at (" << Coord[0] << ", " << Coord[1];
            if (nDim == 3) cout << ", " << Coord[2];
            cout <<   ")." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            cout << "     Res[Rho]" << "     Res[RhoE]" << "   Res[RhoEve]" << "   CDrag(Total)";
            if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES)
              cout << "   HeatLoad(Total)" << endl;
            else cout << endl;
            break;
            
          case RANS : case FLUID_STRUCTURE_RANS:
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;
            cout <<"Max residual point " << iPointMaxResid << " is located at (" << Coord[0] << ", " << Coord[1];
            if (nDim == 3) cout << ", " << Coord[2];
            cout <<   ")." << endl;
            
            /*--- Print out the number of non-physical points and reconstructions ---*/
            if (config[val_iZone]->GetNonphysical_Points() > 0)
              cout << " There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
              cout << " There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            if (incompressible || freesurface) cout << "   Res[Press]";
            else cout << "      Res[Rho]";
            
            switch (config[val_iZone]->GetKind_Turb_Model()){
              case SA:	cout << "       Res[nu]"; break;
              case ML:	cout << "       Res[nu]"; break;
              case SST:	cout << "     Res[kine]" << "     Res[omega]"; break;
            }
            
            switch (config[val_iZone]->GetKind_Trans_Model()){
              case LM:	cout << "      Res[Int]" << "       Res[Re]"; break;
            }
            
            if (rotating_frame && nDim == 3 ) cout << "   CThrust(Total)" << "   CTorque(Total)" << endl;
            else if (aeroelastic) cout << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
            else if (equiv_area) cout << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
            else cout << "   CLift(Total)"   << "   CDrag(Total)"   << endl;
            
            break;
            
          case WAVE_EQUATION :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            cout << "      Res[Wave]" << "   CWave(Total)"<<  endl;
            break;
            
          case HEAT_EQUATION :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            cout << "      Res[Heat]" << "   CHeat(Total)"<<  endl;
            break;
            
          case LINEAR_ELASTICITY :
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (nDim == 2) cout << "    Res[Displx]" << "    Res[Disply]" << "   CFEA(Total)"<<  endl;
            if (nDim == 3) cout << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "   CFEA(Total)"<<  endl;
            break;
            
          case ADJ_EULER :              case ADJ_NAVIER_STOKES :
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
            cout <<"Max residual point " << iPointMaxResid << " is located at (" << Coord[0] << ", " << Coord[1];
            if (nDim == 3) cout << ", " << Coord[2];
            cout <<   ")." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible || freesurface) cout << "   Res[Psi_Press]" << "   Res[Psi_Velx]";
            else cout << "   Res[Psi_Rho]" << "     Res[Psi_E]";
            cout << "      Sens_Geo" << "     Sens_Mach" << endl;
            
            if (freesurface) {
              cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
            }
            break;
            
          case ADJ_RANS :
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
            cout <<"Max residual point " << iPointMaxResid << " is located at (" << Coord[0] << ", " << Coord[1];
            if (nDim == 3) cout << ", " << Coord[2];
            cout <<   ")." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            if (incompressible || freesurface) cout << "     Res[Psi_Press]";
            else cout << "     Res[Psi_Rho]";
            
            if (!config[val_iZone]->GetFrozen_Visc()) {
              cout << "      Res[Psi_nu]";
            }
            else {
              if (incompressible || freesurface) cout << "   Res[Psi_Velx]";
              else cout << "     Res[Psi_E]";
            }
            cout << "     Sens_Geo" << "    Sens_Mach" << endl;
            
            if (freesurface) {
              cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
            }
            break;
            
          case ADJ_TNE2_EULER :              case ADJ_TNE2_NAVIER_STOKES :
            
            /*--- Visualize the maximum residual ---*/
            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetPoint_Max(0);
            Coord = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetPoint_Max_Coord(0);
            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetRes_Max(0)) << "." << endl;
            cout <<"Max residual point " << iPointMaxResid << " is located at (" << Coord[0] << ", " << Coord[1];
            if (nDim == 3) cout << ", " << Coord[2];
            cout <<   ")." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << "  ExtIter";
            
            cout << "   Res[Psi_Rho]" << "     Res[Psi_E]" << "   Res[Psi_Eve]" << "     Sens_Geo" << "    Sens_Mach" << endl;
            
            break;
            
        }
        
      }
      
      /*--- Write the solution on the screen and history file ---*/
      cout.precision(6);
      cout.setf(ios::fixed,ios::floatfield);
      
      if (!Unsteady) {
        cout.width(5); cout << iExtIter;
        cout.width(11); cout << timeiter;
        
      } else {
        cout.width(8); cout << iIntIter;
        cout.width(8); cout << iExtIter;
      }
      
      switch (config[val_iZone]->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES:
        case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            if (compressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (incompressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
            if (freesurface) ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
            if (fluid_structure) ConvHist_file[0] << fea_resid;
            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            if (output_1d) ConvHist_file[0] << oneD_outputs;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
          if (!fluid_structure && !equiv_area) {
            if (compressible) {
              if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
              else { cout.width(14); cout << log10(residual_flow[4]); }
            }
            if (incompressible) { cout.width(14); cout << log10(residual_flow[1]); }
            if (freesurface) { cout.width(14); cout << log10(residual_levelset[0]); }
          }
          else if (fluid_structure) { cout.width(14); cout << log10(residual_fea[0]); }
          
          if (rotating_frame && nDim == 3 ) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << Total_CT;
            cout.width(15); cout << Total_CQ;
            cout.unsetf(ios_base::floatfield);
          }
          else if (equiv_area) { cout.width(15); cout << min(10000.0,max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0,max(-10000.0, Total_CDrag)); cout.width(15);
            cout.precision(4);
            cout.setf(ios::scientific,ios::floatfield);
            cout << Total_CNearFieldOF; }
          else if (freesurface) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CFreeSurface; }
          else { cout.width(15); cout << min(10000.0,max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0,max(-10000.0, Total_CDrag)); }
          if (aeroelastic) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
            cout.width(15); cout << aeroelastic_pitch[0];
            cout.unsetf(ios_base::floatfield);
          }
          cout << endl;
          
          break;
          
        case RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid;
            if (transition) ConvHist_file[0] << trans_resid;
            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            if (output_1d) ConvHist_file[0] << oneD_outputs;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          
          if (incompressible || freesurface) cout.width(13);
          else  cout.width(14);
          cout << log10(residual_flow[0]);
          
          switch(nVar_Turb) {
            case 1: cout.width(14); cout << log10(residual_turbulent[0]); break;
            case 2: cout.width(14); cout << log10(residual_turbulent[0]);
              cout.width(15); cout << log10(residual_turbulent[1]); break;
          }
          
          if (transition) { cout.width(14); cout << log10(residual_transition[0]); cout.width(14); cout << log10(residual_transition[1]); }
          
          if (rotating_frame && nDim == 3 ) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << Total_CT; cout.width(15);
            cout << Total_CQ;
            cout.unsetf(ios_base::floatfield);
          }
          else if (equiv_area) { cout.width(15); cout << min(10000.0,max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0,max(-10000.0, Total_CDrag)); cout.width(15);
            cout.precision(4);
            cout.setf(ios::scientific,ios::floatfield);
            cout << Total_CNearFieldOF; }
          else { cout.width(15); cout << min(10000.0,max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0,max(-10000.0, Total_CDrag)); }
          
          if (aeroelastic) {
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
            cout.width(15); cout << aeroelastic_pitch[0];
            cout.unsetf(ios_base::floatfield);
          }
          cout << endl;
          
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed,ios::floatfield);
            cout.width(13); cout << log10(residual_flow[0]);
            cout.width(14); cout << log10(residual_levelset[0]);
            cout.width(15); cout << Total_CLift;
            cout.width(14); cout << Total_CFreeSurface;
            
            cout << endl;
          }
          
          break;
          
        case TNE2_EULER : case TNE2_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_TNE2[0]);
          cout.width(14); cout << log10(residual_TNE2[nSpecies+nDim]);
          cout.width(14); cout << log10(residual_TNE2[nSpecies+nDim+1]);
          cout.width(15); cout << Total_CDrag;
          if (config[val_iZone]->GetKind_Solver()==TNE2_NAVIER_STOKES) {
            cout.precision(1);
            cout.width(11); cout << Total_MaxHeat;
          }
          cout << endl;
          break;
          
        case WAVE_EQUATION:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << wave_coeff << wave_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(14); cout << log10(residual_wave[0]);
          cout.width(14); cout << Total_CWave;
          cout << endl;
          break;
          
        case HEAT_EQUATION:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << heat_coeff << heat_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(14); cout << log10(residual_heat[0]);
          cout.width(14); cout << Total_CHeat;
          cout << endl;
          break;
          
        case LINEAR_ELASTICITY:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << fea_coeff << fea_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(15); cout << log10(residual_fea[0]);
          cout.width(15); cout << log10(residual_fea[1]);
          if (nDim == 3) { cout.width(15); cout << log10(residual_fea[2]); }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_CFEA;
          cout << endl;
          break;
          
        case ADJ_EULER :              case ADJ_NAVIER_STOKES :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          if (compressible) {
            cout.width(15); cout << log10(residual_adjflow[0]);
            cout.width(15); cout << log10(residual_adjflow[nDim+1]);
          }
          if (incompressible || freesurface) {
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjflow[1]);
          }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed,ios::floatfield);
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjlevelset[0]);
            cout.precision(3);
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(12); cout << Total_Sens_Geo;
            cout.width(12); cout << Total_Sens_Mach;
            cout.unsetf(ios_base::floatfield);
            cout << endl;
          }
          
          break;
          
        case ADJ_RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid;
            if (!config[val_iZone]->GetFrozen_Visc())
              ConvHist_file[0] << adj_turb_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(17); cout << log10(residual_adjflow[0]);
          if (!config[val_iZone]->GetFrozen_Visc()) {
            cout.width(17); cout << log10(residual_adjturbulent[0]);
          }
          else {
            if (compressible) {
              if (geometry[val_iZone][FinestMesh]->GetnDim() == 2 ) { cout.width(15); cout << log10(residual_adjflow[3]); }
              else { cout.width(15); cout << log10(residual_adjflow[4]); }
            }
            if (incompressible || freesurface) {
              cout.width(15); cout << log10(residual_adjflow[1]);
            }
          }
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          if (freesurface) {
            if (!DualTime_Iteration) {
              ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid;
              ConvHist_file[0] << end;
              ConvHist_file[0].flush();
            }
            
            cout.precision(6);
            cout.setf(ios::fixed,ios::floatfield);
            cout.width(17); cout << log10(residual_adjflow[0]);
            cout.width(16); cout << log10(residual_adjlevelset[0]);
            
            cout.precision(4);
            cout.setf(ios::scientific,ios::floatfield);
            cout.width(12); cout << Total_Sens_Geo;
            cout.width(12); cout << Total_Sens_Mach;
            cout << endl;
            cout.unsetf(ios_base::floatfield);
          }
          
          break;
          
        case ADJ_TNE2_EULER :              case ADJ_TNE2_NAVIER_STOKES :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(15); cout << log10(residual_adjTNE2[0]);
          cout.width(15); cout << log10(residual_adjTNE2[nSpecies+nDim]);
          cout.width(15); cout << log10(residual_adjTNE2[nSpecies+nDim+1]);
          
          cout.precision(4);
          cout.setf(ios::scientific,ios::floatfield);
          cout.width(14); cout << Total_Sens_Geo;
          cout.width(14); cout << Total_Sens_Mach;
          cout << endl;
          cout.unsetf(ios_base::floatfield);
          
          break;
          
          
      }
      cout.unsetf(ios::fixed);
      
      delete [] residual_flow;
      delete [] residual_turbulent;
      delete [] residual_transition;
      delete [] residual_TNE2;
      delete [] residual_levelset;
      delete [] residual_wave;
      delete [] residual_fea;
      delete [] residual_heat;
      
      delete [] residual_adjflow;
      delete [] residual_adjturbulent;
      delete [] residual_adjTNE2;
      delete [] residual_adjlevelset;
      
      delete [] Surface_CLift;
      delete [] Surface_CDrag;
      delete [] Surface_CSideForce;
      delete [] Surface_CFx;
      delete [] Surface_CFy;
      delete [] Surface_CFz;
      delete [] Surface_CMx;
      delete [] Surface_CMy;
      delete [] Surface_CMz;
      delete [] aeroelastic_pitch;
      delete [] aeroelastic_plunge;
    }
  }
}

void COutput::SetResult_Files(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
                              unsigned long iExtIter, unsigned short val_nZone) {
  
  int rank = MASTER_NODE;
  
#ifdef HAVE_MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    
#ifdef HAVE_MPI
    /*--- Do not merge the volume solutions if we are running in parallel.
     Force the use of SU2_SOL to merge the volume sols in this case. ---*/
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > SINGLE_NODE) {
      Wrt_Vol = false;
      Wrt_Srf = false;
    }
#endif
    
    bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();
    bool Wrt_Rst = config[iZone]->GetWrt_Restart();
    
    switch (config[iZone]->GetKind_Solver()) {
        
      case EULER : case NAVIER_STOKES : case RANS :
      case FLUID_STRUCTURE_EULER : case FLUID_STRUCTURE_NAVIER_STOKES : case FLUID_STRUCTURE_RANS:
        
        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
        break;
        
      case TNE2_EULER : case TNE2_NAVIER_STOKES :
        
        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][TNE2_SOL], iExtIter, iZone);
        break;
        
      case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][ADJFLOW_SOL], solver_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
        break;
        
      case ADJ_TNE2_EULER : case ADJ_TNE2_NAVIER_STOKES :
        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][ADJTNE2_SOL], solver_container[iZone][MESH_0][TNE2_SOL], iExtIter, iZone);
        break;
        
      case LIN_EULER : case LIN_NAVIER_STOKES :
        if (Wrt_Csv) SetSurfaceCSV_Linearized(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0][LINFLOW_SOL], config[iZone]->GetSurfLinCoeff_FileName(), iExtIter);
        break;
    }
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    bool dynamic_mesh = (config[iZone]->GetUnsteady_Simulation() &&
                         config[iZone]->GetGrid_Movement());
    
    /*--- Merge the node coordinates and connectivity, if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf)
      MergeConnectivity(config[iZone], geometry[iZone][MESH_0], iZone);
    
    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/
    
    MergeCoordinates(config[iZone], geometry[iZone][MESH_0]);
    
    if (rank == MASTER_NODE) {
      
      if (FileFormat == CGNS_SOL) {
        SetCGNS_Coordinates(config[iZone], geometry[iZone][MESH_0], iZone);
        if (!wrote_base_file || dynamic_mesh)
          DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
      } else if (FileFormat == TECPLOT_BINARY) {
        SetTecplot_Mesh(config[iZone], geometry[iZone][MESH_0], iZone);
        SetTecplot_SurfaceMesh(config[iZone], geometry[iZone][MESH_0], iZone);
        if (!wrote_base_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
        if (!wrote_surf_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], wrote_surf_file);
      }
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (Wrt_Vol || Wrt_Rst)
      MergeSolution(config[iZone], geometry[iZone][MESH_0],
                    solver_container[iZone][MESH_0], iZone);
    
    /*--- Write restart, CGNS, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- Write a native restart file ---*/
      if (Wrt_Rst)
        SetRestart(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] ,iZone);
      
      if (Wrt_Vol) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0],iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_Solution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case CGNS_SOL:
            
            /*--- Write a CGNS solution file ---*/
            SetCGNS_Solution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
            break;
            
          default:
            break;
        }
        
      }
      
      if (Wrt_Srf) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone][MESH_0],solver_container[iZone][MESH_0] ,iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_SurfaceSolution(config[iZone], geometry[iZone][MESH_0], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
            break;
            
          default:
            break;
        }
        
      }
      
      /*--- Release memory needed for merging the solution data. ---*/
      if (((Wrt_Vol) || (Wrt_Srf)) && (FileFormat == TECPLOT ||
                                       FileFormat == TECPLOT_BINARY ||
                                       FileFormat == PARAVIEW))
        DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
      
      if (Wrt_Vol || Wrt_Rst)
        DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);
      
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written) & barrier to sync up after master node writes
     output files. ---*/
    
#ifdef HAVE_MPI
    MPI_Bcast(&wrote_base_file, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  }
}

void COutput::SetBaselineResult_Files(CSolver **solver, CGeometry **geometry, CConfig **config,
                                      unsigned long iExtIter, unsigned short val_nZone) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short iZone;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    bool Wrt_Rst = config[iZone]->GetWrt_Restart();
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf) {
      if (rank == MASTER_NODE) cout <<"Merging grid connectivity." << endl;
      MergeConnectivity(config[iZone], geometry[iZone], iZone);
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (Wrt_Vol || Wrt_Rst) {
      if (rank == MASTER_NODE) cout <<"Merging solution." << endl;
      MergeBaselineSolution(config[iZone], geometry[iZone], solver[iZone], iZone);
    }
    
    /*--- Write restart, CGNS, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      if (Wrt_Vol) {
        
        if (rank == MASTER_NODE)
          cout <<"Writing volume solution file." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone], solver,iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone], false);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_Mesh(config[iZone], geometry[iZone], iZone);
            SetTecplot_Solution(config[iZone], geometry[iZone], iZone);
            break;
            
          case CGNS_SOL:
            
            /*--- Write a CGNS solution file ---*/
            SetCGNS_Solution(config[iZone], geometry[iZone], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone], false);
            break;
            
          default:
            break;
        }
        
      }
      
      if (Wrt_Srf) {
        
        if (rank == MASTER_NODE) cout <<"Writing surface solution file." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config[iZone], geometry[iZone],solver, iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone], true);
            break;
            
          case TECPLOT_BINARY:
            
            /*--- Write a Tecplot binary solution file ---*/
            SetTecplot_SurfaceMesh(config[iZone], geometry[iZone], iZone);
            SetTecplot_SurfaceSolution(config[iZone], geometry[iZone], iZone);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone], true);
            break;
            
          default:
            break;
        }
      }
      
      if (FileFormat == TECPLOT_BINARY) {
        if (!wrote_base_file)
          DeallocateConnectivity(config[iZone], geometry[iZone], false);
        if (!wrote_surf_file)
          DeallocateConnectivity(config[iZone], geometry[iZone], wrote_surf_file);
      }
      
      if (Wrt_Vol || Wrt_Srf)
        DeallocateSolution(config[iZone], geometry[iZone]);
    }
    
    /*--- Final broadcast (informing other procs that the base output
     file was written) & barrier to sync up after master node writes
     output files. ---*/
    
#ifdef HAVE_MPI
    MPI_Bcast(&wrote_base_file, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  }
}

void COutput::OneDimensionalOutput(CSolver *solver_container, CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Out1D;
  double *Normal = NULL, Area = 0.0, OverArea = 0.0, *Coord = NULL, UnitaryNormal[3],
  Stag_Pressure, Mach, Temperature, Pressure = 0.0, Density = 0.0, Velocity2, Enthalpy, RhoU, U,// local values at each node (Velocity2 = V^2). U = normal velocity
  SumPressure = 0.0, SumStagPressure = 0.0, SumArea = 0.0, SumMach = 0.0, SumTemperature = 0.0, SumForUref = 0.0, SumRhoU = 0.0, SumEnthalpy = 0.0,// sum of (local value ) * (dA) (integral)
  AveragePressure = 0.0, AverageMach = 0.0, AverageTemperature = 0.0, MassFlowRate = 0.0, // Area Averaged value ( sum / A )
  VelocityRef = 0.0, EnthalpyRef = 0.0, DensityRef = 0.0, PressureRef = 0.0; // Flux conserved values. TemperatureRef follows ideal gas
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  double Gamma = config->GetGamma();
  unsigned short nDim = geometry->GetnDim();
  
  
  /*--- Loop over the markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    Out1D = config->GetMarker_All_Out_1D(iMarker);
    
    /*--- Loop over the vertices to compute the output ---*/
    
    
    if (Out1D == YES) {
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Find the normal direction ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          
          /*--- Compute area, and unitary normal ---*/
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = -Normal[iDim]/Area;
          
          Coord = geometry->node[iPoint]->GetCoord();
          
          if (compressible){
            Pressure = solver_container->node[iPoint]->GetPressure();
            Density = solver_container->node[iPoint]->GetDensity();
          }
          if (incompressible || freesurface){
            Pressure = solver_container->node[iPoint]->GetPressureInc();
            Density = solver_container->node[iPoint]->GetDensityInc();
          }
          
          /*-- Find velocity normal to the marked surface/opening --*/
          
          U = 0.0; Enthalpy = 0.0; Mach = 0.0; Velocity2 = 0.0; Stag_Pressure = 0.0;
          Temperature = 0.0; RhoU = 0.0;
          for (iDim = 0; iDim < geometry->GetnDim(); iDim++){
            U += UnitaryNormal[iDim]*solver_container->node[iPoint]->GetVelocity(iDim);
          }
          
          Enthalpy = solver_container->node[iPoint]->GetEnthalpy();
          Velocity2 = solver_container->node[iPoint]->GetVelocity2();
          Temperature = solver_container->node[iPoint]->GetTemperature();
          
          Mach = (sqrt(Velocity2))/ solver_container->node[iPoint]->GetSoundSpeed();
          Stag_Pressure = Pressure*pow((1.0+((Gamma-1.0)/2.0)*pow(Mach, 2.0)),(Gamma/(Gamma-1.0)));
          
          RhoU = U*Density;
          SumStagPressure += Stag_Pressure * Area;
          SumArea += Area;
          SumMach += Mach*Area;
          SumPressure += Pressure * Area;
          SumTemperature += Temperature*Area;
          SumRhoU += RhoU*Area;
          SumForUref+=RhoU*U*U*Area;
          SumEnthalpy+=RhoU*Enthalpy*Area;
          
        }
      }
      
      if (SumRhoU != 0.0) { // To avoid division by 0
        
        OverArea = 1.0/SumArea;
        AveragePressure += abs(SumStagPressure*OverArea);
        AverageMach += abs(SumMach*OverArea);
        AverageTemperature += abs(SumTemperature*OverArea);
        MassFlowRate += SumRhoU;
        PressureRef += abs(SumPressure*OverArea);
        VelocityRef += abs(sqrt(abs(SumForUref/SumRhoU)));
        EnthalpyRef +=abs(SumEnthalpy/SumRhoU);
        DensityRef +=abs(PressureRef*Gamma/(Gamma-1)/(EnthalpyRef-0.5*VelocityRef*VelocityRef));
        
      }
      
    }
    
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  double My_AveragePressure     = AveragePressure;    AveragePressure = 0.0;
  double My_AverageMach         = AverageMach;        AverageMach = 0.0;
  double My_AverageTemperature  = AverageTemperature; AverageTemperature = 0.0;
  double My_MassFlowRate        = MassFlowRate;       MassFlowRate = 0.0;
  double My_PressureRef         = PressureRef;        PressureRef = 0.0;
  double My_VelocityRef         = VelocityRef;        VelocityRef = 0.0;
  double My_EnthalpyRef         = EnthalpyRef;        EnthalpyRef = 0.0;
  double My_DensityRef          = DensityRef;         DensityRef = 0.0;
  
  MPI_Allreduce(&My_AveragePressure, &AveragePressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_AverageMach, &AverageMach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_AverageTemperature, &AverageTemperature, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_MassFlowRate, &MassFlowRate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_PressureRef, &PressureRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_VelocityRef, &VelocityRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_EnthalpyRef , &EnthalpyRef , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&My_DensityRef , &DensityRef , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif
  
  /*--- Set the 1D output ---*/
  
  solver_container->SetOneD_TotalPress(AveragePressure);
  solver_container->SetOneD_Mach(AverageMach);
  solver_container->SetOneD_Temp(AverageTemperature);
  solver_container->SetOneD_MassFlowRate(MassFlowRate);
  
  solver_container->SetOneD_FluxAvgPress(PressureRef);
  solver_container->SetOneD_FluxAvgDensity(DensityRef);
  solver_container->SetOneD_FluxAvgVelocity(VelocityRef);
  solver_container->SetOneD_FluxAvgEntalpy(EnthalpyRef);
  
}

void COutput::SetForceSections(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  short iSection, nSection;
  unsigned long iVertex, iPoint;
  double *Plane_P0, *Plane_Normal, MinPlane, MaxPlane, *CPressure, MinXCoord, MaxXCoord, Force[3], ForceInviscid[3],
  MomentInviscid[3], MomentDist[3], RefDensity, RefPressure, RefAreaCoeff, Pressure_Inf, *Velocity_Inf, Gas_Constant, Mach2Vel, Mach_Motion, Gamma, RefVel2 = 0.0, factor, NDPressure, *Origin, RefLengthMoment, Alpha, Beta, CDrag_Inv, CLift_Inv, CMy_Inv;
  vector<double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, Pressure_Airfoil;
  string Marker_Tag, Slice_Filename, Slice_Ext;
  ofstream Cp_File;
  unsigned short iDim;
  
  bool grid_movement = config->GetGrid_Movement();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  
  Plane_P0 = new double [3];
  Plane_Normal = new double [3];
  CPressure = new double[geometry->GetnPoint()];
  
  /*--- Compute some reference quantities and necessary values ---*/
  RefDensity = solver_container->GetDensity_Inf();
  RefPressure = solver_container->GetPressure_Inf();
  RefAreaCoeff = config->GetRefAreaCoeff();
  Velocity_Inf = solver_container->GetVelocity_Inf();
  Pressure_Inf = solver_container->GetPressure_Inf();
  Gamma = config->GetGamma();
  Origin = config->GetRefOriginMoment(0);
  RefLengthMoment  = config->GetRefLengthMoment();
  Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  Beta             = config->GetAoS()*PI_NUMBER/180.0;
  
  if (grid_movement) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (geometry->GetnDim() == 3) {
    
    /*--- Copy the pressure to an auxiliar structure ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      if (compressible) {
        CPressure[iPoint] = (solver_container->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff;
      }
      if (incompressible || freesurface) {
        CPressure[iPoint] = (solver_container->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff;
      }
    }
    
    nSection = config->GetnSections();
    
    for (iSection = 0; iSection < nSection; iSection++) {
      
      /*--- Read the values from the config file ---*/
      
      MinPlane = config->GetSection_Location(0); MaxPlane = config->GetSection_Location(1);
      MinXCoord = -1E6; MaxXCoord = 1E6;
      
      Plane_Normal[0] = 0.0;    Plane_P0[0] = 0.0;
      Plane_Normal[1] = 0.0;    Plane_P0[1] = 0.0;
      Plane_Normal[2] = 0.0;    Plane_P0[2] = 0.0;
      
      Plane_Normal[config->GetAxis_Orientation()] = 1.0;
      Plane_P0[config->GetAxis_Orientation()] = MinPlane + iSection*(MaxPlane - MinPlane)/double(nSection-1);
      
      /*--- Compute the airfoil sections (note that we feed in the Cp) ---*/
      
      geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, iSection,
                                       MinXCoord, MaxXCoord, CPressure,
                                       Xcoord_Airfoil, Ycoord_Airfoil,
                                       Zcoord_Airfoil, Pressure_Airfoil, true,
                                       config);
      
      if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() == 0)) {
        cout << "Please check the config file, the section "<< iSection+1 <<" has not been detected." << endl;
      }
      
      /*--- Output the pressure on each section (tecplot format) ---*/
      
      if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() != 0)) {
        
        /*--- Write Cp at each section ---*/
        
        ofstream Cp_File;
        if (iSection == 0) {
          Cp_File.open("Cp_Sections.dat", ios::out);
          Cp_File << "TITLE = \"Airfoil sections\"" << endl;
          Cp_File << "VARIABLES = \"X\",\"Y\",\"Z\",\"Cp\"" << endl;
        }
        else Cp_File.open("Cp_Sections.dat", ios::app);
        
        Cp_File << "ZONE T=\"SECTION_"<< (iSection+1) << "\", NODES= "<< Xcoord_Airfoil.size() << ", ELEMENTS= " << Xcoord_Airfoil.size()-1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
        
        /*--- Coordinates and pressure value ---*/
        
        if (config->GetSystemMeasurements() == SI) {
          for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
            Cp_File << Xcoord_Airfoil[iVertex] <<" "<< Ycoord_Airfoil[iVertex] <<" "<< Zcoord_Airfoil[iVertex] <<" "<< Pressure_Airfoil[iVertex] <<  endl;
          }
        }
        if (config->GetSystemMeasurements() == US) {
          for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
            Cp_File << Xcoord_Airfoil[iVertex]*12.0 <<" "<< Ycoord_Airfoil[iVertex]*12.0 <<" "<< Zcoord_Airfoil[iVertex]*12.0 <<" "<< Pressure_Airfoil[iVertex] <<  endl;
          }
        }
        
        /*--- Basic conectivity ---*/
        
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
          Cp_File << iVertex << "\t" << iVertex+1 << "\n";
        }
        
        Cp_File.close();
        
        
        /*--- Compute load distribution ---*/
        
        ForceInviscid[0] = 0.0; ForceInviscid[1] = 0.0; ForceInviscid[2] = 0.0;
        
        for (iVertex = 0; iVertex < Xcoord_Airfoil.size()-1; iVertex++) {
          
          NDPressure = 0.5*(Pressure_Airfoil[iVertex]+Pressure_Airfoil[iVertex+1]);
          
          Force[0] = -(Zcoord_Airfoil[iVertex+1] - Zcoord_Airfoil[iVertex])*NDPressure;
          Force[1] = 0.0;
          Force[2] = (Xcoord_Airfoil[iVertex+1] - Xcoord_Airfoil[iVertex])*NDPressure;
          
          ForceInviscid[0] += Force[0];
          ForceInviscid[1] += Force[1];
          ForceInviscid[2] += Force[2];
          
          MomentDist[0] = 0.5*(Xcoord_Airfoil[iVertex] + Xcoord_Airfoil[iVertex+1]) - Origin[0];
          MomentDist[1] = 0.5*(Ycoord_Airfoil[iVertex] + Ycoord_Airfoil[iVertex+1]) - Origin[1];
          MomentDist[2] = 0.5*(Zcoord_Airfoil[iVertex] + Zcoord_Airfoil[iVertex+1]) - Origin[3];
          
          MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
          
        }
        
        CLift_Inv = fabs( -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha));
        CDrag_Inv = fabs( ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta));
        CMy_Inv = MomentInviscid[1];
        
        
        /*--- Write load distribution ---*/
        
        ofstream Load_File;
        if (iSection == 0) {
          Load_File.open("Load_Distribution.dat", ios::out);
          Load_File << "TITLE = \"Load distribution\"" << endl;
          Load_File << "VARIABLES = \"Y\",\"C<sub>L</sub>\",\"C<sub>D</sub>\",\"C<supb>My</sub>\"" << endl;
          Load_File << "ZONE T=\"Wing load distribution\", NODES= "<< nSection << ", ELEMENTS= " << nSection-1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
        }
        else Load_File.open("Load_Distribution.dat", ios::app);
        
        /*--- Coordinates and pressure value ---*/
        
        Load_File << Ycoord_Airfoil[0] <<" "<< CLift_Inv <<" "<< CDrag_Inv  <<" "<< CMy_Inv << endl;
        
        /*--- Basic conectivity ---*/
        
        if (iSection == nSection-1) {
          for (iSection = 1; iSection < nSection; iSection++) {
            Load_File << iSection << "\t" << iSection+1 << "\n";
          }
        }
        
        Load_File.close();
        
        
      }
      
    }
    
    
  }
  
  /*--- Delete dynamically allocated memory ---*/
  
  delete [] Plane_P0;
  delete [] Plane_Normal;
  delete [] CPressure;
  
}

void COutput::SetCp_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  unsigned short iMarker, icommas, Boundary, Monitoring, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = NULL, Area, PressDiff;
  bool *PointInDomain;
  string text_line, surfCp_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];
  
  
  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  MPI_Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nPointGlobal = nPointLocal;
#endif
  
  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];
  
  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
        (Boundary == ISOTHERMAL_NONCATALYTIC) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- The Pressure file uses the global numbering ---*/
        
#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
        
        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetCPressureTarget(iMarker, iVertex, 0.0);
        }
        
      }
    }
  }
  
  /*--- Prepare to read the surface pressure files (CSV) ---*/
  
  surfCp_filename = "TargetCp";
  strcpy (cstr, surfCp_filename.c_str());
  
  /*--- Write file name with extension if unsteady or steady ---*/
  
  if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
      (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.dat", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.dat",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.dat",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
  }
  else
    sprintf (buffer, ".dat");
  
  strcat (cstr, buffer);
  
  /*--- Read the surface pressure file ---*/
  
  string::size_type position;
  
  Surface_file.open(cstr, ios::in);
  
  if (!(Surface_file.fail())) {
    
    getline(Surface_file,text_line);
    
    while (getline(Surface_file,text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if(position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      
      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        solver_container->SetCPressureTarget(iMarker, iVertex, PressureCoeff);
        
      }
      
    }
    
    Surface_file.close();
    
  }
  
  /*--- Compute the pressure difference ---*/
  
  PressDiff = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
        (Boundary == ISOTHERMAL_NONCATALYTIC) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        Cp = solver_container->GetCPressure(iMarker, iVertex);
        CpTarget = solver_container->GetCPressureTarget(iMarker, iVertex);
        
        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);
        
        PressDiff += Area * (CpTarget - Cp) * (CpTarget - Cp);
      }
      
    }
  }
  
#ifdef HAVE_MPI
  double MyPressDiff = PressDiff;   PressDiff = 0.0;
  MPI_Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Update the total Cp difference coeffient ---*/
  
  solver_container->SetTotal_CpDiff(PressDiff);
  
  delete[] Point2Vertex;
  
}

void COutput::SetHeat_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  unsigned short iMarker, icommas, Boundary, Monitoring, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  double XCoord, YCoord, ZCoord, PressureCoeff, HeatFlux = 0.0, HeatFluxDiff, HeatFluxTarget, *Normal = NULL, Area,
  Pressure, Cf;
  bool *PointInDomain;
  string text_line, surfHeatFlux_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];
  
  
  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  MPI_Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nPointGlobal = nPointLocal;
#endif
  
  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];
  
  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
        (Boundary == ISOTHERMAL_NONCATALYTIC) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- The Pressure file uses the global numbering ---*/
        
#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
        
        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetHeatFluxTarget(iMarker, iVertex, 0.0);
        }
      }
    }
  }
  
  /*--- Prepare to read the surface pressure files (CSV) ---*/
  
  surfHeatFlux_filename = "TargetHeatFlux";
  strcpy (cstr, surfHeatFlux_filename.c_str());
  
  /*--- Write file name with extension if unsteady or steady ---*/
  
  if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
      (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.dat", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.dat",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.dat",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
  }
  else
    sprintf (buffer, ".dat");
  
  strcat (cstr, buffer);
  
  /*--- Read the surface pressure file ---*/
  
  string::size_type position;
  
  Surface_file.open(cstr, ios::in);
  
  if (!(Surface_file.fail())) {
    
    getline(Surface_file,text_line);
    
    while (getline(Surface_file,text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if(position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      
      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        solver_container->SetHeatFluxTarget(iMarker, iVertex, HeatFlux);
        
      }
      
    }
    
    Surface_file.close();
  }
  
  /*--- Compute the pressure difference ---*/
  
  HeatFluxDiff = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
        (Boundary == ISOTHERMAL_NONCATALYTIC) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        HeatFlux = solver_container->GetHeatFlux(iMarker, iVertex);
        HeatFluxTarget = solver_container->GetHeatFluxTarget(iMarker, iVertex);
        
        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);
        
        HeatFluxDiff += Area * (HeatFluxTarget - HeatFlux) * (HeatFluxTarget - HeatFlux);
        
      }
      
    }
  }
  
#ifdef HAVE_MPI
  double MyHeatFluxDiff = HeatFluxDiff;   HeatFluxDiff = 0.0;
  MPI_Allreduce(&MyHeatFluxDiff, &HeatFluxDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Update the total HeatFlux difference coeffient ---*/
  
  solver_container->SetTotal_HeatFluxDiff(HeatFluxDiff);
  
  delete[] Point2Vertex;
  
}

void COutput::SetEquivalentArea(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {
  
  ofstream EquivArea_file, FuncGrad_file;
  unsigned short iMarker = 0, iDim;
  short *AzimuthalAngle = NULL;
  double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign = 0.0, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion,
  *Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf, Density_Inf,
  RefAreaCoeff, ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
  *Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
  *Weight = NULL, jFunction, jp1Function;
  unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint,
  *IdPoint = NULL, *IdDomain = NULL, auxDomain;
  unsigned short iPhiAngle;
  ofstream NearFieldEA_file; ifstream TargetEA_file;
  
  double XCoordBegin_OF = config->GetEA_IntLimit(0);
  double XCoordEnd_OF = config->GetEA_IntLimit(1);
  
  unsigned short nDim = geometry->GetnDim();
  double AoA = -(config->GetAoA()*PI_NUMBER/180.0);
  double EAScaleFactor = config->GetEA_ScaleFactor(); // The EA Obj. Func. should be ~ force based Obj. Func.
  
  int rank = MESH_0;
  
  Mach  = config->GetMach();
  Gamma = config->GetGamma();
  Beta = sqrt(Mach*Mach-1.0);
  R_Plane = fabs(config->GetEA_IntLimit(2));
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Density_Inf = config->GetDensity_FreeStreamND();
  RefAreaCoeff = config->GetRefAreaCoeff();
  Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
  Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
  Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
  ModVelocity_Inf = 0;
  for (iDim = 0; iDim < 3; iDim++)
    ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];
  ModVelocity_Inf = sqrt(ModVelocity_Inf);
  
  factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);
  
#ifndef HAVE_MPI
  
  /*--- Compute the total number of points on the near-field ---*/
  
  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        /*--- Using Face_Normal(z), and Coord(z) we identify only a surface,
         note that there are 2 NEARFIELD_BOUNDARY surfaces ---*/
        
        if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) nVertex_NearField ++;
      }
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  
  Xcoord = new double[nVertex_NearField];
  Ycoord = new double[nVertex_NearField];
  Zcoord = new double[nVertex_NearField];
  AzimuthalAngle = new short[nVertex_NearField];
  IdPoint = new unsigned long[nVertex_NearField];
  IdDomain = new unsigned long[nVertex_NearField];
  Pressure = new double[nVertex_NearField];
  FaceArea = new double[nVertex_NearField];
  EquivArea = new double[nVertex_NearField];
  TargetArea = new double[nVertex_NearField];
  NearFieldWeight = new double[nVertex_NearField];
  Weight = new double[nVertex_NearField];
  
  /*--- Copy the boundary information to an array ---*/
  
  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
          
          IdPoint[nVertex_NearField] = iPoint;
          Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
          Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
          
          if (nDim ==2) {
            AzimuthalAngle[nVertex_NearField] = 0;
          }
          
          if (nDim == 3) {
            Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
            
            /*--- Rotate the nearfield cylinder (AoA) only 3D ---*/
            
            double YcoordRot = Ycoord[nVertex_NearField];
            double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);
            
            /*--- Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/
            
            double AngleDouble; short AngleInt;
            AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);
            
            /*--- Fix an azimuthal line due to misalignments of the near-field ---*/
            
            double FixAzimuthalLine = config->GetFixAzimuthalLine();
            
            if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1)) AngleDouble = FixAzimuthalLine - 0.1;
            
            AngleInt = (short) floor(AngleDouble + 0.5);
            if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
            else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
          }
          
          if (AzimuthalAngle[nVertex_NearField] <= 60) {
            Pressure[nVertex_NearField] = solver_container->node[iPoint]->GetPressure();
            FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nVertex_NearField ++;
          }
          
        }
      }
  
#else
  
  int nProcessor;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
  int iProcessor;
  
  unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  
  /*--- Compute the total number of points of the near-field ghost nodes ---*/
  
  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
            nLocalVertex_NearField ++;
      }
  
  Buffer_Send_nVertex[0] = nLocalVertex_NearField;
  
  /*--- Send Near-Field vertex information --*/
  
  MPI_Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
  double *Buffer_Send_Ycoord = new double[MaxLocalVertex_NearField];
  double *Buffer_Send_Zcoord = new double[MaxLocalVertex_NearField];
  unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
  double *Buffer_Send_Pressure = new double [MaxLocalVertex_NearField];
  double *Buffer_Send_FaceArea = new double[MaxLocalVertex_NearField];
  
  double *Buffer_Receive_Xcoord = new double[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_Ycoord = new double[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_Zcoord = new double[nProcessor*MaxLocalVertex_NearField];
  unsigned long *Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_Pressure = new double[nProcessor*MaxLocalVertex_NearField];
  double *Buffer_Receive_FaceArea = new double[nProcessor*MaxLocalVertex_NearField];
  
  unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_Ycoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_Zcoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
  unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
  unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;
  
  for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
    Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
    Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
    Buffer_Send_Ycoord[iVertex] = 0.0; Buffer_Send_Zcoord[iVertex] = 0.0;
  }
  
  /*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/
  
  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();
        
        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
            Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
            Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
            Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
            Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
            Buffer_Send_Pressure[nLocalVertex_NearField] = solver_container->node[iPoint]->GetPressure();
            Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nLocalVertex_NearField++;
          }
      }
  
  /*--- Send all the information --*/
  
  MPI_Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI_DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
  if (rank == MASTER_NODE) {
    
    Xcoord = new double[nVertex_NearField];
    Ycoord = new double[nVertex_NearField];
    Zcoord = new double[nVertex_NearField];
    AzimuthalAngle = new short[nVertex_NearField];
    IdPoint = new unsigned long[nVertex_NearField];
    IdDomain = new unsigned long[nVertex_NearField];
    Pressure = new double[nVertex_NearField];
    FaceArea = new double[nVertex_NearField];
    EquivArea = new double[nVertex_NearField];
    TargetArea = new double[nVertex_NearField];
    NearFieldWeight = new double[nVertex_NearField];
    Weight = new double[nVertex_NearField];
    
    nVertex_NearField = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
        Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField+iVertex];
        
        if (nDim == 2) {
          AzimuthalAngle[nVertex_NearField] = 0;
        }
        
        if (nDim == 3) {
          Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
          
          /*--- Rotate the nearfield cylinder  ---*/
          
          double YcoordRot = Ycoord[nVertex_NearField];
          double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);
          
          /*--- Compute the Azimuthal angle ---*/
          
          double AngleDouble; short AngleInt;
          AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);
          
          /*--- Fix an azimuthal line due to misalignments of the near-field ---*/
          
          double FixAzimuthalLine = config->GetFixAzimuthalLine();
          
          if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1))
            AngleDouble = FixAzimuthalLine - 0.1;
          
          AngleInt = (short) floor(AngleDouble + 0.5);
          
          if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
          else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
        }
        
        if (AzimuthalAngle[nVertex_NearField] <= 60) {
          IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
          Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
          FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
          IdDomain[nVertex_NearField] = iProcessor;
          nVertex_NearField++;
        }
        
      }
  }
  
  delete [] Buffer_Receive_nVertex;
  delete [] Buffer_Send_nVertex;
  
  delete [] Buffer_Send_Xcoord;
  delete [] Buffer_Send_Ycoord;
  delete [] Buffer_Send_Zcoord;
  delete [] Buffer_Send_IdPoint;
  delete [] Buffer_Send_Pressure;
  delete [] Buffer_Send_FaceArea;
  
  delete [] Buffer_Receive_Xcoord;
  delete [] Buffer_Receive_IdPoint;
  delete [] Buffer_Receive_Pressure;
  delete [] Buffer_Receive_FaceArea;
  
#endif
  
  if (rank == MASTER_NODE) {
    
    vector<short> PhiAngleList;
    vector<short>::iterator IterPhiAngleList;
    
    for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      PhiAngleList.push_back(AzimuthalAngle[iVertex]);
    
    sort( PhiAngleList.begin(), PhiAngleList.end());
    IterPhiAngleList = unique( PhiAngleList.begin(), PhiAngleList.end());
    PhiAngleList.resize( IterPhiAngleList - PhiAngleList.begin() );
    
    /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
    
    vector<vector<double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size());
    vector<vector<double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());
    
    /*--- Distribute the values among the different PhiAngles ---*/
    
    for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
          Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
          Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
          Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
          IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
          IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
          Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
          FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
          EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
          TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
          NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
          Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
        }
    
    /*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
        for (jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
          if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex+1]) {
            auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex+1]; Xcoord_PhiAngle[iPhiAngle][jVertex+1] = auxXCoord;
            auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex+1]; Ycoord_PhiAngle[iPhiAngle][jVertex+1] = auxYCoord;
            auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex+1]; Zcoord_PhiAngle[iPhiAngle][jVertex+1] = auxZCoord;
            auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex+1]; Pressure_PhiAngle[iPhiAngle][jVertex+1] = auxPress;
            auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex+1]; FaceArea_PhiAngle[iPhiAngle][jVertex+1] = auxArea;
            auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex+1]; IdPoint_PhiAngle[iPhiAngle][jVertex+1] = auxPoint;
            auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex+1]; IdDomain_PhiAngle[iPhiAngle][jVertex+1] = auxDomain;
          }
    
    
    /*--- Check that all the azimuth lists have the same size ---*/
    
    unsigned short nVertex = Xcoord_PhiAngle[0].size();
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      unsigned short nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
      if (nVertex_aux != nVertex) cout <<"Be careful!!! one azimuth list is shorter than the other"<< endl;
      nVertex = min(nVertex, nVertex_aux);
    }
    
    /*--- Compute equivalent area distribution at each azimuth angle ---*/
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
      for (iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
        
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex]*sin(AoA);
        
        for (jVertex = 0; jVertex < iVertex-1; jVertex++) {
          
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex]*sin(AoA);
          jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex+1]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex+1]*sin(AoA);
          
          jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
          jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);
          
          DeltaX = (jp1Coord-Coord_j);
          MeanFuntion = 0.5*(jp1Function + jFunction);
          EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFuntion;
        }
      }
    }
    
    /*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/
    
    NearFieldEA_file.precision(15);
    NearFieldEA_file.open("Equivalent_Area.dat", ios::out);
    NearFieldEA_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << endl;
    
    if (config->GetSystemMeasurements() == US)
      NearFieldEA_file << "VARIABLES = \"Height (in) at radial distance "<< R_Plane*12.0 << " (cylindrical coordinate system)\"";
    else
      NearFieldEA_file << "VARIABLES = \"Height (m) at radial distance "<< R_Plane << " (cylindrical coordinate system)\"";
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << ", \"Equivalent Area (ft<sup>2</sup>), Phi= " << PhiAngleList[iPhiAngle] << " deg.\"";
      else
        NearFieldEA_file << ", \"Equivalent Area (m<sup>2</sup>), Phi= " << PhiAngleList[iPhiAngle] << " deg.\"";
    }
    
    NearFieldEA_file << endl;
    for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {
      
      double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
      double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
      
      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
      else
        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init);
      
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
        NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
      }
      
      NearFieldEA_file << endl;
      
    }
    NearFieldEA_file.close();
    
    /*--- Read target equivalent area from the configuration file,
     this first implementation requires a complete table (same as the original
     EA table). so... no interpolation. ---*/
    
    vector<vector<double> > TargetArea_PhiAngle_Trans;
    TargetEA_file.open("TargetEA.dat", ios::in);
    
    if (TargetEA_file.fail()) {
      if (iExtIter == 0) { cout << "There is no Target Equivalent Area file (TargetEA.dat)!!"<< endl;
        cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
      }
      /*--- Set the table to 0 ---*/
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
    }
    else {
      
      /*--- skip header lines ---*/
      
      string line;
      getline(TargetEA_file, line);
      getline(TargetEA_file, line);
      
      while (TargetEA_file) {
        
        string line;
        getline(TargetEA_file, line);
        istringstream is(line);
        vector<double> row;
        unsigned short iter = 0;
        
        while (is.good()) {
          string token;
          getline(is,token,',');
          
          istringstream js(token);
          
          double data;
          js >> data;
          
          /*--- The first element in the table is the coordinate (in or m)---*/
          
          if (iter != 0) row.push_back(data);
          iter++;
          
        }
        TargetArea_PhiAngle_Trans.push_back(row);
      }
      
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];
      
    }
    
    /*--- Divide by the number of Phi angles in the nearfield ---*/
    
    double PhiFactor = 1.0/double(PhiAngleList.size());
    
    /*--- Evaluate the objective function ---*/
    
    InverseDesign = 0;
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        
        double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
        double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][iVertex]);
        
        if ((percentage < 0.1) || (Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;
        
        InverseDesign += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;
        
      }
    
    /*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
        for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
          Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
          
          double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
          double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][jVertex]);
          
          if ((percentage < 0.1) || (Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;
          
          NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
        }
      }
    
    /*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/
    
    EquivArea_file.precision(15);
    EquivArea_file.open("nearfield_flow.dat", ios::out);
    EquivArea_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << endl;
    
    if (config->GetSystemMeasurements() == US)
      EquivArea_file << "VARIABLES = \"Height (in) at radial distance "<< R_Plane*12.0 << " (cylindrical coordinate system)\",\"Equivalent Area (ft<sup>2</sup>)\",\"Target Equivalent Area (ft<sup>2</sup>)\",\"Cp\"" << endl;
    else
      EquivArea_file << "VARIABLES = \"Height (m) at radial distance "<< R_Plane << " (cylindrical coordinate system)\",\"Equivalent Area (m<sup>2</sup>)\",\"Target Equivalent Area (m<sup>2</sup>)\",\"Cp\"" << endl;
    
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_file << fixed << "ZONE T= \"Azimuthal angle " << PhiAngleList[iPhiAngle] << " deg.\"" << endl;
      for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {
        
        double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
        double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
        
        if (config->GetSystemMeasurements() == US)
          EquivArea_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
        else
          EquivArea_file << scientific << (XcoordRot - XcoordRot_init);
        
        EquivArea_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
        << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << endl;
      }
    }
    
    EquivArea_file.close();
    
    /*--- Write Weight file for adjoint computation ---*/
    
    FuncGrad_file.precision(15);
    FuncGrad_file.open("WeightNF.dat", ios::out);
    
    FuncGrad_file << scientific << "-1.0";
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      FuncGrad_file << scientific << "\t" << PhiAngleList[iPhiAngle];
    FuncGrad_file << endl;
    
    for (iVertex = 0; iVertex < NearFieldWeight_PhiAngle[0].size(); iVertex++) {
      double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
      FuncGrad_file << scientific << XcoordRot;
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        FuncGrad_file << scientific << "\t" << NearFieldWeight_PhiAngle[iPhiAngle][iVertex];
      FuncGrad_file << endl;
    }
    FuncGrad_file.close();
    
    /*--- Delete structures ---*/
    
    delete [] Xcoord; delete [] Ycoord; delete [] Zcoord;
    delete [] AzimuthalAngle; delete [] IdPoint; delete [] IdDomain;
    delete [] Pressure; delete [] FaceArea;
    delete [] EquivArea; delete [] TargetArea;
    delete [] NearFieldWeight; delete [] Weight;
    
  }
  
#ifndef HAVE_MPI
  
  /*--- Store the value of the NearField coefficient ---*/
  
  solver_container->SetTotal_CEquivArea(InverseDesign);
  
#else
  
  /*--- Send the value of the NearField coefficient to all the processors ---*/
  
  MPI_Bcast(&InverseDesign, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- Store the value of the NearField coefficient ---*/
  
  solver_container->SetTotal_CEquivArea(InverseDesign);
  
#endif
  
}

