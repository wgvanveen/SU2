/*!
 * \file fluid_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
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

#pragma once

inline double CFluidModel::GetPressure () { return Pressure; }
inline double CFluidModel::GetSoundSpeed () { return sqrt(SoundSpeed2); }
inline double CFluidModel::GetSoundSpeed2 () { return SoundSpeed2; }
inline double CFluidModel::GetDensity () { return Density; }
inline double CFluidModel::GetEntropy () { return Entropy; }
inline double CFluidModel::GetStaticEnergy () { return StaticEnergy; }
inline double CFluidModel::GetTemperature () { return Temperature; }
inline double CFluidModel::GetdPdrho_e () { return dPdrho_e; }
inline double CFluidModel::GetdPde_rho () { return dPde_rho; }
inline double CFluidModel::GetdTdrho_e () { return dTdrho_e; }
inline double CFluidModel::GetdTde_rho () { return dTde_rho; }

inline double CFluidModel::GetLaminarViscosity (double T, double rho) {
        DynamicViscosity->SetViscosity(T, rho);
        return DynamicViscosity->GetViscosity();
}
inline double CFluidModel::Getdmudrho_T () {
        return DynamicViscosity->Getdmudrho_T();
}
inline double CFluidModel::GetdmudT_rho () {
        return DynamicViscosity->GetdmudT_rho();
}

inline double CFluidModel::GetThermalConductivity (double par1, double par2) {
        ThermalConductivity->SetThermalConductivity(par1, par2);
        return ThermalConductivity->GetThermalConductivity();
}
inline double CFluidModel::Getdktdrho_T () {
        return ThermalConductivity->GetDerThermalConductivity_rho_T();
}
inline double CFluidModel::GetdktdT_rho () {
        return ThermalConductivity->GetDerThermalConductivity_T_rho();
}

inline void CFluidModel::SetTDState_rhoe (double rho, double e ) { }
inline void CFluidModel::SetTDState_PT (double P, double T ) { }
inline void CFluidModel::SetTDState_Prho (double P, double rho ) { }
inline void CFluidModel::SetTDState_hs (double h, double s ) { }
inline void CFluidModel::SetTDState_rhoT (double rho, double T ) { }
inline void CFluidModel::SetEnergy_Prho (double P, double rho ) { }
