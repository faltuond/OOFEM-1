/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "../sm/CrossSections/simplecrosssection.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalmaterialextensioninterface.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"

namespace oofem {
REGISTER_CrossSection(SimpleElectroMechanicsCrossSection);



void
SimpleElectroMechanicsCrossSection :: give_FirstPKStressVector_ElectricalDisplacementVector(FloatArray &P, FloatArray &D, GaussPoint *gp, const FloatArray &F, const FloatArray &E, TimeStep *tStep)
{
    // This function returns the first Piola-Kirchoff stress in vector format and vector of electrical displacement
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    if ( mode == _3dMat ) {
      mat->give_FirstPKStressVector_ElectricalDisplacementVector_3d(P, D, gp, F, E, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}




void
SimpleElectroMechanicsCrossSection :: giveStiffnessMatrix_dPdF(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        mat->give3dMaterialStiffnessMatrix_dPdF(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
SimpleElectroMechanicsCrossSection :: giveStiffnessMatrix_dPdE(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        mat->give3dMaterialStiffnessMatrix_dPdE(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

  
void
SimpleElectroMechanicsCrossSection :: giveStiffnessMatrix_dDdE(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        mat->give3dMaterialStiffnessMatrix_dDdE(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

  

IRResultType
SimpleElectroMechanicsCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;

    double thick = 0.0;
    if ( ir->hasField(_IFT_SimpleElectroMechanicsCrossSection_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thick, _IFT_SimpleElectroMechanicsCrossSection_thick);
        propertyDictionary.add(CS_Thickness, thick);
    }

    double width = 0.0;
    if ( ir->hasField(_IFT_SimpleElectroMechanicsCrossSection_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_SimpleElectroMechanicsCrossSection_width);
        propertyDictionary.add(CS_Width, width);
    }

    double area = 0.0;
    if ( ir->hasField(_IFT_SimpleElectroMechanicsCrossSection_area) ) {
        IR_GIVE_FIELD(ir, area, _IFT_SimpleElectroMechanicsCrossSection_area);
    } else {
        area = thick * width;
    }
    propertyDictionary.add(CS_Area, area);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleElectroMechanicsCrossSection_iy);
    propertyDictionary.add(CS_InertiaMomentY, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleElectroMechanicsCrossSection_iz);
    propertyDictionary.add(CS_InertiaMomentZ, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleElectroMechanicsCrossSection_ik);
    propertyDictionary.add(CS_TorsionMomentX, value);

    this->materialNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleElectroMechanicsCrossSection_MaterialNumber);
  
    return CrossSection :: initializeFrom(ir);
}


void
SimpleElectroMechanicsCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    Material *mat = domain->giveMaterial(materialNumber);
    MaterialStatus *matStat = mat->CreateStatus(& iGP);
    iGP.setMaterialStatus(matStat);
}


bool
SimpleElectroMechanicsCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    if ( this->giveMaterialNumber() ) {
        return this->domain->giveMaterial( this->giveMaterialNumber() )->isCharacteristicMtrxSymmetric(rMode);
    } else {
        return false; // Bet false...
    }
}


Material
*SimpleElectroMechanicsCrossSection :: giveMaterial(IntegrationPoint *ip)
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


double
SimpleElectroMechanicsCrossSection :: give(int aProperty, GaussPoint *gp)
{
    return this->giveMaterial(gp)->give(aProperty, gp);
}


int
SimpleElectroMechanicsCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return this->giveMaterial(ip)->giveIPValue(answer, ip, type, tStep);
}



int
SimpleElectroMechanicsCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    Material *mat = this->giveDomain()->giveMaterial(this->materialNumber);
    if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
        OOFEM_WARNING( "material %s without structural support", mat->giveClassName() );
        result = 0;
    }

    return result;
}


Interface
*SimpleElectroMechanicsCrossSection :: giveMaterialInterface(InterfaceType t, IntegrationPoint *ip)
{
    return this->giveMaterial(ip)->giveInterface(t);
}



int
SimpleElectroMechanicsCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->packUnknowns(buff, tStep, gp);
}

int
SimpleElectroMechanicsCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
SimpleElectroMechanicsCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveMaterial(gp)->estimatePackSize(buff, gp);
}
} // end namespace oofem
