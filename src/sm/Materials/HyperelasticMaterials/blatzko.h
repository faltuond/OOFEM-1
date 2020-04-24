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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef mooneyrivlinmaterial_h
#define mooneyrivlinmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

 ///@name Input fields for MooneyRivlinMaterial
 //@{
#define _IFT_BlatzKoMaterial_Name "blatzkomat"
#define _IFT_BlatzKoMaterial_mu "mu"
//@}

namespace oofem {
    /**
     * This class implements the Blatz-Ko model for porous materials.
     *
     * @author Ondrej Faltus
     *
     */
    class BlatzKoMaterial : public StructuralMaterial
    {
    protected:
        // Material parameters
        double mu;

    public:
        BlatzKoMaterial(int n, Domain *d) : StructuralMaterial(n, d) { };

        virtual IRResultType initializeFrom(InputRecord *ir);

        virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
            MatResponseMode mode, GaussPoint *gp,
            TimeStep *tStep) {
            OOFEM_ERROR("not implemented, this material is designed for large strains only");
        }


        virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
            MatResponseMode mode,
            GaussPoint *gp, TimeStep *tStep);


        virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
        virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

        virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

        virtual const char *giveInputRecordName() const { return _IFT_BlatzKoMaterial_Name; }
        virtual const char *giveClassName() const { return "BlatzKoMaterial"; }
    };

} // end namespace oofem
#endif
