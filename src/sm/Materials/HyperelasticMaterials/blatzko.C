#include "blatzko.h"

namespace oofem {

    IRResultType BlatzKoMaterial::initializeFrom(InputRecord* ir) {
        IRResultType result;           // Required by IR_GIVE_FIELD macro

        result = StructuralMaterial::initializeFrom(ir);
        if ( result != IRRT_OK ) return result;

        IR_GIVE_FIELD(ir, mu, _IFT_BlatzKoMaterial_mu);

        return IRRT_OK;
    }

    void BlatzKoMaterial::give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) {
        
    }

    void BlatzKoMaterial::giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep) {

        //invariants of deformation (invariants of C)
        double i2 = compute_I2_C_from_F(vF);
        double i3 = compute_I3_C_from_F(vF);

        //their first derivatives
        FloatArray di2dF, di3dF;
        compute_dI2_C_dF(di2dF, vF);
        compute_dI3_C_dF(di3dF, vF);

        //assembling the stress vector
        // P = mu/(2I_3) * (dI2dF + (sqrt(I_3) - (I_2/I_3)) * dI3dF)
        
        answer = di3dF;

        double dI3coeff = - sqrt(i3) + (i2 / i3);
        answer.times(dI3coeff);

        answer.add(di2dF);
        
        double mainCoeff = mu / (2 * i3);
        answer.times(mainCoeff);
    }

    MaterialStatus* BlatzKoMaterial::CreateStatus(GaussPoint* gp) const {
        return new StructuralMaterialStatus(1, this->giveDomain(), gp);
    }
};