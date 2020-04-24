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

    }

    MaterialStatus* BlatzKoMaterial::CreateStatus(GaussPoint* gp) const {
        return new StructuralMaterialStatus(1, this->giveDomain(), gp);
    }
};