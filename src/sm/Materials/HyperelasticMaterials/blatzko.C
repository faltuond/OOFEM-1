#include "blatzko.h"

namespace oofem {

    IRResultType BlatzKoMaterial::initializeFrom(InputRecord* ir) {
        return IRResultType();
    }

    void BlatzKoMaterial::give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) {

    }

    void BlatzKoMaterial::giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep) {

    }

    MaterialStatus* BlatzKoMaterial::CreateStatus(GaussPoint* gp) const {
        return nullptr;
    }
};