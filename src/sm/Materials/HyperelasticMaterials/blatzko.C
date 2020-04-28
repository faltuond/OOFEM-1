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

        //retrieve deformation from status and convert to matrix form
        StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>(this->giveStatus(gp));
        FloatArray vF = status->giveTempFVector();
        FloatMatrix F;
        F.beMatrixForm(vF);

        //invariants of deformation (invariants of C)
        double i2 = compute_I2_C_from_F(F);
        double i3 = compute_I3_C_from_F(F);

        //their first derivatives
        FloatArray di2dF, di3dF;
        compute_dI2_C_dF(di2dF, F);
        compute_dI3_C_dF(di3dF, F);

        //their second derivatives
        FloatMatrix d2i2dF2, d2i3dF2;
        compute_d2I2_C_dF2(d2i2dF2, F);
        compute_d2I3_C_dF2(d2i3dF2, F);

        //assembling stiffness
        //A = coeff1 * d2i2dF2
        //  + coeff2 * (di2dF x di3dF)
        //  + coeff3 * (di3dF x di3dF)
        //  + coeff4 * d2i3dF2

        double coeff1 = mu / (2 * i3);
        double coeff2 = -mu / (i3 * i3);
        double coeff3 = mu / (2 * sqrt(i3*i3*i3)) * (2 * i2 / sqrt(i3*i3*i3) - 0.5);
        double coeff4 = 0.5* mu * (1 / sqrt(i3) - i2 / (i3*i3));

        answer = d2i2dF2;
        answer.times(coeff1);

        FloatMatrix di2timesdi3;
        di2timesdi3.beDyadicProductOf(di2dF, di3dF);
        di2timesdi3.times(coeff2);
        answer.add(di2timesdi3);

        FloatMatrix di3timesdi3;
        di3timesdi3.beDyadicProductOf(di3dF, di3dF);
        di3timesdi3.times(coeff3);
        answer.add(di3timesdi3);

        d2i3dF2.times(coeff4);
        answer.add(d2i3dF2);

    }

    void BlatzKoMaterial::giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep) {
        
        //retrieve status
        StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>(this->giveStatus(gp));

        //convert F to matrix form
        FloatMatrix F;
        F.beMatrixForm(vF);

        //invariants of deformation (invariants of C)
        double i2 = compute_I2_C_from_F(F);
        double i3 = compute_I3_C_from_F(F);

        //their first derivatives
        FloatArray di2dF, di3dF;
        compute_dI2_C_dF(di2dF, F);
        compute_dI3_C_dF(di3dF, F);

        //assembling the stress vector
        // P = mu/(2I_3) * (dI2dF + (sqrt(I_3) - (I_2/I_3)) * dI3dF)

        answer = di3dF;

        double dI3coeff = -sqrt(i3) + (i2 / i3);
        answer.times(dI3coeff);

        answer.add(di2dF);

        double mainCoeff = mu / (2 * i3);
        answer.times(mainCoeff);

        //store into status
        status->letTempFVectorBe(vF);
        status->letTempPVectorBe(answer);
    }

    MaterialStatus* BlatzKoMaterial::CreateStatus(GaussPoint* gp) const {
        return new StructuralMaterialStatus(1, this->giveDomain(), gp);
    }
};