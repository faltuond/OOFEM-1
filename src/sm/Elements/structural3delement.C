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

#include "Elements/structural3delement.h"
#include "feinterpol3d.h"
#include "gausspoint.h"
#include "CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"


namespace oofem {
Structural3DElement :: Structural3DElement(int n, Domain *aDomain) :
  NLStructuralElement(n, aDomain), FbarElementExtensionInterface(aDomain), PressureFollowerLoadElementInterface(this),
    matRotation(false)
{
  cellGeometryWrapper = NULL;
}


IRResultType
Structural3DElement :: initializeFrom(InputRecord *ir)
{
    this->matRotation = ir->hasField(_IFT_Structural3DElement_materialCoordinateSystem);

    IRResultType result = FbarElementExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }
    
    return NLStructuralElement :: initializeFrom(ir);
}




void
Structural3DElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int li, int ui)
// Returns the [ 6 x (nno*3) ] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep));
    
    answer.resize(6, dNdx.giveNumberOfRows() * 3);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);

        answer.at(5, 3 * i - 2) = answer.at(4, 3 * i - 1) = dNdx.at(i, 3);
        answer.at(6, 3 * i - 2) = answer.at(4, 3 * i - 0) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 1) = answer.at(5, 3 * i - 0) = dNdx.at(i, 1);
    }
}



void
Structural3DElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, double alpha)
// Returns the [ 9 x (nno * 3) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz, dv/dz, du/dz, du/dy, dw/dy, dw/dx, dv/dx
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),  * this->giveCellGeometryWrapper(tStep, alpha));
    
    answer.resize(9, dNdx.giveNumberOfRows() * 3);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
    }

}


MaterialMode
Structural3DElement :: giveMaterialMode()
{
    return _3dMat;
}


void
Structural3DElement :: giveMaterialOrientationAt(FloatArray &x, FloatArray &y, FloatArray &z, const FloatArray &lcoords)
{
    if ( this->elemLocalCS.isNotEmpty() ) { // User specified orientation
        x.beColumnOf(this->elemLocalCS, 1);
        y.beColumnOf(this->elemLocalCS, 2);
        z.beColumnOf(this->elemLocalCS, 3);
    } else {
        ///@todo This is subject to change. I'm not sure which is the best way to define a local c.s.
        FloatMatrix jac;
        FloatArray help;
        this->giveInterpolation()->giveJacobianMatrixAt( jac, lcoords, FEIElementGeometryWrapper(this) );
        x.beColumnOf(jac, 1); // This is {dx/dxi, dy/dxi, dz/dxi}
        x.normalize();
        help.beColumnOf(jac, 2);
        z.beVectorProductOf(x, help); // Normal to the xi-eta plane.
        z.normalize();
        y.beVectorProductOf(z, x);
    }
}


void
Structural3DElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y, z;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, z, gp->giveNaturalCoordinates() );
        // Transform from global c.s. to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(5) * x(0) * x(1) + e(1) * x(1) * x(1) + e(4) * x(0) * x(2) + e(3) * x(1) * x(2) + e(2) * x(2) * x(2),
            e(0) * y(0) * y(0) + e(5) * y(0) * y(1) + e(1) * y(1) * y(1) + e(4) * y(0) * y(2) + e(3) * y(1) * y(2) + e(2) * y(2) * y(2),
            e(0) * z(0) * z(0) + e(5) * z(0) * z(1) + e(1) * z(1) * z(1) + e(4) * z(0) * z(2) + e(3) * z(1) * z(2) + e(2) * z(2) * z(2),
            2 * e(0) * y(0) * z(0) + e(4) * y(2) * z(0) + 2 * e(1) * y(1) * z(1) + e(3) * y(2) * z(1) + e(5) * ( y(1) * z(0) + y(0) * z(1) ) + ( e(4) * y(0) + e(3) * y(1) + 2 * e(2) * y(2) ) * z(2),
            2 * e(0) * x(0) * z(0) + e(4) * x(2) * z(0) + 2 * e(1) * x(1) * z(1) + e(3) * x(2) * z(1) + e(5) * ( x(1) * z(0) + x(0) * z(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * z(2),
            2 * e(0) * x(0) * y(0) + e(4) * x(2) * y(0) + 2 * e(1) * x(1) * y(1) + e(3) * x(2) * y(1) + e(5) * ( x(1) * y(0) + x(0) * y(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * y(2)
        };
        this->giveStructuralCrossSection()->giveRealStress_3d(s, gp, rotStrain, tStep);
        answer = {
            s(0) * x(0) * x(0) + 2 * s(5) * x(0) * y(0) + s(1) * y(0) * y(0) + 2 * ( s(4) * x(0) + s(3) * y(0) ) * z(0) + s(2) * z(0) * z(0),
            s(0) * x(1) * x(1) + 2 * s(5) * x(1) * y(1) + s(1) * y(1) * y(1) + 2 * ( s(4) * x(1) + s(3) * y(1) ) * z(1) + s(2) * z(1) * z(1),
            s(0) * x(2) * x(2) + 2 * s(5) * x(2) * y(2) + s(1) * y(2) * y(2) + 2 * ( s(4) * x(2) + s(3) * y(2) ) * z(2) + s(2) * z(2) * z(2),
            y(2) * ( s(5) * x(1) + s(1) * y(1) + s(3) * z(1) ) + x(2) * ( s(0) * x(1) + s(5) * y(1) + s(4) * z(1) ) + ( s(4) * x(1) + s(3) * y(1) + s(2) * z(1) ) * z(2),
            y(2) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(2) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(2),
            y(1) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(1) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(1)
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, e, tStep);
    }
}

void
Structural3DElement :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y, z;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, z, gp->giveNaturalCoordinates() );

        Q = {
            { x(0) * x(0), x(1) * x(1), x(2) * x(2), x(1) * x(2), x(0) * x(2), x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), y(2) * y(2), y(1) * y(2), y(0) * y(2), y(0) * y(1) },
            { z(0) * z(0), z(1) * z(1), z(2) * z(2), z(1) * z(2), z(0) * z(2), z(0) * z(1) },
            { 2 * y(0) * z(0), 2 * y(1) * z(1), 2 * y(2) * z(2), y(2) * z(1) + y(1) * z(2), y(2) * z(0) + y(0) * z(2), y(1) * z(0) + y(0) * z(1) },
            { 2 * x(0) * z(0), 2 * x(1) * z(1), 2 * x(2) * z(2), x(2) * z(1) + x(1) * z(2), x(2) * z(0) + x(0) * z(2), x(1) * z(0) + x(0) * z(1) },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 2 * x(2) * y(2), x(2) * y(1) + x(1) * y(2), x(2) * y(0) + x(0) * y(2), x(1) * y(0) + x(0) * y(1) }
        };
        answer.rotatedWith(Q, 't');
    }
}


void
Structural3DElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w};
}


FEICellGeometry *
Structural3DElement :: giveCellGeometryWrapper(TimeStep *tStep, double alpha)
{
    if ( !cellGeometryWrapper ) {
      if(this->giveDomain()->giveEngngModel()->giveFormulation() == AL) {
	cellGeometryWrapper = new FEIElementDeformedGeometryWrapper(this);
      } else {
        cellGeometryWrapper = new FEIElementGeometryWrapper(this);
      }
    }

    if(this->giveDomain()->giveEngngModel()->giveFormulation() == AL) {
      FEIElementDeformedGeometryWrapper *dgw = static_cast<FEIElementDeformedGeometryWrapper*> (cellGeometryWrapper);      
      if(tStep != NULL) {
	dgw->setTimeStep(tStep);
      }

      dgw->setAlpha(alpha);
    }

    return cellGeometryWrapper;
}



int 
Structural3DElement :: computeNumberOfDofs() 
{
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask; 
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize(); 
}


void Structural3DElement :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 6) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}



double 
Structural3DElement :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper()));

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}



double
Structural3DElement :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveLengthInDir(normalToCrackPlane);
}


void
Structural3DElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType modeType)
{
  if(!FbarFlag) {
    NLStructuralElement :: computeDeformationGradientVector(answer, gp, tStep, modeType);
  } else {
    FbarElementExtensionInterface :: computeFbarDeformationGradientVector(answer, gp, tStep, this, modeType);
  }
}


void
Structural3DElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  if(!FbarFlag) {
    NLStructuralElement :: computeStiffnessMatrix(answer, rMode, tStep);
  } else {
    FbarElementExtensionInterface :: computeFbarStiffnessMatrix(answer, rMode, tStep, this);
  }
}

void
Structural3DElement :: computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    FloatArray vF;
      if(!FbarFlag) {
	NLStructuralElement :: computeFirstPKStressVector(answer, gp, tStep);
      } else {
	double J0_J = computeFbarDeformationGradientVector(vF, gp, tStep, this, VM_Total);
	cs->giveFirstPKStresses(answer, gp, vF, tStep);
	answer.times(pow(J0_J,-2./3.));
      }
}

  

// Surface support
void
Structural3DElement :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    /* Returns the [ 3 x (nno*3) ] shape function matrix {N} of the receiver, 
     * evaluated at the given gp.
     * {u} = {N}*{a} gives the displacements at the integration point.
     */ 
          
    // Evaluate the shape functions at the position of the gp. 
    FloatArray N;
    static_cast< FEInterpolation3d* > ( this->giveInterpolation() )->
        surfaceEvalN( N, iSurf, sgp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );  
    answer.beNMatrixOf(N, 3);
}

void
Structural3DElement :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    IntArray nodes;
    const int ndofsn = 3;

     static_cast< FEInterpolation3d* > ( this->giveInterpolation() )->
        computeLocalSurfaceMapping(nodes, iSurf);

    answer.resize(nodes.giveSize() *3 );

    for ( int i = 1; i <= nodes.giveSize(); i++ ) {
        answer.at(i * ndofsn - 2) = nodes.at(i) * ndofsn - 2;
        answer.at(i * ndofsn - 1) = nodes.at(i) * ndofsn - 1;
        answer.at(i * ndofsn) = nodes.at(i) * ndofsn;
    }
}

double
Structural3DElement :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( static_cast< FEInterpolation3d* > ( this->giveInterpolation() )-> 
        surfaceGiveTransformationJacobian( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}



int
Structural3DElement :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *)
{
    OOFEM_ERROR("surface local coordinate system not supported");
    return 1;
}






// Edge support

void
Structural3DElement :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */
    IntArray eNodes;
    static_cast< FEInterpolation3d* > ( this->giveInterpolation() )->computeLocalEdgeMapping(eNodes,  iEdge);

    answer.resize( eNodes.giveSize() * 3 );
    for ( int i = 1; i <= eNodes.giveSize(); i++ ) {
        answer.at(i * 3 - 2) = eNodes.at(i) * 3 - 2;
        answer.at(i * 3 - 1) = eNodes.at(i) * 3 - 1;
        answer.at(i * 3)     = eNodes.at(i) * 3;
    }
}



double
Structural3DElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    /* Returns the line element ds associated with the given gp on the specific edge.
     * Note: The name is misleading since there is no volume to speak of in this case. 
     * The returned value is used for integration of a line integral (external forces).
     */
    double detJ = static_cast< FEInterpolation3d* > ( this->giveInterpolation() )->
        edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ * gp->giveWeight();
}


int
Structural3DElement :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    ///@todo how should this be supported 
    OOFEM_ERROR("egde local coordinate system not supported");
    return 1;
}

Interface *
Structural3DElement :: giveInterface(InterfaceType interface)
{
    if ( interface == PressureFollowerLoadElementInterfaceType) {
        return static_cast< PressureFollowerLoadElementInterface* >(this);
    }
    return NULL;
}

// support for pressure follower load interface
void
Structural3DElement ::  surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{

  FEInterpolation3d *interp3d = static_cast<FEInterpolation3d*> (this->giveInterpolation());
  interp3d->surfaceEvaldNdxi(answer, iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
}


void
Structural3DElement ::  surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{
  IntArray bNodes;
  FloatArray gcoords, lcoords, vU;
  FloatMatrix dNdxi, dxdxi, x;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundarySurfaceNodes (bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(nNodes,3);
  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v, D_w}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1) + vU.at( (i-1) * 3 + 1);
      x.at(i,2) = node->giveCoordinate(2) + vU.at( (i-1) * 3 + 2);
      x.at(i,3) = node->giveCoordinate(3) + vU.at( (i-1) * 3 + 3);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1);
      x.at(i,2) = node->giveCoordinate(2);
      x.at(i,3) = node->giveCoordinate(3);
    }
    
  }

  /* testing stuff
  FloatMatrix dNdk(4,1), xx, dxdx;
  for(int i = 1; i <= 4; i++) {
    dNdk.at(i,1) = dNdxi.at(i,1);
  }

  dxdx.beTProductOf(x, dNdk);
  */
  
  dxdxi.beTProductOf(dNdxi,x);
  dxdxi.copyRow(dxdksi, 2);
  dxdxi.copyRow(dxdeta, 1);
  answer.beVectorProductOf(dxdksi, dxdeta);

}



  

} // end namespace oofem
