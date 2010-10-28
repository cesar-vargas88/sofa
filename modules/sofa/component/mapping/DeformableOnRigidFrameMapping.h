/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_MAPPING_DEFORMABLEONRIGIDFRAME_H
#define SOFA_COMPONENT_MAPPING_DEFORMABLEONRIGIDFRAME_H

#include <sofa/core/Multi2Mapping.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <sofa/component/component.h>

#include <sofa/defaulttype/VecTypes.h>

#include <vector>

namespace sofa
{

namespace component
{

namespace mapping
{

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class InDataTypes, class OutDataTypes>
class DeformableOnRigidFrameMappingInternalData
{
public:
};

template <class TIn, class TInRoot, class TOut>
class DeformableOnRigidFrameMapping : public core::Multi2Mapping<TIn, TInRoot, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(DeformableOnRigidFrameMapping, TIn, TInRoot, TOut), SOFA_TEMPLATE3(core::Multi2Mapping, TIn, TInRoot, TOut));

    typedef core::Multi2Mapping<TIn, TInRoot, TOut> Inherit;

    typedef TIn In;
    typedef TInRoot InRoot;
    typedef TOut Out;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename OutCoord::value_type OutReal;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::Real InReal;
    typedef Data<InVecCoord> InDataVecCoord;
    typedef Data<InVecDeriv> InDataVecDeriv;
    typedef Data<InMatrixDeriv> InDataMatrixDeriv;

    typedef typename InRoot::VecCoord InRootVecCoord;
    typedef typename InRoot::VecDeriv InRootVecDeriv;
    typedef typename InRoot::MatrixDeriv InRootMatrixDeriv;
    typedef typename InRoot::Coord InRootCoord;
    typedef typename InRoot::Deriv InRootDeriv;
    typedef typename InRoot::Real InRootReal;
    typedef Data<InRootVecCoord> InRootDataVecCoord;
    typedef Data<InRootVecDeriv> InRootDataVecDeriv;
    typedef Data<InRootMatrixDeriv> InRootDataMatrixDeriv;

    typedef typename OutCoord::value_type Real;
    typedef OutCoord Coord;
    typedef OutDeriv Deriv;
    enum { N=Out::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;
    typedef defaulttype::Vec<N,Real> Vector ;

    OutVecCoord rotatedPoints;
    DeformableOnRigidFrameMappingInternalData<In, Out> data;
    Data<unsigned int> index;
    sofa::core::objectmodel::DataFileName fileDeformableOnRigidFrameMapping;
    Data< bool > useX0;
    Data< bool > indexFromEnd;
    Data<sofa::helper::vector<unsigned int> >  repartition;
    Data< bool > globalToLocalCoords;

    helper::ParticleMask* maskFrom;
    helper::ParticleMask* maskTo;

    DeformableOnRigidFrameMapping(helper::vector< core::State<In>* > from,
            helper::vector< core::State<InRoot>* > fromRoot,
            helper::vector< core::State<Out>* > to);

    virtual ~DeformableOnRigidFrameMapping()
    {}

    int addPoint ( const OutCoord& c );
    int addPoint ( const OutCoord& c, int indexFrom );

    void init();

    //Apply
    void apply( OutVecCoord& out, const InVecCoord& in, const InRootVecCoord* inroot  );
    void apply( const helper::vector<OutDataVecCoord*>& dataVecOutPos,
            const helper::vector<const InDataVecCoord*>& dataVecInPos ,
            const helper::vector<const InRootDataVecCoord*>& dataVecInRootPos,
            const core::MechanicalParams* /* mparams */)
    {
        if(dataVecOutPos.empty() || dataVecInPos.empty())
            return;

        const InRootVecCoord* inroot = NULL;

        //We need only one input In model and input Root model (if present)
        OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
        const InVecCoord& in = dataVecInPos[0]->getValue();

        if (!dataVecInRootPos.empty())
            inroot = &dataVecInRootPos[0]->getValue();

        apply(out, in, inroot);

        dataVecOutPos[0]->endEdit();
    }

    //ApplyJ
    void applyJ( OutVecDeriv& out, const InVecDeriv& in, const InRootVecDeriv* inroot );
    void applyJ( const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
            const helper::vector<const InDataVecDeriv*>& dataVecInVel,
            const helper::vector<const InRootDataVecDeriv*>& dataVecInRootVel,
            const core::MechanicalParams* /* mparams */)
    {
        if(dataVecOutVel.empty() || dataVecInVel.empty())
            return;

        const InRootVecDeriv* inroot = NULL;

        //We need only one input In model and input Root model (if present)
        OutVecDeriv& out = *dataVecOutVel[0]->beginEdit();
        const InVecDeriv& in = dataVecInVel[0]->getValue();

        if (!dataVecInRootVel.empty())
            inroot = &dataVecInRootVel[0]->getValue();

        applyJ(out,in, inroot);

        dataVecOutVel[0]->endEdit();
    }

    //ApplyJT Force
    void applyJT( InVecDeriv& out, const OutVecDeriv& in, InRootVecDeriv* outroot );
    void applyJT( const helper::vector< InDataVecDeriv*>& dataVecOutForce,
            const helper::vector< InRootDataVecDeriv*>& dataVecOutRootForce,
            const helper::vector<const OutDataVecDeriv*>& dataVecInForce,
            const core::MechanicalParams* /* mparams */)
    {
        if(dataVecOutForce.empty() || dataVecInForce.empty())
            return;

        InRootVecDeriv* outroot = NULL;

        //We need only one input In model and input Root model (if present)
        InVecDeriv& out = *dataVecOutForce[0]->beginEdit();
        const OutVecDeriv& in = dataVecInForce[0]->getValue();

        if (!dataVecOutRootForce.empty())
            outroot = dataVecOutRootForce[0]->beginEdit();

        applyJT(out,in, outroot);

        dataVecOutForce[0]->endEdit();
        if (outroot != NULL)
            dataVecOutRootForce[0]->endEdit();
    }

    //ApplyJT Constraint
    void applyJT( InMatrixDeriv& out, const OutMatrixDeriv& in, InRootMatrixDeriv* outroot );
    void applyJT( const helper::vector< InDataMatrixDeriv*>& dataMatOutConst ,
            const helper::vector< InRootDataMatrixDeriv*>&  dataMatOutRootConst ,
            const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst ,
            const core::MechanicalParams* /* mparams */)
    {
        if(dataMatOutConst.empty() || dataMatInConst.empty())
            return;

        InRootMatrixDeriv* outroot = NULL;

        //We need only one input In model and input Root model (if present)
        InMatrixDeriv& out = *dataMatOutConst[0]->beginEdit();
        const OutMatrixDeriv& in = dataMatInConst[0]->getValue();

        if (!dataMatOutRootConst.empty())
            outroot = dataMatOutRootConst[0]->beginEdit();

        applyJT(out,in, outroot);

        dataMatOutConst[0]->endEdit();
        if (outroot != NULL)
            dataMatOutRootConst[0]->endEdit();
    }
//
///**
//	 * @name
//	 */
//	//@{
//	/**
//	 * @brief
//	 */
//	void propagateX();
//
//	/**
//	 * @brief
//	 */
//	void propagateXfree();
//
//
//	/**
//	 * @brief
//	 */
//	void propagateV();
//
//	/**
//	 * @brief
//	 */
//	void propagateDx();
//
//	/**
//	 * @brief
//	 */
//	void accumulateForce();
//
//	/**
//	 * @brief
//	 */
//	void accumulateDf();
//
//	/**
//	 * @brief
//	 */
//	void accumulateConstraint();

    /**
      * @brief
      MAP the mass: this function recompute the rigid mass (gravity center position and inertia) of the object
          based on its deformed shape
      */
    void recomputeRigidMass() {}

    //@}

    void draw();

    void clear ( int reserve=0 );

    void setRepartition ( unsigned int value );
    void setRepartition ( sofa::helper::vector<unsigned int> values );


    //We have to overload canCreate & create because input2 can be empty
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        bool createResult = true;

        helper::vector< core::State<In>* > stin1;
        helper::vector< core::State<InRoot>* > stin2;
        helper::vector< core::State<Out>* > stout;

        createResult = sofa::core::objectmodel::VectorObjectRef::parseAll< core::State<In> >("input1", arg, stin1);
        createResult = createResult && sofa::core::objectmodel::VectorObjectRef::parseAll< core::State<Out> >("output", arg, stout);

        //If one of the parsing failed
        if(!createResult)
            return false;

        return core::BaseMapping::canCreate(obj, context, arg);
    }

    template<class T>
    static void create(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        bool createResult = true;
        helper::vector< core::State<In>* > stin1;
        helper::vector< core::State<InRoot>* > stin2;
        helper::vector< core::State<Out>* > stout;

        createResult = sofa::core::objectmodel::VectorObjectRef::parseAll< core::State<In> >("input1", arg, stin1);
        createResult = createResult && sofa::core::objectmodel::VectorObjectRef::parseAll< core::State<Out> >("output", arg, stout);

        //If one of the parsing failed
        if(!createResult)
            return;

        sofa::core::objectmodel::VectorObjectRef::parseAll< core::State<InRoot> >("input2", arg, stin2);

        obj = new T( stin1, stin2, stout);

        if (context)
            context->addObject(obj);

        if (arg)
            obj->parse(arg);
    }

protected:
    core::State<In>* m_fromModel;
    core::State<Out>* m_toModel;
    core::State<InRoot>* m_fromRootModel;

    class Loader;
    void load ( const char* filename );  /// SUPRESS ? ///
    //const VecCoord& getPoints();         /// SUPRESS ? ///
    InRootCoord rootX;
};


using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::ExtVec3fTypes;

#if defined(WIN32) && !defined(SOFA_COMPONENT_MAPPING_DEFORMABLEONRIGIDFRAME_CPP)  //// ATTENTION PB COMPIL WIN3Z
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
extern template class SOFA_COMPONENT_MAPPING_API DeformableOnRigidFrameMapping< Vec3dTypes, Rigid3dTypes, Vec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_COMPONENT_MAPPING_API DeformableOnRigidFrameMapping< Vec3fTypes, Rigid3fTypes, Vec3fTypes >;
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
