/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program. If not, see <http://www.gnu.org/licenses/>.              *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/


#include <sstream>
using std::ostringstream;
#include <fstream>

#include <string>
using std::string;

#include <vector>

#include <runSofaValidation.h>

#include <sofa/helper/ArgumentParser.h>
#include <SofaSimulationCommon/config.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/simulation/config.h> 
// #defines SOFA_HAVE_DAG (or not)
#include <SofaSimulationCommon/init.h>

#ifdef SOFA_HAVE_DAG

#include <SofaSimulationGraph/init.h>
#include <SofaSimulationGraph/DAGSimulation.h>

#endif

#include <SofaSimulationTree/init.h>
#include <SofaSimulationTree/TreeSimulation.h>
using sofa::simulation::Node;
//#include <sofa/simulation/SceneLoaderFactory.h>
#include <SofaGraphComponent/SceneCheckerListener.h>
using sofa::simulation::scenechecking::SceneCheckerListener;

#include <SofaCommon/initSofaCommon.h>
#include <SofaBase/initSofaBase.h>
#include <SofaGeneral/initSofaGeneral.h>
#include <SofaMisc/initSofaMisc.h>

//#include <SofaGeneralLoader/ReadState.h>
//#include <sofa/helper/Factory.h>
//#include <sofa/helper/cast.h>
#include <sofa/helper/BackTrace.h>
//#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/FileSystem.h>
using sofa::helper::system::FileSystem;
//#include <sofa/helper/system/SetDirectory.h>
#include <sofa/helper/Utils.h>
#include <sofa/gui/GUIManager.h>
using sofa::gui::GUIManager;

#include <sofa/gui/Main.h>
#include <sofa/gui/BatchGUI.h>  // For the default number of iterations
//#include <sofa/helper/system/gl.h>

//using sofa::core::ExecParams ;

#include <sofa/helper/system/console.h>
using sofa::helper::Utils;

//using sofa::component::misc::ReadStateActivator;
using sofa::simulation::tree::TreeSimulation;
using sofa::simulation::graph::DAGSimulation;
using sofa::helper::system::SetDirectory;
//using sofa::core::objectmodel::BaseNode ;
//using sofa::gui::BatchGUI;
using sofa::gui::BaseGUI;

//#include <sofa/helper/logging/Messaging.h>

#include <sofa/helper/logging/ConsoleMessageHandler.h>
using sofa::helper::logging::ConsoleMessageHandler;

#include <sofa/core/logging/RichConsoleStyleMessageFormatter.h>
using  sofa::helper::logging::RichConsoleStyleMessageFormatter;

#include <sofa/core/logging/PerComponentLoggingMessageHandler.h>
using  sofa::helper::logging::MainPerComponentLoggingMessageHandler;



#include <sofa/helper/AdvancedTimer.h>

#include <sofa/gui/GuiDataRepository.h>
using sofa::gui::GuiDataRepository;

using sofa::helper::system::DataRepository;

using sofa::helper::system::PluginRepository;
using sofa::helper::system::PluginManager;

#include <sofa/helper/logging/MessageDispatcher.h>
using sofa::helper::logging::MessageDispatcher;

#include <sofa/helper/logging/ClangMessageHandler.h>
using sofa::helper::logging::ClangMessageHandler;

#include <sofa/helper/logging/ExceptionMessageHandler.h>
using sofa::helper::logging::ExceptionMessageHandler;




//////////////////////
//      Scene       //
//////////////////////

//using Coord3 = sofa::defaulttype::Vector3;
//using VecCoord3 = sofa::helper::vector<Coord3>;


#include <sofa/defaulttype/VecTypes.h>
using sofa::defaulttype::StdVectorTypes;
using sofa::defaulttype::Vec;
using sofa::defaulttype::Vec3Types;


using sofa::core::objectmodel::New;
using namespace sofa;


// collision pipeline

#include <SofaBaseCollision/DefaultPipeline.h>
using sofa::component::collision::DefaultPipeline;

#include <SofaBaseCollision/BruteForceDetection.h>
using sofa::component::collision::BruteForceDetection;

#include <SofaBaseCollision/NewProximityIntersection.h>
using sofa::component::collision::NewProximityIntersection;

#include <SofaBaseCollision/DefaultContactManager.h>
using sofa::component::collision::DefaultContactManager;

#include <SofaMeshCollision/TriangleModel.h>
using TriangleCollisionModel3 = sofa::component::collision::TriangleCollisionModel<Vec3Types>;

#include <SofaMeshCollision/LineModel.h>
using LineCollisionModel3 = sofa::component::collision::LineCollisionModel<Vec3Types>;

#include <SofaMeshCollision/PointModel.h>
using PointCollisionModel3 = sofa::component::collision::PointCollisionModel<Vec3Types>;

#include <SofaValidation/EvalPointsDistance.h>
using EvalPointsDistance3 = sofa::component::misc::EvalPointsDistance<Vec3Types>;

// solvers
#include <SofaBaseLinearSolver/CGLinearSolver.h>
using sofa::component::linearsolver::CGLinearSolver;
using sofa::component::linearsolver::GraphScatteredMatrix;
using sofa::component::linearsolver::GraphScatteredVector;

#include <SofaImplicitOdeSolver/EulerImplicitSolver.h>
using sofa::component::odesolver::EulerImplicitSolver;

#include <SofaConstraint/GenericConstraintSolver.h>
using sofa::component::constraintset::GenericConstraintSolver;

#include <SofaConstraint/GenericConstraintCorrection.h>
using sofa::component::constraintset::GenericConstraintCorrection;

#include <SofaSparseSolver/SparseLDLSolver.h>
using namespace sofa::component::linearsolver;
using SparseLDLSolverGraph = SparseLDLSolver<CompressedRowSparseMatrix<double>, FullVector<double>>;

// mechanical object

//#include <SofaGraphComponent/Gravity.h>

#include <SofaBaseMechanics/MechanicalObject.h>
using MechanicalObject3 = sofa::component::container::MechanicalObject<Vec3Types>;

#include <SofaBaseMechanics/UniformMass.h>
using UniformMass3 = sofa::component::mass::UniformMass<Vec3Types, SReal>;

#include <SofaBoundaryCondition/ConstantForceField.h>
using sofa::component::forcefield::ConstantForceField;

#include <SofaSimpleFem/HexahedronFEMForceField.h>
using HexahedronFEMForceField3 = sofa::component::forcefield::HexahedronFEMForceField<Vec3Types>;

#include <SofaGeneralSimpleFem/TetrahedralCorotationalFEMForceField.h>
using TetrahedralCorotationalFEMForceField3 = sofa::component::forcefield::TetrahedralCorotationalFEMForceField<Vec3Types>;

#include <SofaSimpleFem/TetrahedronFEMForceField.h>
using TetrahedronFEMForceField3 = sofa::component::forcefield::TetrahedronFEMForceField<Vec3Types>;

#include <SofaBoundaryCondition/FixedConstraint.h>
using FixedConstraint3 = sofa::component::projectiveconstraintset::FixedConstraint< StdVectorTypes<Vec<3, double>, Vec<3, double>, double>>;

#include <SofaBaseMechanics/BarycentricMapping.h>
using BarycentricMapping3 = sofa::component::mapping::BarycentricMapping<Vec3Types, Vec3Types>;

#include <SofaBaseCollision/DiscreteIntersection.h>


//////////////////////////////////////

#include <SofaBaseTopology/SparseGridTopology.h>
using sofa::component::topology::SparseGridTopology;

#include <SofaGeneralLoader/MeshGmshLoader.h>
#include <SofaLoader/MeshObjLoader.h>
using namespace sofa::component::loader;

#include <SofaBaseTopology/MeshTopology.h>
using sofa::component::topology::MeshTopology;

#include <SofaOpenglVisual/OglModel.h>
using sofa::component::visualmodel::OglModel;

#include <SofaBaseVisual/VisualStyle.h>
using sofa::component::visualmodel::VisualStyle;

#include <sofa/core/visual/DisplayFlags.h>
using sofa::core::visual::DisplayFlags;

#include <SofaGeneralLoader/SphereLoader.h>
using SphereCollisionModel3 = sofa::component::collision::SphereCollisionModel<Vec3Types>;

#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
using namespace sofa::component::topology;

#include <SofaGeneralSimpleFem/TetrahedralCorotationalFEMForceField.h>
using TetrahedralCorotationalFEMForceField3 = sofa::component::forcefield::TetrahedralCorotationalFEMForceField<Vec3Types>;

#include <sofa/defaulttype/VecTypes.h>
using sofa::defaulttype::Vec3Types;

#include <SofaBaseMechanics/SubsetMapping.h>
using SubsetMapping3 = sofa::component::mapping::SubsetMapping<Vec3Types, Vec3Types>;

#include <SofaOpenglVisual/DataDisplay.h>
using sofa::component::visualmodel::DataDisplay;

#include <sofa/helper/ColorMap.h>
using sofa::helper::ColorMap;

#include <SofaBaseCollision/MinProximityIntersection.h>

#include <SofaConstraint/FreeMotionAnimationLoop.h>
using sofa::component::animationloop::FreeMotionAnimationLoop;

#include <sofa/simulation/DefaultAnimationLoop.h>
using sofa::simulation::DefaultAnimationLoop;

#include <SofaBaseTopology/RegularGridTopology.h>



#include <SofaBaseMechanics/IdentityMapping.h>
using IdentityMapping3 = sofa::component::mapping::IdentityMapping<Vec3Types, Vec3Types>;

#include <SofaOpenglVisual/OglColorMap.h>
using sofa::component::visualmodel::OglColorMap;


#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
using TetrahedronSetGeometryAlgorithms3 = sofa::component::topology::TetrahedronSetGeometryAlgorithms<Vec3Types>;


#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>
using sofa::component::topology::HexahedronSetTopologyContainer;

#include <SofaBaseTopology/HexahedronSetTopologyModifier.h>
using sofa::component::topology::HexahedronSetTopologyModifier;

#include <SofaBaseTopology/HexahedronSetTopologyAlgorithms.h>
using HexahedronSetTopologyAlgorithms3 = sofa::component::topology::HexahedronSetTopologyAlgorithms<Vec3Types>;

#include <SofaBaseTopology/HexahedronSetGeometryAlgorithms.h>
using HexahedronSetGeometryAlgorithms3 = sofa::component::topology::HexahedronSetGeometryAlgorithms<Vec3Types>;



#include <SofaBaseTopology/QuadSetTopologyContainer.h>
using sofa::component::topology::QuadSetTopologyContainer;

#include <SofaBaseTopology/QuadSetTopologyModifier.h>
using sofa::component::topology::QuadSetTopologyModifier;

#include <SofaBaseTopology/QuadSetTopologyAlgorithms.h>
using QuadSetTopologyAlgorithms3 = sofa::component::topology::QuadSetTopologyAlgorithms<Vec3Types>;

#include <SofaBaseTopology/QuadSetGeometryAlgorithms.h>
using QuadSetGeometryAlgorithms3 = sofa::component::topology::QuadSetGeometryAlgorithms<Vec3Types>;


#include <SofaTopologyMapping/Hexa2QuadTopologicalMapping.h>
using sofa::component::topology::Hexa2QuadTopologicalMapping;


// CUDA
//#include <SofaCUDA/sofa/gpu/cuda/CudaMechanicalObject.h>
//using CudaMechanicalObject3 = sofa::component::container::MechanicalObject<Vec3Types>;



/// helper for more compact component creation
template<class Component> typename Component::SPtr addNew(Node::SPtr parentNode, std::string name = "")
{
	typename Component::SPtr component = New<Component>();
	parentNode->addObject(component);
	component->setName(name);
	return component;
}

/// Create musculoskeletic system
simulation::Node::SPtr createScene(std::string project_path)
{
	#pragma region scene

	std::cout << "Project path: " << project_path << std::endl;
	std::string mesh_path = "mesh";// +project_path[project_path.length() - 1];
	mesh_path+= project_path[project_path.length() - 1];
	std::cout << "Mesh path: " << project_path + mesh_path << std::endl;

	Node::SPtr root = simulation::getSimulation()->createNewGraph("root");
	root->setGravity(sofa::defaulttype::Vector3(0, -9.81, 0));
	root->setAnimate(false);
	root->setDt(0.02);

	VisualStyle::SPtr visualStyle = addNew<VisualStyle>(root);
	visualStyle->setName("VisualStyle");
	DisplayFlags& flags = *visualStyle->displayFlags.beginEdit();
	flags.setShowVisualModels(true);
	flags.setShowBehaviorModels(true);
	flags.setShowForceFields(true);
	flags.setShowInteractionForceFields(false);
	flags.setShowCollisionModels(true);
	flags.setShowBoundingCollisionModels(true);
	flags.setShowMappings(false);
	flags.setShowMechanicalMappings(false);
	flags.setShowWireFrame(false);
	flags.setShowNormals(false);
	visualStyle->displayFlags.endEdit();

	FreeMotionAnimationLoop::SPtr freeMotionAnimationLoop = New<FreeMotionAnimationLoop>(root.get());
	freeMotionAnimationLoop->setName("FreeMotionAnimationLoop");
	root->addObject(freeMotionAnimationLoop);

	GenericConstraintSolver::SPtr genericConstraintSolver = addNew<GenericConstraintSolver>(root);
	genericConstraintSolver->setName("GenericConstraintSolver");
	genericConstraintSolver->maxIt.setValue(200);
	genericConstraintSolver->tolerance.setValue(1.0e-3);

	DefaultPipeline::SPtr defaultPipeline = addNew<DefaultPipeline>(root);
	defaultPipeline->setName("DefaultPipeline");

	BruteForceDetection::SPtr bruteForceDetection = addNew<BruteForceDetection>(root);
	bruteForceDetection->setName("BruteForceDetection");

	component::collision::MinProximityIntersection::SPtr minProximityIntersection = addNew<component::collision::MinProximityIntersection>(root);
	minProximityIntersection->setName("MinProximityIntersection");
	minProximityIntersection->setAlarmDistance(1.2);
	minProximityIntersection->setContactDistance(0.1);

	DefaultContactManager::SPtr defaultContactManager = addNew<DefaultContactManager>(root);
	defaultContactManager->setName("DefaultContactManager");
	defaultContactManager->setDefaultResponseType("FrictionContact");

	#pragma endregion

	#pragma region tp

	Node::SPtr tp = root->createChild("tp");

	EulerImplicitSolver::SPtr tp_EulerImplicitSolver = addNew<EulerImplicitSolver>(tp);
	SparseLDLSolverGraph::SPtr tp_SparseLDLSolverGraph = addNew<SparseLDLSolverGraph>(tp);

	MeshGmshLoader::SPtr tp_MeshGmshLoader = addNew<MeshGmshLoader>(tp);
	tp_MeshGmshLoader->setName("loader");
	tp_MeshGmshLoader->setFilename(project_path + mesh_path + "tp.msh");
	tp_MeshGmshLoader->load();

	TetrahedronSetTopologyContainer::SPtr tp_TetrahedronSetTopologyContainer = addNew<TetrahedronSetTopologyContainer>(tp);
	tp_TetrahedronSetTopologyContainer->setSrc("@loader");

	MechanicalObject3::SPtr tp_MechanicalObject = addNew<MechanicalObject3>(tp);
	tp_MechanicalObject->setTranslation(-5, 5, 0);
	tp_MechanicalObject->setRotation(0, 90, -90);
	tp_MechanicalObject->setScale(10, 10, 10);

	TetrahedronSetGeometryAlgorithms3::SPtr tp_TetrahedronSetGeometryAlgorithms = addNew<TetrahedronSetGeometryAlgorithms3>(tp);

	UniformMass3::SPtr tp_uniformMass = addNew<UniformMass3>(tp);
	tp_uniformMass->setTotalMass(2);

	TetrahedronFEMForceField3::SPtr tp_TetrahedronFEMForceField = addNew<TetrahedronFEMForceField3>(tp);
	tp_TetrahedronFEMForceField->setMethod(HexahedronFEMForceField3::SMALL);
	tp_TetrahedronFEMForceField->setPoissonRatio(0.495);
	tp_TetrahedronFEMForceField->setYoungModulus(3000);
	tp_TetrahedronFEMForceField->_computeVonMisesStress.setValue(1);
	tp_TetrahedronFEMForceField->_updateStiffness.setValue(true);
	tp_TetrahedronFEMForceField->_showVonMisesStressPerNode.setValue(false);
	tp_TetrahedronFEMForceField->_showStressColorMap.setValue("green  0   0.5 0   1  0.5  1   0   1");

	GenericConstraintCorrection::SPtr tp_GenericConstraintCorrection = addNew<GenericConstraintCorrection>(tp);

	FixedConstraint3::SPtr tp_FixedConstraints = addNew<FixedConstraint3>(tp);
	tp_FixedConstraints->addConstraint(3);
	tp_FixedConstraints->addConstraint(39);
	tp_FixedConstraints->addConstraint(64);

		#pragma region tp_surface

		Node::SPtr tp_surface = tp->createChild("surface");
		/*
		MeshObjLoader::SPtr tp_surface_MeshObjLoader = addNew<MeshObjLoader>(tp_surface);
		tp_surface_MeshObjLoader->setName("loader");
		tp_surface_MeshObjLoader->setFilename(project_path + mesh_path + "tp.obj");
		tp_surface_MeshObjLoader->load();

		MeshTopology::SPtr tp_surface_MeshTopology = addNew<MeshTopology>(tp_surface);
		tp_surface_MeshTopology->setSrc("@loader");

		MechanicalObject3::SPtr tp_surface_MechanicalObject = addNew<MechanicalObject3>(tp_surface);
		tp_surface_MechanicalObject->setTranslation(-5, 5, 0);
		tp_surface_MechanicalObject->setRotation(0, 90, -90);
		tp_surface_MechanicalObject->setScale(10, 10, 10);

		BarycentricMapping3::SPtr tp_surface_BarycentricMapping = addNew<BarycentricMapping3>(tp_surface);
		tp_surface_BarycentricMapping->setPathInputObject("@..");
		tp_surface_BarycentricMapping->setPathOutputObject("@.");

		TriangleCollisionModel3::SPtr tp_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(tp_surface);
		LineCollisionModel3::SPtr tp_surface_LineCollisionModel = addNew<LineCollisionModel3>(tp_surface);
		PointCollisionModel3::SPtr tp_surface_PointCollisionModel = addNew<PointCollisionModel3>(tp_surface);
		*/
		
		HexahedronSetTopologyContainer::SPtr tp_surface_HexahedronSetTopologyContainer = addNew<HexahedronSetTopologyContainer>(tp_surface);
		tp_surface_HexahedronSetTopologyContainer->setName("Container1");
		tp_surface_HexahedronSetTopologyContainer->setSrc("@../loader");
		HexahedronSetTopologyModifier::SPtr tp_surface_HexahedronSetTopologyModifier = addNew<HexahedronSetTopologyModifier>(tp_surface);
		HexahedronSetTopologyAlgorithms3::SPtr tp_surface_HexahedronSetTopologyAlgorithms = addNew<HexahedronSetTopologyAlgorithms3>(tp_surface);
		HexahedronSetGeometryAlgorithms3::SPtr tp_surface_HexahedronSetGeometryAlgorithms = addNew<HexahedronSetGeometryAlgorithms3>(tp_surface);

		QuadSetTopologyContainer::SPtr tp_surface_QuadSetTopologyContainer = addNew<QuadSetTopologyContainer>(tp_surface);
		tp_surface_QuadSetTopologyContainer->setName("Container2");
		QuadSetTopologyModifier::SPtr tp_surface_QuadSetTopologyModifier = addNew<QuadSetTopologyModifier>(tp_surface);
		QuadSetTopologyAlgorithms3::SPtr tp_surface_QuadSetTopologyAlgorithms = addNew<QuadSetTopologyAlgorithms3>(tp_surface);
		QuadSetGeometryAlgorithms3::SPtr tp_surface_QuadSetGeometryAlgorithms = addNew<QuadSetGeometryAlgorithms3>(tp_surface);

		Hexa2QuadTopologicalMapping::SPtr tp_surface_Hexa2QuadTopologicalMapping = addNew<Hexa2QuadTopologicalMapping>(tp_surface);
		tp_surface_Hexa2QuadTopologicalMapping->setPathInputObject("@Container1");
		tp_surface_Hexa2QuadTopologicalMapping->setPathOutputObject("@Container2");

		TriangleCollisionModel3::SPtr tp_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(tp_surface);
		LineCollisionModel3::SPtr tp_surface_LineCollisionModel = addNew<LineCollisionModel3>(tp_surface);
		PointCollisionModel3::SPtr tp_surface_PointCollisionModel = addNew<PointCollisionModel3>(tp_surface);

		#pragma endregion

	#pragma endregion


	#pragma region piper

	Node::SPtr piper = root->createChild("piper");

	EulerImplicitSolver::SPtr piper_EulerImplicitSolver = addNew<EulerImplicitSolver>(piper);
	SparseLDLSolverGraph::SPtr piper_SparseLDLSolverGraph = addNew<SparseLDLSolverGraph>(piper);

	MeshGmshLoader::SPtr piper_MeshGmshLoader = addNew<MeshGmshLoader>(piper);
	piper_MeshGmshLoader->setName("loader");

	piper_MeshGmshLoader->setFilename(project_path + mesh_path + "piper.msh");
	piper_MeshGmshLoader->load();

	TetrahedronSetTopologyContainer::SPtr piper_TetrahedronSetTopologyContainer = addNew<TetrahedronSetTopologyContainer>(piper);
	piper_TetrahedronSetTopologyContainer->setSrc("@loader");

	MechanicalObject3::SPtr piper_MechanicalObject = addNew<MechanicalObject3>(piper);
	piper_MechanicalObject->setTranslation(7, 1, 0);
	piper_MechanicalObject->setRotation(0, 0, 0);
	piper_MechanicalObject->setScale(1, 1, 1);

	TetrahedronSetGeometryAlgorithms3::SPtr piper_TetrahedronSetGeometryAlgorithms = addNew<TetrahedronSetGeometryAlgorithms3>(piper);

	UniformMass3::SPtr piper_uniformMass = addNew<UniformMass3>(piper);
	piper_uniformMass->setTotalMass(10);

	TetrahedronFEMForceField3::SPtr piper_TetrahedronFEMForceField = addNew<TetrahedronFEMForceField3>(piper);
	piper_TetrahedronFEMForceField->setMethod(HexahedronFEMForceField3::SMALL);
	piper_TetrahedronFEMForceField->setPoissonRatio(0.495);
	piper_TetrahedronFEMForceField->setYoungModulus(3000);
	piper_TetrahedronFEMForceField->_computeVonMisesStress.setValue(1);
	piper_TetrahedronFEMForceField->_updateStiffness.setValue(true);
	piper_TetrahedronFEMForceField->_showVonMisesStressPerNode.setValue(false);
	piper_TetrahedronFEMForceField->_showStressColorMap.setValue("green  0   0.5 0   1  0.5  1   0   1");

	GenericConstraintCorrection::SPtr piper_GenericConstraintCorrection = addNew<GenericConstraintCorrection>(piper);

	FixedConstraint3::SPtr piper_FixedConstraints = addNew<FixedConstraint3>(piper);
	piper_FixedConstraints->addConstraint(3);
	piper_FixedConstraints->addConstraint(39);
	piper_FixedConstraints->addConstraint(64);

		#pragma region piper_surface

		Node::SPtr piper_surface = piper->createChild("surface");

		HexahedronSetTopologyContainer::SPtr piper_surface_HexahedronSetTopologyContainer = addNew<HexahedronSetTopologyContainer>(piper_surface);
		piper_surface_HexahedronSetTopologyContainer->setName("Container1");
		piper_surface_HexahedronSetTopologyContainer->setSrc("@../loader");
		HexahedronSetTopologyModifier::SPtr piper_surface_HexahedronSetTopologyModifier = addNew<HexahedronSetTopologyModifier>(piper_surface);
		HexahedronSetTopologyAlgorithms3::SPtr piper_surface_HexahedronSetTopologyAlgorithms = addNew<HexahedronSetTopologyAlgorithms3>(piper_surface);
		HexahedronSetGeometryAlgorithms3::SPtr piper_surface_HexahedronSetGeometryAlgorithms = addNew<HexahedronSetGeometryAlgorithms3>(piper_surface);

		QuadSetTopologyContainer::SPtr piper_surface_QuadSetTopologyContainer = addNew<QuadSetTopologyContainer>(piper_surface);
		piper_surface_QuadSetTopologyContainer->setName("Container2");
		QuadSetTopologyModifier::SPtr piper_surface_QuadSetTopologyModifier = addNew<QuadSetTopologyModifier>(piper_surface);
		QuadSetTopologyAlgorithms3::SPtr piper_surface_QuadSetTopologyAlgorithms = addNew<QuadSetTopologyAlgorithms3>(piper_surface);
		QuadSetGeometryAlgorithms3::SPtr piper_surface_QuadSetGeometryAlgorithms = addNew<QuadSetGeometryAlgorithms3>(piper_surface);

		Hexa2QuadTopologicalMapping::SPtr piper_surface_Hexa2QuadTopologicalMapping = addNew<Hexa2QuadTopologicalMapping>(piper_surface);
		piper_surface_Hexa2QuadTopologicalMapping->setPathInputObject("@Container1");
		piper_surface_Hexa2QuadTopologicalMapping->setPathOutputObject("@Container2");

		TriangleCollisionModel3::SPtr piper_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(piper_surface);
		LineCollisionModel3::SPtr piper_surface_LineCollisionModel = addNew<LineCollisionModel3>(piper_surface);
		PointCollisionModel3::SPtr piper_surface_PointCollisionModel = addNew<PointCollisionModel3>(piper_surface);

		#pragma endregion

	#pragma endregion 


	#pragma region mattress

	Node::SPtr mattress = root->createChild("mattress");

	EulerImplicitSolver::SPtr mattress_EulerImplicitSolver = addNew<EulerImplicitSolver>(mattress);
	SparseLDLSolverGraph::SPtr mattress_SparseLDLSolverGraph = addNew<SparseLDLSolverGraph>(mattress);

	RegularGridTopology::SPtr mattress_RegularGridTopology = addNew<RegularGridTopology>(mattress);
	mattress_RegularGridTopology->setSize(3, 3, 3);
	mattress_RegularGridTopology->setPos(-1, 1, 2, 2.2, -2, 2);

	MechanicalObject3::SPtr mattress_MechanicalObject = addNew<MechanicalObject3>(mattress);
	mattress_MechanicalObject->setTranslation(0,0,0);
	mattress_MechanicalObject->setRotation(0,0,0);
	mattress_MechanicalObject->setScale(1,1,1);

	UniformMass3::SPtr mattress_uniformMass = addNew<UniformMass3>(mattress);
	mattress_uniformMass->setTotalMass(10);

	TetrahedronFEMForceField3::SPtr mattress_TetrahedronFEMForceField = addNew<TetrahedronFEMForceField3>(mattress);
	mattress_TetrahedronFEMForceField->setMethod(HexahedronFEMForceField3::SMALL);
	mattress_TetrahedronFEMForceField->setPoissonRatio(0.495);
	mattress_TetrahedronFEMForceField->setYoungModulus(3000);
	mattress_TetrahedronFEMForceField->_computeVonMisesStress.setValue(1);
	mattress_TetrahedronFEMForceField->_updateStiffness.setValue(true);
	mattress_TetrahedronFEMForceField->_showVonMisesStressPerNode.setValue(false);
	mattress_TetrahedronFEMForceField->_showStressColorMap.setValue("green  0   0.5 0   1  0.5  1   0   1");

	GenericConstraintCorrection::SPtr mattress_GenericConstraintCorrection = addNew<GenericConstraintCorrection>(mattress);

		#pragma region mattress_surface

		Node::SPtr mattress_surface = mattress->createChild("surface");

		RegularGridTopology::SPtr mattress_surface_RegularGridTopology = addNew<RegularGridTopology>(mattress_surface);
		mattress_surface_RegularGridTopology->setSize(3, 3, 3);
		mattress_surface_RegularGridTopology->setPos(-1, 1, 2, 2.2, -2, 2);

		MechanicalObject3::SPtr mattress_surface_MechanicalObject = addNew<MechanicalObject3>(mattress_surface);
		mattress_surface_MechanicalObject->setTranslation(0,0,0);
		mattress_surface_MechanicalObject->setRotation(0,0,0);
		mattress_surface_MechanicalObject->setScale(1,1,1);

		TriangleCollisionModel3::SPtr mattress_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(mattress_surface);
		LineCollisionModel3::SPtr mattress_surface_LineCollisionModel = addNew<LineCollisionModel3>(mattress_surface);
		PointCollisionModel3::SPtr mattress_surface_PointCollisionModel = addNew<PointCollisionModel3>(mattress_surface);

		SubsetMapping3::SPtr mattress_surface_SubsetMapping = addNew<SubsetMapping3>(mattress_surface);
		mattress_surface_SubsetMapping->setPathInputObject("@../.");
		mattress_surface_SubsetMapping->setPathOutputObject("@.");

		#pragma endregion	

	#pragma endregion
	

	#pragma region bed

	Node::SPtr bed = root->createChild("bed");

	EulerImplicitSolver::SPtr bed_EulerImplicitSolver = addNew<EulerImplicitSolver>(bed);
	SparseLDLSolverGraph::SPtr bed_SparseLDLSolverGraph = addNew<SparseLDLSolverGraph>(bed);

	RegularGridTopology::SPtr bed_RegularGridTopology = addNew<RegularGridTopology>(bed);
	bed_RegularGridTopology->setSize(3, 3, 3);
	bed_RegularGridTopology->setPos(-1, 1, 0, 1, -2, 2);

	MechanicalObject3::SPtr bed_MechanicalObject = addNew<MechanicalObject3>(bed);
	bed_MechanicalObject->setTranslation(0, 0.2, 0);
	bed_MechanicalObject->setRotation(0, 0, 0);
	bed_MechanicalObject->setScale(1, 1, 1);

	UniformMass3::SPtr bed_uniformMass = addNew<UniformMass3>(bed);
	bed_uniformMass->setTotalMass(50);

	TetrahedronFEMForceField3::SPtr bed_TetrahedronFEMForceField = addNew<TetrahedronFEMForceField3>(bed);
	bed_TetrahedronFEMForceField->setMethod(HexahedronFEMForceField3::SMALL);
	bed_TetrahedronFEMForceField->setPoissonRatio(0.495);
	bed_TetrahedronFEMForceField->setYoungModulus(3000);
	bed_TetrahedronFEMForceField->_computeVonMisesStress.setValue(1);
	bed_TetrahedronFEMForceField->_updateStiffness.setValue(true);
	bed_TetrahedronFEMForceField->_showVonMisesStressPerNode.setValue(false);
	bed_TetrahedronFEMForceField->_showStressColorMap.setValue("green  0   0.5 0   1  0.5  1   0   1");

	GenericConstraintCorrection::SPtr bed_GenericConstraintCorrection = addNew<GenericConstraintCorrection>(bed);

		#pragma region bed_surface

		Node::SPtr bed_surface = bed->createChild("surface");

		RegularGridTopology::SPtr bed_surface_RegularGridTopology = addNew<RegularGridTopology>(bed_surface);
		bed_surface_RegularGridTopology->setSize(3, 3, 3);
		bed_surface_RegularGridTopology->setPos(-1, 1, 0, 1, -2, 2);

		MechanicalObject3::SPtr bed_surface_MechanicalObject = addNew<MechanicalObject3>(bed_surface);
		bed_surface_MechanicalObject->setTranslation(0, 0.2, 0);
		bed_surface_MechanicalObject->setRotation(0, 0, 0);
		bed_surface_MechanicalObject->setScale(1, 1, 1);

		TriangleCollisionModel3::SPtr bed_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(bed_surface);
		LineCollisionModel3::SPtr bed_surface_LineCollisionModel = addNew<LineCollisionModel3>(bed_surface);
		PointCollisionModel3::SPtr bed_surface_PointCollisionModel = addNew<PointCollisionModel3>(bed_surface);

		SubsetMapping3::SPtr bed_surface_SubsetMapping = addNew<SubsetMapping3>(bed_surface);
		bed_surface_SubsetMapping->setPathInputObject("@../.");
		bed_surface_SubsetMapping->setPathOutputObject("@.");

		#pragma endregion

	#pragma endregion


	#pragma region floor

	Node::SPtr floor = root->createChild("floor");

	RegularGridTopology::SPtr floor_RegularGridTopology = addNew<RegularGridTopology>(floor);
	floor_RegularGridTopology->setSize(2, 2, 2);
	floor_RegularGridTopology->setPos(-10, 10, -0.5, 0, -10, 10);

	MechanicalObject3::SPtr floor_MechanicalObject = addNew<MechanicalObject3>(floor);
	floor_MechanicalObject->setTranslation(0, 0, 0);
	floor_MechanicalObject->setRotation(0, 0, 0);
	floor_MechanicalObject->setScale(1, 1, 1);

	OglModel::SPtr floor_oglModel = addNew<OglModel>(floor);
	floor_oglModel->loadTexture(sofa::helper::system::DataRepository.getFile(project_path + "/mesh/floor.bmp"));

	TriangleCollisionModel3::SPtr floor_TriangleCollisionModel = addNew<TriangleCollisionModel3>(floor);
	LineCollisionModel3::SPtr floor_LineCollisionModel = addNew<LineCollisionModel3>(floor);
	PointCollisionModel3::SPtr floor_PointCollisionModel = addNew<PointCollisionModel3>(floor);

	#pragma endregion 

	return root;
}

/// Create musculoskeletic system
simulation::Node::SPtr createSceneCUDA(std::string project_path)
{
#pragma region scene

	std::cout << "Project path: " << project_path << std::endl;
	std::string mesh_path = "mesh";// +project_path[project_path.length() - 1];
	mesh_path += project_path[project_path.length() - 1];
	std::cout << "Mesh path: " << project_path + mesh_path << std::endl;

	Node::SPtr root = simulation::getSimulation()->createNewGraph("root");
	root->setGravity(sofa::defaulttype::Vector3(0, -9.81, 0));
	root->setAnimate(false);
	root->setDt(0.02);

	VisualStyle::SPtr visualStyle = addNew<VisualStyle>(root);
	visualStyle->setName("VisualStyle");
	DisplayFlags& flags = *visualStyle->displayFlags.beginEdit();
	flags.setShowVisualModels(true);
	flags.setShowBehaviorModels(true);
	flags.setShowForceFields(true);
	flags.setShowInteractionForceFields(false);
	flags.setShowCollisionModels(true);
	flags.setShowBoundingCollisionModels(true);
	flags.setShowMappings(false);
	flags.setShowMechanicalMappings(false);
	flags.setShowWireFrame(false);
	flags.setShowNormals(false);
	visualStyle->displayFlags.endEdit();

	FreeMotionAnimationLoop::SPtr freeMotionAnimationLoop = New<FreeMotionAnimationLoop>(root.get());
	freeMotionAnimationLoop->setName("FreeMotionAnimationLoop");
	root->addObject(freeMotionAnimationLoop);

	GenericConstraintSolver::SPtr genericConstraintSolver = addNew<GenericConstraintSolver>(root);
	genericConstraintSolver->setName("GenericConstraintSolver");
	genericConstraintSolver->maxIt.setValue(200);
	genericConstraintSolver->tolerance.setValue(1.0e-3);

	DefaultPipeline::SPtr defaultPipeline = addNew<DefaultPipeline>(root);
	defaultPipeline->setName("DefaultPipeline");

	BruteForceDetection::SPtr bruteForceDetection = addNew<BruteForceDetection>(root);
	bruteForceDetection->setName("BruteForceDetection");

	component::collision::MinProximityIntersection::SPtr minProximityIntersection = addNew<component::collision::MinProximityIntersection>(root);
	minProximityIntersection->setName("MinProximityIntersection");
	minProximityIntersection->setAlarmDistance(1.2);
	minProximityIntersection->setContactDistance(0.1);

	DefaultContactManager::SPtr defaultContactManager = addNew<DefaultContactManager>(root);
	defaultContactManager->setName("DefaultContactManager");
	defaultContactManager->setDefaultResponseType("FrictionContact");

#pragma endregion

#pragma region tp

	Node::SPtr tp = root->createChild("tp");

	EulerImplicitSolver::SPtr tp_EulerImplicitSolver = addNew<EulerImplicitSolver>(tp);
	SparseLDLSolverGraph::SPtr tp_SparseLDLSolverGraph = addNew<SparseLDLSolverGraph>(tp);

	MeshGmshLoader::SPtr tp_MeshGmshLoader = addNew<MeshGmshLoader>(tp);
	tp_MeshGmshLoader->setName("loader");
	tp_MeshGmshLoader->setFilename(project_path + mesh_path + "tp.msh");
	tp_MeshGmshLoader->load();

	TetrahedronSetTopologyContainer::SPtr tp_TetrahedronSetTopologyContainer = addNew<TetrahedronSetTopologyContainer>(tp);
	tp_TetrahedronSetTopologyContainer->setSrc("@loader");

	MechanicalObject3::SPtr tp_MechanicalObject = addNew<MechanicalObject3>(tp);
	tp_MechanicalObject->setTranslation(-5, 5, 0);
	tp_MechanicalObject->setRotation(0, 90, -90);
	tp_MechanicalObject->setScale(10, 10, 10);

	TetrahedronSetGeometryAlgorithms3::SPtr tp_TetrahedronSetGeometryAlgorithms = addNew<TetrahedronSetGeometryAlgorithms3>(tp);

	UniformMass3::SPtr tp_uniformMass = addNew<UniformMass3>(tp);
	tp_uniformMass->setTotalMass(2);

	TetrahedronFEMForceField3::SPtr tp_TetrahedronFEMForceField = addNew<TetrahedronFEMForceField3>(tp);
	tp_TetrahedronFEMForceField->setMethod(HexahedronFEMForceField3::SMALL);
	tp_TetrahedronFEMForceField->setPoissonRatio(0.495);
	tp_TetrahedronFEMForceField->setYoungModulus(3000);
	tp_TetrahedronFEMForceField->_computeVonMisesStress.setValue(1);
	tp_TetrahedronFEMForceField->_updateStiffness.setValue(true);
	tp_TetrahedronFEMForceField->_showVonMisesStressPerNode.setValue(false);
	tp_TetrahedronFEMForceField->_showStressColorMap.setValue("green  0   0.5 0   1  0.5  1   0   1");

	GenericConstraintCorrection::SPtr tp_GenericConstraintCorrection = addNew<GenericConstraintCorrection>(tp);

	FixedConstraint3::SPtr tp_FixedConstraints = addNew<FixedConstraint3>(tp);
	tp_FixedConstraints->addConstraint(3);
	tp_FixedConstraints->addConstraint(39);
	tp_FixedConstraints->addConstraint(64);

#pragma region tp_surface

	Node::SPtr tp_surface = tp->createChild("surface");
	/*
	MeshObjLoader::SPtr tp_surface_MeshObjLoader = addNew<MeshObjLoader>(tp_surface);
	tp_surface_MeshObjLoader->setName("loader");
	tp_surface_MeshObjLoader->setFilename(project_path + mesh_path + "tp.obj");
	tp_surface_MeshObjLoader->load();

	MeshTopology::SPtr tp_surface_MeshTopology = addNew<MeshTopology>(tp_surface);
	tp_surface_MeshTopology->setSrc("@loader");

	MechanicalObject3::SPtr tp_surface_MechanicalObject = addNew<MechanicalObject3>(tp_surface);
	tp_surface_MechanicalObject->setTranslation(-5, 5, 0);
	tp_surface_MechanicalObject->setRotation(0, 90, -90);
	tp_surface_MechanicalObject->setScale(10, 10, 10);

	BarycentricMapping3::SPtr tp_surface_BarycentricMapping = addNew<BarycentricMapping3>(tp_surface);
	tp_surface_BarycentricMapping->setPathInputObject("@..");
	tp_surface_BarycentricMapping->setPathOutputObject("@.");

	TriangleCollisionModel3::SPtr tp_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(tp_surface);
	LineCollisionModel3::SPtr tp_surface_LineCollisionModel = addNew<LineCollisionModel3>(tp_surface);
	PointCollisionModel3::SPtr tp_surface_PointCollisionModel = addNew<PointCollisionModel3>(tp_surface);
	*/

	HexahedronSetTopologyContainer::SPtr tp_surface_HexahedronSetTopologyContainer = addNew<HexahedronSetTopologyContainer>(tp_surface);
	tp_surface_HexahedronSetTopologyContainer->setName("Container1");
	tp_surface_HexahedronSetTopologyContainer->setSrc("@../loader");
	HexahedronSetTopologyModifier::SPtr tp_surface_HexahedronSetTopologyModifier = addNew<HexahedronSetTopologyModifier>(tp_surface);
	HexahedronSetTopologyAlgorithms3::SPtr tp_surface_HexahedronSetTopologyAlgorithms = addNew<HexahedronSetTopologyAlgorithms3>(tp_surface);
	HexahedronSetGeometryAlgorithms3::SPtr tp_surface_HexahedronSetGeometryAlgorithms = addNew<HexahedronSetGeometryAlgorithms3>(tp_surface);

	QuadSetTopologyContainer::SPtr tp_surface_QuadSetTopologyContainer = addNew<QuadSetTopologyContainer>(tp_surface);
	tp_surface_QuadSetTopologyContainer->setName("Container2");
	QuadSetTopologyModifier::SPtr tp_surface_QuadSetTopologyModifier = addNew<QuadSetTopologyModifier>(tp_surface);
	QuadSetTopologyAlgorithms3::SPtr tp_surface_QuadSetTopologyAlgorithms = addNew<QuadSetTopologyAlgorithms3>(tp_surface);
	QuadSetGeometryAlgorithms3::SPtr tp_surface_QuadSetGeometryAlgorithms = addNew<QuadSetGeometryAlgorithms3>(tp_surface);

	Hexa2QuadTopologicalMapping::SPtr tp_surface_Hexa2QuadTopologicalMapping = addNew<Hexa2QuadTopologicalMapping>(tp_surface);
	tp_surface_Hexa2QuadTopologicalMapping->setPathInputObject("@Container1");
	tp_surface_Hexa2QuadTopologicalMapping->setPathOutputObject("@Container2");

	TriangleCollisionModel3::SPtr tp_surface_TriangleCollisionModel = addNew<TriangleCollisionModel3>(tp_surface);
	LineCollisionModel3::SPtr tp_surface_LineCollisionModel = addNew<LineCollisionModel3>(tp_surface);

#pragma endregion

#pragma endregion

#pragma region floor

	Node::SPtr floor = root->createChild("floor");

	RegularGridTopology::SPtr floor_RegularGridTopology = addNew<RegularGridTopology>(floor);
	floor_RegularGridTopology->setSize(2, 2, 2);
	floor_RegularGridTopology->setPos(-10, 10, -0.5, 0, -10, 10);

	MechanicalObject3::SPtr floor_MechanicalObject = addNew<MechanicalObject3>(floor);
	floor_MechanicalObject->setTranslation(0, 0, 0);
	floor_MechanicalObject->setRotation(0, 0, 0);
	floor_MechanicalObject->setScale(1, 1, 1);

	OglModel::SPtr floor_oglModel = addNew<OglModel>(floor);
	floor_oglModel->loadTexture(sofa::helper::system::DataRepository.getFile(project_path + "/mesh/floor.bmp"));

	TriangleCollisionModel3::SPtr floor_TriangleCollisionModel = addNew<TriangleCollisionModel3>(floor);
	LineCollisionModel3::SPtr floor_LineCollisionModel = addNew<LineCollisionModel3>(floor);
	PointCollisionModel3::SPtr floor_PointCollisionModel = addNew<PointCollisionModel3>(floor);

#pragma endregion 

	return root;
}

void addGUIParameters(ArgumentParser* argumentParser)
{
	GUIManager::RegisterParameters(argumentParser);
}

// ---------------------------------------------------------------------
// ---
// ---------------------------------------------------------------------

// Get project path
std::string project_path;

int main(int argc, char** argv)
{
	// Add resources dir to GuiDataRepository
	const std::string runSofaIniFilePath = Utils::getSofaPathTo("/etc/runSofa.ini");

	std::map<std::string, std::string> iniFileValues = Utils::readBasicIniFile(runSofaIniFilePath);
	if (iniFileValues.find("RESOURCES_DIR") != iniFileValues.end())
	{
		std::string dir = iniFileValues["RESOURCES_DIR"];
		dir = SetDirectory::GetRelativeFromProcess(dir.c_str());
		if (FileSystem::isDirectory(dir))
		{
			sofa::gui::GuiDataRepository.addFirstPath(dir);
		}

		try
		{
			// Get project path
			project_path = dir;

			if (project_path.find_last_of("/") > 0)
			{
				std::size_t found = project_path.find_last_of("/");
				project_path = project_path.substr(0, found + 1);
			}
			else if (project_path.find_last_of("\\") > 0)
			{
				std::size_t found = project_path.find_last_of("\\");
				project_path = project_path.substr(0, found + 1);
			}

			std::cout << "Project path: " << project_path << std::endl;
		}
		catch (...) {
			std::cout << "Problem parsing the project path: " << dir << std::endl;
		}
	}
	
	// Add plugins dir to PluginRepository
	if (FileSystem::isDirectory(Utils::getSofaPathPrefix() + "/plugins"))
	{
		PluginRepository.addFirstPath(Utils::getSofaPathPrefix() + "/plugins");
	}

	sofa::helper::BackTrace::autodump();

#ifdef WIN32
	{
		HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
		COORD s;
		s.X = 160; s.Y = 10000;
		SetConsoleScreenBufferSize(hStdout, s);
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		if (GetConsoleScreenBufferInfo(hStdout, &csbi))
		{
			SMALL_RECT winfo;
			winfo = csbi.srWindow;
			//winfo.Top = 0;
			winfo.Left = 0;
			//winfo.Bottom = csbi.dwSize.Y-1;
			winfo.Right = csbi.dwMaximumWindowSize.X - 1;
			SetConsoleWindowInfo(hStdout, TRUE, &winfo);
		}

	}
#endif

	sofa::gui::initMain();

	string fileName;
	bool        startAnim = false;
	bool        showHelp = false;
	bool        printFactory = false;
	bool        loadRecent = false;
	bool        temporaryFile = false;
	bool        testMode = false;
	bool        noAutoloadPlugins = false;
	bool        noSceneCheck = false;
	unsigned int nbMSSASamples = 1;
	bool computationTimeAtBegin = false;
	unsigned int computationTimeSampling = 0; ///< Frequency of display of the computation time statistics, in number of animation steps. 0 means never.
	string    computationTimeOutputType = "stdout";

	string gui = "";
	string verif = "";

#if defined(SOFA_HAVE_DAG)
	string simulationType = "dag";
#else
	string simulationType = "tree";
#endif

	std::vector<string> plugins;
	std::vector<string> files;

	string colorsStatus = "unset";
	string messageHandler = "auto";
	bool enableInteraction = false;
	int width = 800;
	int height = 600;

	string gui_help = "choose the UI (";
	gui_help += GUIManager::ListSupportedGUI('|');
	gui_help += ")";

	ArgumentParser* argParser = new ArgumentParser(argc, argv);
	argParser->addArgument(
		boost::program_options::value<bool>(&showHelp)
		->default_value(false)
		->implicit_value(true),
		"help,h",
		"Display this help message"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&startAnim)
		->default_value(false)
		->implicit_value(true),
		"start,a",
		"start the animation loop"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&computationTimeAtBegin)
		->default_value(false)
		->implicit_value(true),
		"computationTimeAtBegin,b",
		"Output computation time statistics of the init (at the begin of the simulation)"
	);
	argParser->addArgument(
		boost::program_options::value<unsigned int>(&computationTimeSampling)
		->default_value(0),
		"computationTimeSampling",
		"Frequency of display of the computation time statistics, in number of animation steps. 0 means never."
	);
	argParser->addArgument(
		boost::program_options::value<std::string>(&computationTimeOutputType)
		->default_value("stdout"),
		"computationTimeOutputType,o",
		"Output type for the computation time statistics: either stdout, json or ljson"
	);
	argParser->addArgument(
		boost::program_options::value<std::string>(&gui)->default_value(""),
		"gui,g",
		gui_help.c_str()
	);
	argParser->addArgument(
		boost::program_options::value<std::vector<std::string>>(&plugins),
		"load,l",
		"load given plugins"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&noAutoloadPlugins)
		->default_value(false)
		->implicit_value(true),
		"noautoload",
		"disable plugins autoloading"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&noSceneCheck)
		->default_value(false)
		->implicit_value(true),
		"noscenecheck",
		"disable scene checking for each scene loading"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&printFactory)
		->default_value(false)
		->implicit_value(true),
		"factory,p",
		"print factory logs"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&loadRecent)
		->default_value(false)->implicit_value(true),
		"recent,r",
		"load most recently opened file"
	);
	argParser->addArgument(
		boost::program_options::value<std::string>(&simulationType),
		"simu,s", "select the type of simulation (bgl, dag, tree)"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&temporaryFile)
		->default_value(false)->implicit_value(true),
		"tmp",
		"the loaded scene won't appear in history of opened files"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&testMode)
		->default_value(false)->implicit_value(true),
		"test",
		"select test mode with xml output after N iteration"
	);
	argParser->addArgument(
		boost::program_options::value<std::string>(&verif)
		->default_value(""),
		"verification,v",
		"load verification data for the scene"
	);
	argParser->addArgument(
		boost::program_options::value<std::string>(&colorsStatus)
		->default_value("unset", "auto")
		->implicit_value("yes"),
		"colors,c",
		"use colors on stdout and stderr (yes, no, auto)"
	);
	argParser->addArgument(
		boost::program_options::value<std::string>(&messageHandler)
		->default_value("auto"),
		"formatting,f",
		"select the message formatting to use (auto, clang, sofa, rich, test)"
	);
	argParser->addArgument(
		boost::program_options::value<bool>(&enableInteraction)
		->default_value(false)
		->implicit_value(true),
		"interactive,i",
		"enable interactive mode for the GUI which includes idle and mouse events (EXPERIMENTAL)"
	);
	argParser->addArgument(
		boost::program_options::value<std::vector<std::string> >()
		->multitoken(),
		"argv",
		"forward extra args to the python interpreter"
	);

	// example of an option using lambda function which ensure the value passed is > 0
	argParser->addArgument(
		boost::program_options::value<unsigned int>(&nbMSSASamples)
		->default_value(1)
		->notifier([](unsigned int value) {
		if (value < 1) {
			msg_error("runSofa") << "msaa sample cannot be lower than 1";
			exit(EXIT_FAILURE);
		}
	}),
		"msaa,m",
		"number of samples for MSAA (Multi Sampling Anti Aliasing ; value < 2 means disabled"
		);

	addGUIParameters(argParser);
	argParser->parse();
	files = argParser->getInputFileList();

	if (showHelp)
	{
		argParser->showHelp();
		exit(EXIT_SUCCESS);
	}

	// Note that initializations must be done after ArgumentParser that can exit the application (without cleanup)
	// even if everything is ok e.g. asking for help
	sofa::simulation::tree::init();
#ifdef SOFA_HAVE_DAG
	sofa::simulation::graph::init();
#endif
	sofa::component::initSofaBase();
	sofa::component::initSofaCommon();
	sofa::component::initSofaGeneral();
	sofa::component::initSofaMisc();

#ifdef SOFA_HAVE_DAG
	if (simulationType == "tree")
		sofa::simulation::setSimulation(new TreeSimulation());
	else
		sofa::simulation::setSimulation(new DAGSimulation());
#else //SOFA_HAVE_DAG
	sofa::simulation::setSimulation(new TreeSimulation());
#endif

	if (colorsStatus == "unset") {
		// If the parameter is unset, check the environment variable
		const char * colorStatusEnvironment = std::getenv("SOFA_COLOR_TERMINAL");
		if (colorStatusEnvironment != nullptr) {
			const std::string status(colorStatusEnvironment);
			if (status == "yes" || status == "on" || status == "always")
				sofa::helper::console::setStatus(sofa::helper::console::Status::On);
			else if (status == "no" || status == "off" || status == "never")
				sofa::helper::console::setStatus(sofa::helper::console::Status::Off);
			else
				sofa::helper::console::setStatus(sofa::helper::console::Status::Auto);
		}
	}
	else if (colorsStatus == "auto")
		sofa::helper::console::setStatus(sofa::helper::console::Status::Auto);
	else if (colorsStatus == "yes")
		sofa::helper::console::setStatus(sofa::helper::console::Status::On);
	else if (colorsStatus == "no")
		sofa::helper::console::setStatus(sofa::helper::console::Status::Off);

	//TODO(dmarchal): Use smart pointer there to avoid memory leaks !!
	if (messageHandler == "auto")
	{
		MessageDispatcher::clearHandlers();
		MessageDispatcher::addHandler(new ConsoleMessageHandler());
	}
	else if (messageHandler == "clang")
	{
		MessageDispatcher::clearHandlers();
		MessageDispatcher::addHandler(new ClangMessageHandler());
	}
	else if (messageHandler == "sofa")
	{
		MessageDispatcher::clearHandlers();
		MessageDispatcher::addHandler(new ConsoleMessageHandler());
	}
	else if (messageHandler == "rich")
	{
		MessageDispatcher::clearHandlers();
		MessageDispatcher::addHandler(new ConsoleMessageHandler(&RichConsoleStyleMessageFormatter::getInstance()));
	}
	else if (messageHandler == "test") {
		MessageDispatcher::addHandler(new ExceptionMessageHandler());
	}
	else {
		msg_warning("") << "Invalid argument '" << messageHandler << "' for '--formatting'";
	}
	MessageDispatcher::addHandler(&MainPerComponentLoggingMessageHandler::getInstance());

	// Output FileRepositories
	msg_info("runSofa") << "PluginRepository paths = " << PluginRepository.getPathsJoined();
	msg_info("runSofa") << "DataRepository paths = " << DataRepository.getPathsJoined();
	msg_info("runSofa") << "GuiDataRepository paths = " << GuiDataRepository.getPathsJoined();

	// Initialise paths
	BaseGUI::setConfigDirectoryPath(Utils::getSofaPathPrefix() + "/config", true);
	BaseGUI::setScreenshotDirectoryPath(Utils::getSofaPathPrefix() + "/screenshots", true);

	if (!files.empty())
		fileName = files[0];

	for (unsigned int i = 0; i < plugins.size(); i++)
		PluginManager::getInstance().loadPlugin(plugins[i]);

	std::string configPluginPath = sofa_tostring(CONFIG_PLUGIN_FILENAME);
	std::string defaultConfigPluginPath = sofa_tostring(DEFAULT_CONFIG_PLUGIN_FILENAME);

	if (!noAutoloadPlugins)
	{
		if (PluginRepository.findFile(configPluginPath, "", nullptr))
		{
			msg_info("runSofa") << "Loading automatically plugin list in " << configPluginPath;
			PluginManager::getInstance().readFromIniFile(configPluginPath);
		}
		else if (PluginRepository.findFile(defaultConfigPluginPath, "", nullptr))
		{
			msg_info("runSofa") << "Loading automatically plugin list in " << defaultConfigPluginPath;
			PluginManager::getInstance().readFromIniFile(defaultConfigPluginPath);
		}
		else
			msg_info("runSofa") << "No plugin list found. No plugin will be automatically loaded.";
	}
	else
		msg_info("runSofa") << "Automatic plugin loading disabled.";

	PluginManager::getInstance().init();

	if (int err = GUIManager::Init(argv[0], gui.c_str()))
		return err;

	if (fileName.empty())
	{
		if (loadRecent) // try to reload the latest scene
		{
			string scenes = BaseGUI::getConfigDirectoryPath() + "/runSofa.ini";
			std::ifstream mrulist(scenes.c_str());
			std::getline(mrulist, fileName);
			mrulist.close();
		}
		else
			fileName = "Demos/aescape.scn";

		//fileName = DataRepository.getFile(fileName);
		fileName = project_path + "scenes/VonMisesStress_CUDA.scn";
	}

	if (int err = GUIManager::createGUI(nullptr))
		return err;

	//To set a specific resolution for the viewer, use the component ViewerSetting in you scene graph
	GUIManager::SetDimension(width, height);

	// Create and register the SceneCheckerListener before scene loading
	if (!noSceneCheck)
	{
		sofa::simulation::SceneLoader::addListener(SceneCheckerListener::getInstance());
	}

	const std::vector<std::string> sceneArgs = sofa::helper::ArgumentParser::extra_args();

	sofa::simulation::Node::SPtr groot;
	//groot = sofa::simulation::getSimulation()->load(fileName, false, sceneArgs);
	if (!groot) {

		//=================================================
		//groot = createSceneCUDA(project_path);
		groot = createScene(project_path);
		//=================================================
	}

	if (!verif.empty())
	{
		runSofa::Validation::execute(verif, fileName, groot.get());
	}

	if (computationTimeAtBegin)
	{
		sofa::helper::AdvancedTimer::setEnabled("Init", true);
		sofa::helper::AdvancedTimer::setInterval("Init", 1);
		sofa::helper::AdvancedTimer::setOutputType("Init", computationTimeOutputType);
		sofa::helper::AdvancedTimer::begin("Init");
	}

	sofa::simulation::getSimulation()->init(groot.get());

	if (computationTimeAtBegin)
	{
		msg_info("") << sofa::helper::AdvancedTimer::end("Init", groot.get());
	}

	//=======================================
	//Apply Options

	// start anim option
	if (startAnim)
		groot->setAnimate(true);

	// set scene and animation root to the gui
	GUIManager::SetScene(groot, fileName.c_str(), temporaryFile);

	if (printFactory)
	{
		msg_info("") << "////////// FACTORY //////////";
		sofa::helper::printFactoryLog();
		msg_info("") << "//////// END FACTORY ////////";
	}

	if (computationTimeSampling > 0)
	{
		sofa::helper::AdvancedTimer::setEnabled("Animate", true);
		sofa::helper::AdvancedTimer::setInterval("Animate", computationTimeSampling);
		sofa::helper::AdvancedTimer::setOutputType("Animate", computationTimeOutputType);
	}

	//=======================================
	// Run the main loop
	if (int err = GUIManager::MainLoop(groot, fileName.c_str()))
		return err;

	groot = dynamic_cast<Node*>(GUIManager::CurrentSimulation());

	if (testMode)
	{
		string xmlname = fileName.substr(0, fileName.length() - 4) + "-scene.scn";
		msg_info("") << "Exporting to XML " << xmlname;
		sofa::simulation::getSimulation()->exportXML(groot.get(), xmlname.c_str());
	}

	if (groot != nullptr)
		sofa::simulation::getSimulation()->unload(groot);


	GUIManager::closeGUI();

	sofa::simulation::common::cleanup();
	sofa::simulation::tree::cleanup();
#ifdef SOFA_HAVE_DAG
	sofa::simulation::graph::cleanup();
#endif
	return 0;
}



/*
#include <sstream>
using std::ostringstream ;
#include <fstream>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <runSofaValidation.h>

#include <sofa/helper/ArgumentParser.h>
#include <SofaSimulationCommon/config.h>
#include <sofa/simulation/Node.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/simulation/config.h> // #defines SOFA_HAVE_DAG (or not)
#include <SofaSimulationCommon/init.h>
#ifdef SOFA_HAVE_DAG
#include <SofaSimulationGraph/init.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#endif
#include <SofaSimulationTree/init.h>
#include <SofaSimulationTree/TreeSimulation.h>
using sofa::simulation::Node;
#include <sofa/simulation/SceneLoaderFactory.h>
#include <SofaGraphComponent/SceneCheckerListener.h>
using sofa::simulation::scenechecking::SceneCheckerListener;

#include <SofaCommon/initSofaCommon.h>
#include <SofaBase/initSofaBase.h>
#include <SofaGeneral/initSofaGeneral.h>
#include <SofaMisc/initSofaMisc.h>

#include <SofaGeneralLoader/ReadState.h>
#include <sofa/helper/Factory.h>
#include <sofa/helper/cast.h>
#include <sofa/helper/BackTrace.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/FileSystem.h>
using sofa::helper::system::FileSystem;
#include <sofa/helper/system/SetDirectory.h>
#include <sofa/helper/Utils.h>
#include <sofa/gui/GUIManager.h>
using sofa::gui::GUIManager;

#include <sofa/gui/Main.h>
#include <sofa/gui/BatchGUI.h>  // For the default number of iterations
#include <sofa/helper/system/gl.h>

using sofa::core::ExecParams ;

#include <sofa/helper/system/console.h>
using sofa::helper::Utils;

using sofa::component::misc::ReadStateActivator;
using sofa::simulation::tree::TreeSimulation;
using sofa::simulation::graph::DAGSimulation;
using sofa::helper::system::SetDirectory;
using sofa::core::objectmodel::BaseNode ;
using sofa::gui::BatchGUI;
using sofa::gui::BaseGUI;

#include <sofa/helper/logging/Messaging.h>

#include <sofa/helper/logging/ConsoleMessageHandler.h>
using sofa::helper::logging::ConsoleMessageHandler ;

#include <sofa/core/logging/RichConsoleStyleMessageFormatter.h>
using  sofa::helper::logging::RichConsoleStyleMessageFormatter ;

#include <sofa/core/logging/PerComponentLoggingMessageHandler.h>
using  sofa::helper::logging::MainPerComponentLoggingMessageHandler ;

#include <sofa/helper/AdvancedTimer.h>

#include <sofa/gui/GuiDataRepository.h>
using sofa::gui::GuiDataRepository ;

using sofa::helper::system::DataRepository;
using sofa::helper::system::PluginRepository;
using sofa::helper::system::PluginManager;

#include <sofa/helper/logging/MessageDispatcher.h>
using sofa::helper::logging::MessageDispatcher ;

#include <sofa/helper/logging/ClangMessageHandler.h>
using sofa::helper::logging::ClangMessageHandler ;

#include <sofa/helper/logging/ExceptionMessageHandler.h>
using sofa::helper::logging::ExceptionMessageHandler;

#include <boost/program_options.hpp>



void addGUIParameters(ArgumentParser* argumentParser)
{
    GUIManager::RegisterParameters(argumentParser);
}


// ---------------------------------------------------------------------
// ---
// ---------------------------------------------------------------------
int main(int argc, char** argv)
{
    // Add resources dir to GuiDataRepository
    const std::string runSofaIniFilePath = Utils::getSofaPathTo("/etc/runSofa.ini");
    std::map<std::string, std::string> iniFileValues = Utils::readBasicIniFile(runSofaIniFilePath);
    if (iniFileValues.find("RESOURCES_DIR") != iniFileValues.end())
    {
        std::string dir = iniFileValues["RESOURCES_DIR"];
        dir = SetDirectory::GetRelativeFromProcess(dir.c_str());
        if(FileSystem::isDirectory(dir))
        {
            sofa::gui::GuiDataRepository.addFirstPath(dir);
        }
    }

    // Add plugins dir to PluginRepository
    if ( FileSystem::isDirectory(Utils::getSofaPathPrefix()+"/plugins") )
    {
        PluginRepository.addFirstPath(Utils::getSofaPathPrefix()+"/plugins");
    }

    sofa::helper::BackTrace::autodump();

#ifdef WIN32
    {
        HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
        COORD s;
        s.X = 160; s.Y = 10000;
        SetConsoleScreenBufferSize(hStdout, s);
        CONSOLE_SCREEN_BUFFER_INFO csbi;
        if (GetConsoleScreenBufferInfo(hStdout, &csbi))
        {
            SMALL_RECT winfo;
            winfo = csbi.srWindow;
            //winfo.Top = 0;
            winfo.Left = 0;
            //winfo.Bottom = csbi.dwSize.Y-1;
            winfo.Right = csbi.dwMaximumWindowSize.X-1;
            SetConsoleWindowInfo(hStdout, TRUE, &winfo);
        }

    }
#endif

    sofa::gui::initMain();

    string fileName ;
    bool        startAnim = false;
    bool        showHelp = false;
    bool        printFactory = false;
    bool        loadRecent = false;
    bool        temporaryFile = false;
    bool        testMode = false;
    bool        noAutoloadPlugins = false;
    bool        noSceneCheck = false;
    unsigned int nbMSSASamples = 1;
    bool computationTimeAtBegin = false;
    unsigned int computationTimeSampling=0; ///< Frequency of display of the computation time statistics, in number of animation steps. 0 means never.
    string    computationTimeOutputType="stdout";

    string gui = "";
    string verif = "";

#if defined(SOFA_HAVE_DAG)
    string simulationType = "dag";
#else
    string simulationType = "tree";
#endif

    vector<string> plugins;
    vector<string> files;

    string colorsStatus = "unset";
    string messageHandler = "auto";
    bool enableInteraction = false ;
    int width = 800;
    int height = 600;

    string gui_help = "choose the UI (";
    gui_help += GUIManager::ListSupportedGUI('|');
    gui_help += ")";

    ArgumentParser* argParser = new ArgumentParser(argc, argv);
    argParser->addArgument(
        boost::program_options::value<bool>(&showHelp)
        ->default_value(false)
        ->implicit_value(true),
        "help,h",
        "Display this help message"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&startAnim)
        ->default_value(false)
        ->implicit_value(true),
        "start,a",
        "start the animation loop"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&computationTimeAtBegin)
        ->default_value(false)
        ->implicit_value(true),
        "computationTimeAtBegin,b",
        "Output computation time statistics of the init (at the begin of the simulation)"
    );
    argParser->addArgument(
        boost::program_options::value<unsigned int>(&computationTimeSampling)
        ->default_value(0),
        "computationTimeSampling",
        "Frequency of display of the computation time statistics, in number of animation steps. 0 means never."
    );
    argParser->addArgument(
        boost::program_options::value<std::string>(&computationTimeOutputType)
        ->default_value("stdout"),
        "computationTimeOutputType,o",
        "Output type for the computation time statistics: either stdout, json or ljson"
    );
    argParser->addArgument(
        boost::program_options::value<std::string>(&gui)->default_value(""),
        "gui,g",
        gui_help.c_str()
    );
    argParser->addArgument(
        boost::program_options::value<std::vector<std::string>>(&plugins),
        "load,l",
        "load given plugins"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&noAutoloadPlugins)
        ->default_value(false)
        ->implicit_value(true),
        "noautoload",
        "disable plugins autoloading"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&noSceneCheck)
        ->default_value(false)
        ->implicit_value(true),
        "noscenecheck",
        "disable scene checking for each scene loading"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&printFactory)
        ->default_value(false)
        ->implicit_value(true),
        "factory,p",
        "print factory logs"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&loadRecent)
        ->default_value(false)->implicit_value(true),
        "recent,r",
        "load most recently opened file"
    );
    argParser->addArgument(
        boost::program_options::value<std::string>(&simulationType),
        "simu,s", "select the type of simulation (bgl, dag, tree)"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&temporaryFile)
        ->default_value(false)->implicit_value(true),
        "tmp",
        "the loaded scene won't appear in history of opened files"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&testMode)
        ->default_value(false)->implicit_value(true),
        "test",
        "select test mode with xml output after N iteration"
     );
    argParser->addArgument(
        boost::program_options::value<std::string>(&verif)
        ->default_value(""),
        "verification,v",
        "load verification data for the scene"
    );
    argParser->addArgument(
        boost::program_options::value<std::string>(&colorsStatus)
        ->default_value("unset", "auto")
        ->implicit_value("yes"),
        "colors,c",
        "use colors on stdout and stderr (yes, no, auto)"
    );
    argParser->addArgument(
        boost::program_options::value<std::string>(&messageHandler)
        ->default_value("auto"),
        "formatting,f",
        "select the message formatting to use (auto, clang, sofa, rich, test)"
    );
    argParser->addArgument(
        boost::program_options::value<bool>(&enableInteraction)
        ->default_value(false)
        ->implicit_value(true),
        "interactive,i",
        "enable interactive mode for the GUI which includes idle and mouse events (EXPERIMENTAL)"
    );
    argParser->addArgument(
        boost::program_options::value<std::vector<std::string> >()
        ->multitoken(),
        "argv",
        "forward extra args to the python interpreter"
    );

    // example of an option using lambda function which ensure the value passed is > 0
    argParser->addArgument(
        boost::program_options::value<unsigned int>(&nbMSSASamples)
        ->default_value(1)
        ->notifier([](unsigned int value) {
            if (value < 1) {
                msg_error("runSofa") << "msaa sample cannot be lower than 1";
                exit( EXIT_FAILURE );
            }
        }),
        "msaa,m",
        "number of samples for MSAA (Multi Sampling Anti Aliasing ; value < 2 means disabled"
    );

    addGUIParameters(argParser);
    argParser->parse();
    files = argParser->getInputFileList();

    if(showHelp)
    {
        argParser->showHelp();
        exit( EXIT_SUCCESS );
    }

    // Note that initializations must be done after ArgumentParser that can exit the application (without cleanup)
    // even if everything is ok e.g. asking for help
    sofa::simulation::tree::init();
#ifdef SOFA_HAVE_DAG
    sofa::simulation::graph::init();
#endif
    sofa::component::initSofaBase();
    sofa::component::initSofaCommon();
    sofa::component::initSofaGeneral();
    sofa::component::initSofaMisc();

#ifdef SOFA_HAVE_DAG
    if (simulationType == "tree")
        sofa::simulation::setSimulation(new TreeSimulation());
    else
        sofa::simulation::setSimulation(new DAGSimulation());
#else //SOFA_HAVE_DAG
    sofa::simulation::setSimulation(new TreeSimulation());
#endif

    if (colorsStatus == "unset") {
        // If the parameter is unset, check the environment variable
        const char * colorStatusEnvironment = std::getenv("SOFA_COLOR_TERMINAL");
        if (colorStatusEnvironment != nullptr) {
            const std::string status (colorStatusEnvironment);
            if (status == "yes" || status == "on" || status == "always")
                sofa::helper::console::setStatus(sofa::helper::console::Status::On);
            else if (status == "no" || status == "off" || status == "never")
                sofa::helper::console::setStatus(sofa::helper::console::Status::Off);
            else
                sofa::helper::console::setStatus(sofa::helper::console::Status::Auto);
        }
    } else if (colorsStatus == "auto")
        sofa::helper::console::setStatus(sofa::helper::console::Status::Auto);
    else if (colorsStatus == "yes")
        sofa::helper::console::setStatus(sofa::helper::console::Status::On);
    else if (colorsStatus == "no")
        sofa::helper::console::setStatus(sofa::helper::console::Status::Off);

    //TODO(dmarchal): Use smart pointer there to avoid memory leaks !!
    if (messageHandler == "auto" )
    {
        MessageDispatcher::clearHandlers() ;
        MessageDispatcher::addHandler( new ConsoleMessageHandler() ) ;
    }
    else if (messageHandler == "clang")
    {
        MessageDispatcher::clearHandlers() ;
        MessageDispatcher::addHandler( new ClangMessageHandler() ) ;
    }
    else if (messageHandler == "sofa")
    {
        MessageDispatcher::clearHandlers() ;
        MessageDispatcher::addHandler( new ConsoleMessageHandler() ) ;
    }
    else if (messageHandler == "rich")
    {
        MessageDispatcher::clearHandlers() ;
        MessageDispatcher::addHandler( new ConsoleMessageHandler(&RichConsoleStyleMessageFormatter::getInstance()) ) ;
    }
    else if (messageHandler == "test"){
        MessageDispatcher::addHandler( new ExceptionMessageHandler() ) ;
    }
    else{
        msg_warning("") << "Invalid argument '" << messageHandler << "' for '--formatting'";
    }
    MessageDispatcher::addHandler(&MainPerComponentLoggingMessageHandler::getInstance()) ;

    // Output FileRepositories
    msg_info("runSofa") << "PluginRepository paths = " << PluginRepository.getPathsJoined();
    msg_info("runSofa") << "DataRepository paths = " << DataRepository.getPathsJoined();
    msg_info("runSofa") << "GuiDataRepository paths = " << GuiDataRepository.getPathsJoined();

    // Initialise paths
    BaseGUI::setConfigDirectoryPath(Utils::getSofaPathPrefix() + "/config", true);
    BaseGUI::setScreenshotDirectoryPath(Utils::getSofaPathPrefix() + "/screenshots", true);

    if (!files.empty())
        fileName = files[0];

    for (unsigned int i=0; i<plugins.size(); i++)
        PluginManager::getInstance().loadPlugin(plugins[i]);

    std::string configPluginPath = sofa_tostring(CONFIG_PLUGIN_FILENAME);
    std::string defaultConfigPluginPath = sofa_tostring(DEFAULT_CONFIG_PLUGIN_FILENAME);

    if (!noAutoloadPlugins)
    {
        if (PluginRepository.findFile(configPluginPath, "", nullptr))
        {
            msg_info("runSofa") << "Loading automatically plugin list in " << configPluginPath;
            PluginManager::getInstance().readFromIniFile(configPluginPath);
        }
        else if (PluginRepository.findFile(defaultConfigPluginPath, "", nullptr))
        {
            msg_info("runSofa") << "Loading automatically plugin list in " << defaultConfigPluginPath;
            PluginManager::getInstance().readFromIniFile(defaultConfigPluginPath);
        }
        else
            msg_info("runSofa") << "No plugin list found. No plugin will be automatically loaded.";
    }
    else
        msg_info("runSofa") << "Automatic plugin loading disabled.";

    PluginManager::getInstance().init();

    if (int err = GUIManager::Init(argv[0],gui.c_str()))
        return err;

    if (fileName.empty())
    {
        if (loadRecent) // try to reload the latest scene
        {
            string scenes = BaseGUI::getConfigDirectoryPath() + "/runSofa.ini";
            std::ifstream mrulist(scenes.c_str());
            std::getline(mrulist,fileName);
            mrulist.close();
        }
        else
            fileName = "Demos/caduceus.scn";

        fileName = DataRepository.getFile(fileName);
    }


    if (int err=GUIManager::createGUI(nullptr))
        return err;

    //To set a specific resolution for the viewer, use the component ViewerSetting in you scene graph
    GUIManager::SetDimension(width, height);

    // Create and register the SceneCheckerListener before scene loading
    if(!noSceneCheck)
    {
        sofa::simulation::SceneLoader::addListener( SceneCheckerListener::getInstance() );
    }

    const std::vector<std::string> sceneArgs = sofa::helper::ArgumentParser::extra_args();
    Node::SPtr groot = sofa::simulation::getSimulation()->load(fileName, false, sceneArgs);
    if( !groot )
        groot = sofa::simulation::getSimulation()->createNewGraph("");

    if (!verif.empty())
    {
        runSofa::Validation::execute(verif, fileName, groot.get());
    }

    if( computationTimeAtBegin )
    {
        sofa::helper::AdvancedTimer::setEnabled("Init", true);
        sofa::helper::AdvancedTimer::setInterval("Init", 1);
        sofa::helper::AdvancedTimer::setOutputType("Init", computationTimeOutputType);
        sofa::helper::AdvancedTimer::begin("Init");
    }

    sofa::simulation::getSimulation()->init(groot.get());
    if( computationTimeAtBegin )
    {
        msg_info("") << sofa::helper::AdvancedTimer::end("Init", groot.get());
    }

    //=======================================
    //Apply Options

    // start anim option
    if (startAnim)
        groot->setAnimate(true);

    // set scene and animation root to the gui
    GUIManager::SetScene(groot, fileName.c_str(), temporaryFile);

    if (printFactory)
    {
        msg_info("") << "////////// FACTORY //////////" ;
        sofa::helper::printFactoryLog();
        msg_info("") << "//////// END FACTORY ////////" ;
    }

    if( computationTimeSampling>0 )
    {
        sofa::helper::AdvancedTimer::setEnabled("Animate", true);
        sofa::helper::AdvancedTimer::setInterval("Animate", computationTimeSampling);
        sofa::helper::AdvancedTimer::setOutputType("Animate", computationTimeOutputType);
    }

    //=======================================
    // Run the main loop
    if (int err = GUIManager::MainLoop(groot,fileName.c_str()))
        return err;
    groot = dynamic_cast<Node*>( GUIManager::CurrentSimulation() );

    if (testMode)
    {
        string xmlname = fileName.substr(0,fileName.length()-4)+"-scene.scn";
        msg_info("") << "Exporting to XML " << xmlname ;
        sofa::simulation::getSimulation()->exportXML(groot.get(), xmlname.c_str());
    }

    if (groot!=nullptr)
        sofa::simulation::getSimulation()->unload(groot);


    GUIManager::closeGUI();

    sofa::simulation::common::cleanup();
    sofa::simulation::tree::cleanup();
#ifdef SOFA_HAVE_DAG
    sofa::simulation::graph::cleanup();
#endif
    return 0;
}
*/