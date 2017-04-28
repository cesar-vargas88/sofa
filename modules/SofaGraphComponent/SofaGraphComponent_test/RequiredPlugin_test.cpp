#include <SofaTest/Sofa_test.h>
#include <SceneCreator/SceneCreator.h>

#include <SofaTest/TestMessageHandler.h>
using sofa::helper::logging::Message ;
using sofa::helper::logging::ExpectMessage ;
using sofa::helper::logging::MessageAsTestFailure;

#include <SofaSimulationCommon/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;
using sofa::simulation::Node ;

using sofa::core::ExecParams;

namespace sofa
{

struct RequiredPlugin_test : public Sofa_test<>
{
    void testNotExistingPlugin()
    {
        ExpectMessage raii(Message::Error);

        std::stringstream scene ;
        scene << "<?xml version='1.0'?>"
                 "<Node 	name='Root' gravity='0 -9.81 0' time='0' animate='0' >               \n"
                 "   <RequiredPlugin name=\"notExist\" pluginName=\"SofaNotExist\" />            \n"
                 "</Node>                                                                        \n" ;

        Node::SPtr root = SceneLoaderXML::loadFromMemory ("testscene",
                                                          scene.str().c_str(),
                                                          scene.str().size()) ;

        ASSERT_NE(root.get(), nullptr) ;
        root->init(ExecParams::defaultInstance()) ;
    }

    void testNoParameter()
    {
        ExpectMessage raii(Message::Error);

        std::stringstream scene ;
        scene << "<?xml version='1.0'?>"
                 "<Node 	name='Root' gravity='0 -9.81 0' time='0' animate='0' >               \n"
                 "   <RequiredPlugin />            \n"
                 "</Node>                                                                        \n" ;

        Node::SPtr root = SceneLoaderXML::loadFromMemory ("testscene",
                                                          scene.str().c_str(),
                                                          scene.str().size()) ;

        ASSERT_NE(root.get(), nullptr) ;
        root->init(ExecParams::defaultInstance()) ;
    }
};

TEST_F(RequiredPlugin_test, testNotExistingPlugin ) { testNotExistingPlugin(); }
TEST_F(RequiredPlugin_test, testNoParameter ) { testNoParameter(); }

}
