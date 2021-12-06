#include <GlfwApp.h>

#include <SceneGraph.h>
#include <Topology/MixSet.h>
#include <Peridynamics/Dolphin.h>
#include <Peridynamics/ElasticBody.h>
#include <Peridynamics/ElasticityModule.h>
#include <ParticleSystem/StaticBoundary.h>

// Internal OpenGL Renderer
#include <GLRenderEngine.h>
#include <GLPointVisualModule.h>
#include <GLSurfaceVisualModule.h>

using namespace dyno;

int main()
{
	SceneGraph& scene = SceneGraph::getInstance();

	// set scene
	std::shared_ptr<StaticBoundary<DataType3f>> root = scene.createNewScene<StaticBoundary<DataType3f>>();
	root->loadCube(Vec3f(0), Vec3f(1), 0.005f, true);
	// root->loadShpere(Vec3f(0.5, 0.7f, 0.5), 0.08f, 0.005f, false, true);

	// set dolphin
	std::shared_ptr<Dolphin<DataType3f>> dolphin = std::make_shared<Dolphin<DataType3f>>();
	root->addParticleSystem(dolphin); 

	dolphin->setMass(1.0f);
	dolphin->loadMixFile("../../data/dolphin/Dolphin");
	dolphin->scale(1.0f);
	dolphin->translate(Vec3f(0.5f, 0.1f, 0.5f));
	dolphin->setVisible(true);

	
	bool useVTK = true;
	RenderEngine* engine;

	engine = new GLRenderEngine;
	// set point 
	auto pointRenderer = std::make_shared<GLPointVisualModule>();
	pointRenderer->setColor(Vec3f(1, 0.2, 1));
	pointRenderer->setColorMapMode(GLPointVisualModule::PER_OBJECT_SHADER);
	dolphin->getSurfaceNode()->currentTopology()->connect(pointRenderer->inPointSet());
	dolphin->getSurfaceNode()->graphicsPipeline()->pushModule(pointRenderer);
	// dolphin->currentPoints()->connect(pointRenderer->inPointSet());
	// dolphin->currentVelocity()->connect(pointRenderer->inColor());

	// dolphin->graphicsPipeline()->pushModule(pointRenderer);

	// pointRenderer->setVisible(false);

	// set surface
	auto sRender = std::make_shared<GLSurfaceVisualModule>();
	sRender->setColor(Vec3f(1, 1, 0));
	dolphin->getSurfaceNode()->currentTopology()->connect(sRender->inTriangleSet());
	dolphin->getSurfaceNode()->graphicsPipeline()->pushModule(sRender);
	// sRender->setVisible(false);

	GlfwApp window;
	window.setRenderEngine(engine);
	window.createWindow(1024, 768);
	window.mainLoop();

	delete engine;
    return 0;
}