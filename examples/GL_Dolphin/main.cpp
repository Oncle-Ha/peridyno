#include <GlfwApp.h>

#include <SceneGraph.h>
#include <Topology/MixSet.h>
#include <Topology/JointTree.h>
#include <Peridynamics/Dolphin.h>
#include <Peridynamics/ElasticBody.h>
#include <Peridynamics/ElasticityModule.h>
#include <ParticleSystem/StaticBoundary.h>

#include <ofbx.h>

// Internal OpenGL Renderer
#include <GLRenderEngine.h>
#include <GLPointVisualModule.h>
#include <GLSurfaceVisualModule.h>

using namespace dyno;

// TODO: move to Topology
ofbx::IScene* g_scene = nullptr;
bool init(const char* filepath)
{

	FILE* fp = fopen(filepath, "rb");

	if (!fp) return false;

	fseek(fp, 0, SEEK_END);
	long file_size = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	auto* content = new ofbx::u8[file_size];
	fread(content, 1, file_size, fp);

	g_scene = ofbx::load((ofbx::u8*)content, file_size, (ofbx::u64)ofbx::LoadFlags::TRIANGULATE);

	delete[] content;
	fclose(fp);

	return true;
}


void copyVec(DataType3f::Coord &dest, ofbx::Vec3 src)
{
	dest[0] = src.x;
	dest[1] = src.x;
	dest[2] = src.y;
}

void getModelProperties(const ofbx::Object& object, std::shared_ptr<JointTree<DataType3f>> cur)
{
	cur->id = object.id;
	copyVec(cur->PreRotation, object.getPreRotation());
	copyVec(cur->LclTranslation, object.getLocalTranslation());
	copyVec(cur->LclRotation, object.getLocalRotation());
	copyVec(cur->LclScaling, object.getLocalScaling());
}


void getLimbNode(const ofbx::Object& object, std::shared_ptr<JointTree<DataType3f>> parent, std::shared_ptr<JointTree<DataType3f>> limbRoot)
{
	if (object.getType() != ofbx::Object::Type::LIMB_NODE && parent != limbRoot) return;

	std::shared_ptr<JointTree<DataType3f>> cur;

	if (object.getType() == ofbx::Object::Type::LIMB_NODE){

		cur = std::make_shared<JointTree<DataType3f>>();
		parent->children.push_back(cur);
		cur->parent = parent;
		getModelProperties(object, cur);
	}
	int i = 0;
	while (ofbx::Object* child = object.resolveObjectLink(i))
	{
		if (object.getType() == ofbx::Object::Type::LIMB_NODE) getLimbNode(*child, cur, limbRoot);
		else getLimbNode(*child, parent, limbRoot);
		++i;
	}	
}

void getLimbNodes(const ofbx::IScene& scene, std::shared_ptr<JointTree<DataType3f>> limbRoot)
{
	const ofbx::Object* root = scene.getRoot();
	if (root) getLimbNode(*root, limbRoot, limbRoot);
}

void loadFBX(const char* filepath, std::shared_ptr<JointTree<DataType3f>> limbRoot)
{
	init(filepath);
	getLimbNodes(*g_scene, limbRoot);
}


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
	dolphin->scale(0.2f); // 太大会导致粒子间距过大以至于邻居为空
	dolphin->translate(Vec3f(0.5f, 0.1f, 0.5f));
	dolphin->setVisible(true);

	
	bool useVTK = true;
	RenderEngine* engine;

	engine = new GLRenderEngine;
	// set point 
	/*
	//auto pointRenderer = std::make_shared<GLPointVisualModule>();
	//pointRenderer->setColor(Vec3f(1, 0.2, 1));
	//pointRenderer->setColorMapMode(GLPointVisualModule::PER_OBJECT_SHADER);

	//dolphin->getSurfaceNode()->currentTopology()->connect(pointRenderer->inPointSet());
	//dolphin->getSurfaceNode()->graphicsPipeline()->pushModule(pointRenderer);

	//dolphin->currentPoints()->connect(pointRenderer->inPointSet());
	//dolphin->currentVelocity()->connect(pointRenderer->inColor());

	//dolphin->graphicsPipeline()->pushModule(pointRenderer);

	// pointRenderer->setVisible(false);
	*/

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