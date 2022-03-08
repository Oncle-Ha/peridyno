#include <GlfwApp.h>

#include <SceneGraph.h>
#include <Topology/MixSet.h>
#include <Topology/JointTree.h>
#include <Topology/Cluster.h>
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
std::shared_ptr<Dolphin<DataType3f>> temp_Dolphin;

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


void copyVec(DataType3f::Coord &dest, ofbx::Vec3 src){dest = DataType3f::Coord(src.x, src.y, src.z);}
void copyVecR(DataType3f::Coord &dest, ofbx::Vec3 src){dest = DataType3f::Coord(src.x, -src.y, -src.z);}
void copyVecT(DataType3f::Coord &dest, ofbx::Vec3 src){dest = DataType3f::Coord(src.x, src.y, -src.z);}

void getModelProperties(const ofbx::Object& object, std::shared_ptr<JointTree<DataType3f>> cur)
{
	cur->id = object.id;

	copyVec(cur->PreRotation, object.getPreRotation());
	copyVecT(cur->LclTranslation, object.getLocalTranslation());
	copyVecR(cur->LclRotation, object.getLocalRotation());
	copyVec(cur->LclScaling, object.getLocalScaling());

	temp_Dolphin->m_jointMap.push_back(cur);
}


void getLimbNode(const ofbx::Object& object, std::shared_ptr<JointTree<DataType3f>> parent)
{
	if (object.getType() != ofbx::Object::Type::LIMB_NODE) return;

	std::shared_ptr<JointTree<DataType3f>> cur;

	if (object.getType() == ofbx::Object::Type::LIMB_NODE){

		cur = std::make_shared<JointTree<DataType3f>>();
		if(parent != nullptr) parent->children.push_back(cur);
		cur->parent = parent;
		getModelProperties(object, cur);
	}
	int i = 0;
	while (ofbx::Object* child = object.resolveObjectLink(i))
	{
		if (object.getType() == ofbx::Object::Type::LIMB_NODE) getLimbNode(*child, cur);
		else getLimbNode(*child, parent);
		++i;
	}	
}


// import Skin
// For one mesh
void getClusterProperties(const ofbx::Object& object)
{
	ofbx::Cluster* obj_cluster = (ofbx::Cluster*)&object;
	std::shared_ptr<Cluster<DataType3f>> temp_cluster = std::make_shared<Cluster<DataType3f>>(
		obj_cluster->getIndices(), obj_cluster->getIndicesCount(),
		obj_cluster->getWeights(), obj_cluster->getWeightsCount(),
		obj_cluster->getTransformMatrix().m, obj_cluster->getTransformLinkMatrix().m);
	int num = 0; int size = temp_Dolphin->m_jointMap.size();
	for(; num < size; ++num)
	if(obj_cluster->getLink()->id == temp_Dolphin->m_jointMap[num]->id) break;

	assert(num != size);

	temp_cluster->m_jointIndex = num;

	temp_Dolphin->m_clusters.push_back(temp_cluster);
}

void getCluster(const ofbx::Object& object)
{
	if (object.getType() == ofbx::Object::Type::CLUSTER){
		getClusterProperties(object);
	}else {
		int i = 0;
		while (ofbx::Object* child = object.resolveObjectLink(i))
		{
			getCluster(*child);
			++i;
		}			
	}
}

void getNodes(const ofbx::IScene& scene)
{
	const ofbx::Object* root = scene.getRoot();
	if (root) {
		int i = 0;
		while (ofbx::Object* child = root->resolveObjectLink(i))
		{
			getLimbNode(*child, nullptr);
			++i;
		}			
	}
	// getCluster(*root);
}

void loadFBX(const char* filepath)
{
	init(filepath);
	getNodes(*g_scene);
}

int main()
{
	SceneGraph& scene = SceneGraph::getInstance();

	scene.setGravity(Vec3f(0, 0, 0));
	// set scene
	std::shared_ptr<StaticBoundary<DataType3f>> root = scene.createNewScene<StaticBoundary<DataType3f>>();
	root->loadCube(Vec3f(0), Vec3f(1), 0.005f, true);
	root->varNormalFriction()->setValue(0.0f);
	// root->loadShpere(Vec3f(0.5, 0.7f, 0.5), 0.08f, 0.005f, false, true);

	// set dolphin
	std::shared_ptr<Dolphin<DataType3f>> dolphin = std::make_shared<Dolphin<DataType3f>>();
	temp_Dolphin = dolphin;
	root->addParticleSystem(dolphin); 

	dolphin->setMass(1.0f);

	// dolphin->loadMixFile("../../data/dolphin/Dolphin");
	// loadFBX("../../data/dolphin/Dolphin_Particles.fbx");
	// // 顺序：缩放，平移
	// dolphin->scale(0.2f); 
	// dolphin->translate(Vec3f(0.5f, 0.1f, 0.5f));
	
	//DEBUG 
	dolphin->loadParticles(Vec3f(-2, 0, -0.5), Vec3f(2, 4, 0.5), 0.08);
	loadFBX("../../data/dolphin/BoneBox.fbx");
	//顺序：缩放，平移
	dolphin->scale(0.1f); 
	dolphin->translate(Vec3f(0.5f, 0.4f, 0.25f));

	dolphin->setVisible(true);

	
	bool useVTK = true;
	RenderEngine* engine;

	engine = new GLRenderEngine;
	// set point 
	
	auto pointRenderer = std::make_shared<GLPointVisualModule>();
	pointRenderer->setColor(Vec3f(1, 0, 0));
	pointRenderer->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
	pointRenderer->setColorMapRange(0, 1);
	//dolphin->getSurfaceNode()->currentTopology()->connect(pointRenderer->inPointSet());
	//dolphin->getSurfaceNode()->graphicsPipeline()->pushModule(pointRenderer);

	dolphin->currentTopology()->connect(pointRenderer->inPointSet());
	// dolphin->currentVelocity()->connect(pointRenderer->inColor());
	dolphin->currentColor()->connect(pointRenderer->inColor()); //DEBUG

	dolphin->graphicsPipeline()->pushModule(pointRenderer);

	pointRenderer->setVisible(true);
	

	// set surface
	// auto sRender = std::make_shared<GLSurfaceVisualModule>();
	// sRender->setColor(Vec3f(1, 1, 0));
	// dolphin->getSurfaceNode()->currentTopology()->connect(sRender->inTriangleSet());
	// dolphin->getSurfaceNode()->graphicsPipeline()->pushModule(sRender);
	// sRender->setVisible(false);

	GlfwApp window;
	window.setRenderEngine(engine);
	window.createWindow(1024, 768);
	window.mainLoop();

	delete engine;
    return 0;
}