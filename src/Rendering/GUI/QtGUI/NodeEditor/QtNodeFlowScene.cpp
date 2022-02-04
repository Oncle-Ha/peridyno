#include "QtNodeFlowScene.h"
#include "QtNodeWidget.h"

#include "nodes/Node"

#include "Object.h"
#include "NodeIterator.h"
#include "NodePort.h"
#include "SceneGraph.h"

#include <QtWidgets/QMessageBox>

namespace Qt
{

	QtNodeFlowScene::QtNodeFlowScene(std::shared_ptr<QtDataModelRegistry> registry, QObject* parent)
		: QtFlowScene(registry, parent)
	{

	}

	QtNodeFlowScene::QtNodeFlowScene(QObject* parent)
		: QtFlowScene(parent)
	{
		auto classMap = dyno::Object::getClassMap();

		auto ret = std::make_shared<QtDataModelRegistry>();
		int id = 0;
		for (auto const c : *classMap)
		{
			id++;

			QString str = QString::fromStdString(c.first);
			auto obj = dyno::Object::createObject(str.toStdString());
			std::shared_ptr<Node> node(dynamic_cast<Node*>(obj));

			if (node != nullptr)
			{
				QtDataModelRegistry::RegistryItemCreator creator = [str]() {
					auto node_obj = dyno::Object::createObject(str.toStdString());
					std::shared_ptr<Node> new_node(dynamic_cast<Node*>(node_obj));
					auto dat = std::make_unique<QtNodeWidget>(std::move(new_node));
					return dat;
				};

				QString category = "Default";// QString::fromStdString(module->getModuleType());
				ret->registerModel<QtNodeWidget>(category, creator);
			}
		}

		this->setRegistry(ret);

		dyno::SceneGraph& scn = dyno::SceneGraph::getInstance();
		showSceneGraph(&scn);

		connect(this, &QtFlowScene::nodeMoved, this, &QtNodeFlowScene::moveModulePosition);
	}


	QtNodeFlowScene::~QtNodeFlowScene()
	{
		clearScene();
	}


	void QtNodeFlowScene::showSceneGraph(SceneGraph* scn)
	{
		std::map<dyno::ObjectId, QtNode*> nodeMap;

		auto root = scn->getRootNode();

		SceneGraph::Iterator it_end(nullptr);

		auto addNodeWidget = [&](std::shared_ptr<Node> m) -> void
		{
			auto mId = m->objectId();

			auto type = std::make_unique<QtNodeWidget>(m);

			auto& node = this->createNode(std::move(type));

			nodeMap[mId] = &node;

			QPointF posView(m->bx(), m->by());

			node.nodeGraphicsObject().setPos(posView);

			this->nodePlaced(node);
		};

		for (auto it = scn->begin(); it != it_end; it++)
		{
			addNodeWidget(it.get());
		}

		auto createNodeConnections = [&](std::shared_ptr<Node> nd) -> void
		{
			auto inId = nd->objectId();

			if (nodeMap.find(inId) != nodeMap.end())
			{
				auto inBlock = nodeMap[nd->objectId()];

				auto ports = nd->getAllNodePorts();

				for (int i = 0; i < ports.size(); i++)
				{
					dyno::NodePortType pType = ports[i]->getPortType();
					if (dyno::Single == pType)
					{
						auto node = ports[i]->getNodes()[0];
						if (node != nullptr)
						{
							auto inBlock = nodeMap[node->objectId()];
							createConnection(*inBlock, 0, *inBlock, i);
						}
					}
					else if (dyno::Multiple == pType)
					{
						//TODO: a weird problem exist here, if the expression "auto& nodes = ports[i]->getNodes()" is used,
						//we still have to call clear to avoid memory leak.
						auto& nodes = ports[i]->getNodes();
						//ports[i]->clear();
						for (int j = 0; j < nodes.size(); j++)
						{
							if (nodes[j] != nullptr)
							{
								auto outId = nodes[j]->objectId();
								if (nodeMap.find(outId) != nodeMap.end())
								{
									auto outBlock = nodeMap[outId];
									createConnection(*inBlock, i, *outBlock, 0);
								}
							}
						}
						//nodes.clear();
					}
				}
			}
		};

		for (auto it = scn->begin(); it != it_end; it++)
		{
			createNodeConnections(it.get());
		}

		// 	clearScene();
		// 
		for (auto it = scn->begin(); it != it_end; it++)
		{
			auto node_ptr = it.get();
			std::cout << node_ptr->getClassInfo()->getClassName() << ": " << node_ptr.use_count() << std::endl;
		}
		nodeMap.clear();
	}

	void QtNodeFlowScene::moveModulePosition(QtNode& n, const QPointF& newLocation)
	{

	}
}