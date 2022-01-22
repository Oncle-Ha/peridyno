#pragma once
#include "Module/TopologyModule.h"

namespace dyno
{
	/*!
	*	\class	JointTree
	*	\brief	A JointTree(Skeleton) represents a hierarchical tree structure of joints
	*/    
    template<typename TDataType>
	class JointTree : public TopologyModule
	{ 
        DECLARE_CLASS_1(JointTree, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TDataType::Matrix Matrix;

        JointTree();
        ~JointTree();

        void copyFrom(JointTree<TDataType>& jointTree);
        void getGlobalTransform();
        Mat4f getLocalTransform();
        
        unsigned long long id;
        Coord PreRotation;
        Coord LclTranslation;
        Coord LclRotation;
        Coord LclScaling; 
        Mat4f GlobalTransform;
        bool RotationActive;
        std::vector<std::shared_ptr<JointTree>> children;
        std::shared_ptr<JointTree> parent;
    };
}