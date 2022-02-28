#pragma once
#include "Module/TopologyModule.h"

namespace dyno
{

    struct JCapsule
    {
        typedef typename DataType3f::Coord Coord;
        int id_joint;
        int id_cap;
        Coord v0,v1;

        DYN_FUNC JCapsule(){};
        DYN_FUNC JCapsule(int id1, int id2, Coord a0, Coord a1): 
        id_joint(id1), id_cap(id2), v0(a0), v1(a1){};
        DYN_FUNC ~JCapsule(){};
    };

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

        void copyFrom(JointTree<TDataType> jointTree);

        void getGlobalTransform();
        Mat4f getLocalTransform();
        
        
        unsigned long long id;
        Coord PreRotation;
        Coord LclTranslation;   // Local Joint's coord
        Coord LclRotation;
        Coord LclScaling; 
        Coord GlCoord;          // Global Joint's coord 
        Mat4f GlobalTransform;
        bool RotationActive;
        std::vector<std::shared_ptr<JointTree>> children;
        std::shared_ptr<JointTree> parent;
    };
}