#pragma once
#include "JointTree.h"

namespace dyno
{
    IMPLEMENT_CLASS_1(JointTree, TDataType)

    template<typename TDataType>
    JointTree<TDataType>::JointTree()
    {
        id = -1;
        RotationActive = true;
        PreRotation = Coord(0);
        LclTranslation = Coord(0);
        LclRotation = Coord(0);
        LclScaling = Coord(0);
    }
    
    template<typename TDataType>
    JointTree<TDataType>::~JointTree()
    {

    }

    // TODO Global Local
    // Matrix Object::evalLocal(const Vec3& translation, const Vec3& rotation) const
    // {
    //     return evalLocal(translation, rotation, getLocalScaling());
    // }


    // Matrix Object::evalLocal(const Vec3& translation, const Vec3& rotation, const Vec3& scaling) const
    // {
    //     Vec3 rotation_pivot = getRotationPivot();
    //     Vec3 scaling_pivot = getScalingPivot();
    //     RotationOrder rotation_order = getRotationOrder();

    //     Matrix s = makeIdentity();
    //     s.m[0] = scaling.x;
    //     s.m[5] = scaling.y;
    //     s.m[10] = scaling.z;

    //     Matrix t = makeIdentity();
    //     setTranslation(translation, &t);

    //     Matrix r = getRotationMatrix(rotation, rotation_order);
    //     Matrix r_pre = getRotationMatrix(getPreRotation(), RotationOrder::EULER_XYZ);
    //     Matrix r_post_inv = getRotationMatrix(-getPostRotation(), RotationOrder::EULER_ZYX);

    //     Matrix r_off = makeIdentity();
    //     setTranslation(getRotationOffset(), &r_off);

    //     Matrix r_p = makeIdentity();
    //     setTranslation(rotation_pivot, &r_p);

    //     Matrix r_p_inv = makeIdentity();
    //     setTranslation(-rotation_pivot, &r_p_inv);

    //     Matrix s_off = makeIdentity();
    //     setTranslation(getScalingOffset(), &s_off);

    //     Matrix s_p = makeIdentity();
    //     setTranslation(scaling_pivot, &s_p);

    //     Matrix s_p_inv = makeIdentity();
    //     setTranslation(-scaling_pivot, &s_p_inv);

    //     // http://help.autodesk.com/view/FBX/2017/ENU/?guid=__files_GUID_10CDD63C_79C1_4F2D_BB28_AD2BE65A02ED_htm
    //     return t * r_off * r_p * r_pre * r * r_post_inv * r_p_inv * s_off * s_p * s * s_p_inv;
    // }

    // Matrix Object::getGlobalTransform() const
    // {
    //     const Object* parent = getParent();
    //     if (!parent) return evalLocal(getLocalTranslation(), getLocalRotation());

    //     return parent->getGlobalTransform() * evalLocal(getLocalTranslation(), getLocalRotation());
    // }


    // Matrix Object::getLocalTransform() const
    // {
    //     return evalLocal(getLocalTranslation(), getLocalRotation(), getLocalScaling());
    // }

    // TODO
	template<typename TDataType>
    Mat4f JointTree<TDataType>::getGlobalTransform()
    {
		return Mat4f(0.0);
    }

#ifdef PRECISION_FLOAT
	template class JointTree<DataType3f>;
#else
	template class JointTree<DataType3d>;
#endif    
}