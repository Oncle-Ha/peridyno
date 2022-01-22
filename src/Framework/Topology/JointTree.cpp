#pragma once
#include "JointTree.h"

#define cos_angle(angle) cos(angle / 180.0f * 3.1415926)
#define sin_angle(angle) sin(angle / 180.0f * 3.1415926)

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

    // TODO Local
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

	template<typename TDataType>
    Mat4f JointTree<TDataType>::getLocalTransform()
    {
        Mat4f translation = Mat4f(
            0, 0, 0, LclTranslation[0],
            0, 0, 0, LclTranslation[1],
            0, 0, 0, LclTranslation[2],
            0, 0, 0, 1);
        
        float P = LclRotation[1];
        Mat4f rotation_x = Mat4f(
            1, 0, 0, 0,
            0, cos_angle(P), -sin_angle(P), 0,
            0, sin_angle(P), cos_angle(P), 0,
            0, 0, 0, 1);

        float H = LclRotation[0];
        Mat4f rotation_y = Mat4f(
            cos_angle(H), 0, sin_angle(H), 0,
            0, 1, 0, 0,
            -sin_angle(H), 0, cos_angle(H), 0,
            0, 0, 0, 1);

        float B = LclRotation[2];
        Mat4f rotation_z = Mat4f(
            cos_angle(B), -sin_angle(B), 0, 0,
            sin_angle(B), cos_angle(B), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1);

        Mat4f scaling= Mat4f(
            LclScaling[0], 0, 0, 0,
            0, LclScaling[1], 0, 0,
            0, 0, LclScaling[2], 0,
            0, 0, 0, 1);

        return  translation * scaling * rotation_x * rotation_y * rotation_z;
        
    }
    // 遍历关节层次时，顺便更新
	template<typename TDataType>
    void JointTree<TDataType>::getGlobalTransform()
    {
        // 注意顺序
        this->GlobalTransform = getLocalTransform();
        if(this->parent != nullptr)
            this->GlobalTransform = this->parent->GlobalTransform * this->GlobalTransform;
    }

#ifdef PRECISION_FLOAT
	template class JointTree<DataType3f>;
#else
	template class JointTree<DataType3d>;
#endif    
}