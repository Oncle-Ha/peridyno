#pragma once
#include "JointTree.h"

#define cos_angle(angle) cos(double(angle * 3.1415926f / 180.0))
#define sin_angle(angle) sin(double(angle * 3.1415926f / 180.0))

namespace dyno
{
    IMPLEMENT_CLASS_1(JointTree, TDataType)
    
    template<typename TDataType>
    JointTree<TDataType>::JointTree()
    {
        id = -1;
        RotationActive = true;
        PreRotation = Coord(0);
        PreTranslation = Coord(0);
        PreScaling = Coord(1); 

        LclTranslation = Coord(0);
        LclRotation = Coord(0);
        LclScaling = Coord(1);
		
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
    Mat4f JointTree<TDataType>::getTransform(Coord & T, Coord& R, Coord& S)
    {
		// TODO:改为四元数
        Mat4f translation = Mat4f(
            1, 0, 0, T[0],
            0, 1, 0, T[1],
            0, 0, 1, T[2],
            0, 0, 0, 1);
        
        double X = R[0];
        Mat4f rotation_x = Mat4f(
            1, 0, 0, 0,
            0, cos_angle(X), -sin_angle(X), 0,
            0, sin_angle(X), cos_angle(X), 0,
            0, 0, 0, 1);

		double Y = R[1];
        Mat4f rotation_y = Mat4f(
            cos_angle(Y), 0, sin_angle(Y), 0,
            0, 1, 0, 0,
            -sin_angle(Y), 0, cos_angle(Y), 0,
            0, 0, 0, 1);

		double Z = R[2];
        Mat4f rotation_z = Mat4f(
            cos_angle(Z), -sin_angle(Z), 0, 0,
            sin_angle(Z), cos_angle(Z), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1);

        Mat4f scaling= Mat4f(
            S[0], 0, 0, 0,
            0, S[1], 0, 0,
            0, 0, S[2], 0,
            0, 0, 0, 1);

        return  translation * scaling * rotation_x * rotation_y * rotation_z;
    }

    // 遍历关节层次时，顺便更新
	template<typename TDataType>
    void JointTree<TDataType>::getGlobalTransform()
    {
        // 注意顺序
        this->GlobalTransform = getTransform(this->LclTranslation, this->LclRotation, this->LclScaling);
        if(this->parent != nullptr)
            this->GlobalTransform = this->parent->GlobalTransform * this->GlobalTransform;
        else
        {
            //DEBUG 简单平移
            this->GlobalTransform = getTransform(this->PreTranslation, this->PreRotation, this->PreScaling) * this->GlobalTransform;
             
            Coord tmp_d(0,0,0);
            for (int i = 0; i < 3; ++i)
                tmp_d[i] += dist(generator)  * 0.0005;
            this->PreTranslation += tmp_d;
            
            // float s = 10000;
            // Mat4f scaling= Mat4f(
            // s, 0, 0, 0,
            // 0, s, 0, 0,
            // 0, 0, s, 0,
            // 0, 0, 0, 1);
            // this->GlobalTransform = scaling * this->GlobalTransform;
        }
            
    }


    template<typename TDataType>
    void JointTree<TDataType>::copyFrom(JointTree<TDataType> jointTree)
    {
        this->id = jointTree.id;
        this->PreRotation = jointTree.PreRotation;
        this->LclTranslation = jointTree.LclTranslation;
        this->LclRotation = jointTree.LclRotation;
        this->LclScaling = jointTree.LclScaling;
        this->GlCoord = jointTree.GlCoord;
        this->GlobalTransform = jointTree.GlobalTransform;
        this->RotationActive = jointTree.RotationActive;
        this->children.assign(jointTree.children.begin(), jointTree.children.end());
        this->parent = jointTree.parent;
    }

    template<typename TDataType>
    void JointTree<TDataType>::scale(Real s)
    {
		PreScaling *= s;
    }

    template<typename TDataType>
    void JointTree<TDataType>::translate(Coord t)
    {
        PreTranslation += t;
    }
    
#ifdef PRECISION_FLOAT
	template class JointTree<DataType3f>;
#else
	template class JointTree<DataType3d>;
#endif    
}