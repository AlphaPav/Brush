#pragma once

#include "BrushMPxToolCmd.h"
#include <maya/MPoint.h>
#include <maya/MCursor.h>
#include <maya/M3dView.h>
#include <maya/MFnMesh.h>
#include <maya/MSelectionList.h>
#include <maya/MGlobal.h>
#include <maya/MDagPath.h>
#include <maya/MFloatPoint.h>
#include <maya/MFnCamera.h>
#include <maya/MPxContext.h>
#include <maya/MPxToolCommand.h>
#include <maya/MEvent.h>
#include <maya/MUiMessage.h>
#include <maya/MUIDrawManager.h>
#include <maya/MFrameContext.h>
#include <maya/MPoint.h>
#include <maya/MColor.h>
#include <maya/MAngle.h>
#include <maya/MString.h>
#include <maya/MGlobal.h>
#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MArgList.h>
#include <maya/MItSelectionList.h>
#include <maya/MSelectionList.h>
#include <maya/MIOStream.h>
#include <math.h>
#include <stdlib.h>

#if defined(__APPLE__)

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

//COMPUTE
#include <maya/MFnDagNode.h>
#include <maya/MDagPathArray.h>
#include <maya/MDagPath.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MPointArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFnTransform.h>
#include <maya/MFnMesh.h>
#include <vector>

class BrushMPxCtx : public MPxContext {
public:
    //float radius;//r
    float pressure;//|f|
    float radius_scale_1;//epsilon
    float radius_scale_2;//epsilon
    float radius_scale_3;//epsilon
    float elastic_shear_modulus;//miu
    float poisson_ratio;//v
    double brush_radius;
    MStringArray mesh_to_manip;
    enum BRUSHMODE {
        mode_GRAB, mode_GRAB_SYMMETRY,
        mode_GRAB_BISCALE, mode_GRAB_TRISCALE,
        mode_TWIST, mode_TWIST_SYMMETRY,
        mode_SCALE, mode_SCALE_SYMMETRY,
        mode_PINCH, mode_PINCH_SYMMETRY,
        mode_BRUSHSIZE
    };
    BRUSHMODE brushMode;
    bool is_2D;
    MString brushStringUI;
    //params to be computed
    float a;//1/4/pi/elasticShearModulus
    float b;//a/4/(1-poissonRatio)
    float c;//2/(3a-2b)
    MPoint pressScreenPos;
    MPoint currentScreenPos;
    // MFloatPoint lastPointOnMesh;
    MVector pressRay;
    MFloatPoint x0;//concentracted part's center
    //  MFloatPoint x1; // symmetry point
    M3dView view;
    std::vector<MFnMesh *> fnMesh;
    MFloatPoint camera_pos;
    bool is_ortho;

    void setPressure(float pressure) ;
    void setBrushRadius(double brushRadius) ;
    void setIs2DBrush(bool is2D) ;
    void setBrushMode(int brushMode);
    void setRadiusScale(float radiusScale, int index);
    void setElasticShearModulus(float elasticShearModulus) ;
    void setPoissonRatio(float poissonRatio) ;
    void setMeshToManip(MStringArray meshNames) ;
    void computeParams();
    void drawBrushUI(short mouse_x, short mouse_y, MHWRender::MUIDrawManager &drawMgr,
                     const MHWRender::MFrameContext &context) ;

    BrushMPxCtx();
    virtual ~BrushMPxCtx();

    // 继承Maya类MPxContext中的函数
    virtual void toolOnSetup(MEvent &event) {}

    virtual MStatus doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context);

    virtual MStatus doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context);

    virtual MStatus
    doRelease(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context);

    virtual MStatus
    doEnterRegion(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context) ;

    virtual MStatus doEnterRegion(MEvent &event) ;

    virtual MStatus doHold(MEvent &event) ;

    virtual void getClassName(MString &name) const ;


private:
    MVector computeForce(const MPoint &startScreenPos, const MPoint &endScreenPos) ;
};

