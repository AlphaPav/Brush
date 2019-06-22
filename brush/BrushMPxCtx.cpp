#include <Eigen/Core>
#include <Eigen/Dense>
#include "BrushMPxCtx.h"
#include <maya/MFloatPointArray.h>
#include <maya/MPointArray.h>
#include <maya/MPoint.h>
#include <maya/MFloatVectorArray.h>

BrushMPxCtx:: BrushMPxCtx() {
    pressure = 1.f;
    radius_scale_1 = 3.f;
    radius_scale_2 = radius_scale_1 * 2;
    radius_scale_3 = radius_scale_2 * 2;
    elastic_shear_modulus = 1;
    poisson_ratio = 0.4;
    brush_radius = 50;
    brushMode = mode_GRAB;
    is_2D = true;
    setTitleString("KevinLetBrush");
    setCursor(MCursor::handCursor);
    setImage("KevinLetBrush.xpm", MPxContext::kImage1);
    computeParams();
}

BrushMPxCtx:: ~BrushMPxCtx() {
    for (size_t i = 0; i < fnMesh.size(); i++)
        delete fnMesh[i];
}

void BrushMPxCtx:: setBrushMode(int brushMode) {
    this->brushMode = BRUSHMODE(brushMode);
    switch (brushMode) {
        case mode_GRAB:
            brushStringUI = "grab";
            break;
        case mode_GRAB_SYMMETRY:
            brushStringUI = "grab_symmetry";
            break;
        case mode_GRAB_BISCALE:
            brushStringUI = "grab_biscale";
            break;
        case mode_GRAB_TRISCALE:
            brushStringUI = "grab_triscale";
            break;
        case mode_TWIST:
            brushStringUI = "twist";
            break;
        case mode_TWIST_SYMMETRY:
            brushStringUI = "twist_symmetry";
            break;
        case mode_SCALE:
            brushStringUI = "scale";
            break;
        case mode_SCALE_SYMMETRY:
            brushStringUI = "scale_symmetry";
            break;
        case mode_PINCH:
            brushStringUI = "pinch";
            break;
        case mode_PINCH_SYMMETRY:
            brushStringUI = "pinch_symmetry";
            break;
        case mode_BRUSHSIZE:
            brushStringUI = "size_adjust";
            break;
    }
}


void BrushMPxCtx::setRadiusScale(float radiusScale, int index) {
    if (index == 1)
        radius_scale_1 = radiusScale;
    if (index == 2)
        radius_scale_2 = radiusScale;
    if (index == 3)
        radius_scale_3 = radiusScale;

    if (this->radius_scale_2 < this->radius_scale_1)
        this->radius_scale_2 = this->radius_scale_1 * 2;
    if (this->radius_scale_3 < this->radius_scale_2)
        this->radius_scale_3 = this->radius_scale_2 * 2;
}


void BrushMPxCtx::setElasticShearModulus(float elasticShearModulus) {
    this->elastic_shear_modulus = elasticShearModulus;
    computeParams();
}

void BrushMPxCtx::setPoissonRatio(float poissonRatio) {
    this->poisson_ratio = poissonRatio;
    computeParams();
}

void BrushMPxCtx::setPressure(float pressure) {
    this->pressure = pressure;
}

void BrushMPxCtx::setBrushRadius(double brushRadius) {
    this->brush_radius = brushRadius;
}

void BrushMPxCtx::setIs2DBrush(bool is2D) {
    this->is_2D = is2D;
    cout << "context2d:" << this->is_2D << endl;
}

void BrushMPxCtx::drawBrushUI(short mouse_x, short mouse_y, MHWRender::MUIDrawManager &drawMgr,
                 const MHWRender::MFrameContext &context) {
    drawMgr.beginDrawable();
    drawMgr.setColor(MColor(1.0f, 0.0f, 0.0f));
    drawMgr.circle2d(MPoint(mouse_x, mouse_y), brush_radius, false);
    drawMgr.setLineWidth(2.0);
    if (is_2D)
        drawMgr.text2d(MPoint(mouse_x, mouse_y), MString(brushStringUI + "2D"));
    else
        drawMgr.text2d(MPoint(mouse_x, mouse_y), MString(brushStringUI + "3D"));
    drawMgr.setPointSize(2.0);
    drawMgr.point(x0);
    drawMgr.endDrawable();
    view.refresh();
}


MVector BrushMPxCtx::computeForce(const MPoint &startScreenPos, const MPoint &endScreenPos) {
    MPoint pt1_start, pt2_start, pt1_end, pt2_end;
    MVector v1, v2;
    view.viewToWorld((int) startScreenPos.x, (int) startScreenPos.y, pt1_start, pt2_start);
    view.viewToWorld((int) endScreenPos.x, (int) endScreenPos.y, pt1_end, pt2_end);
    return MVector(pt2_end - pt2_start).normal();
}

void BrushMPxCtx::setMeshToManip(MStringArray meshNames) {
    this->mesh_to_manip = meshNames;
    size_t mesh_num = meshNames.length();
    fnMesh.resize(mesh_num);
    for (size_t i = 0; i < mesh_num; i++) {
        MSelectionList slist;
        MGlobal::getSelectionListByName(mesh_to_manip[i], slist);
        MGlobal::displayInfo(mesh_to_manip[i]);
        MDagPath dag;
        slist.getDagPath(0, dag);
        dag.extendToShape();
        fnMesh[i] = new MFnMesh(dag);
    }
}
void BrushMPxCtx::computeParams() {
    a = 1 / (4 * M_PI * elastic_shear_modulus);
    b = a / (4 * (1 - poisson_ratio));
    c = 2 / (3 * a - 2 * b);
}

MStatus
BrushMPxCtx::doPress(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context) {
    MGlobal::displayInfo("doPress!");
    // eye_position
    // look_direction
    MDagPath cameraDagPath;
    view.getCamera(cameraDagPath);
    MFnCamera fnCam(cameraDagPath);
    camera_pos = fnCam.eyePoint(MSpace::kWorld);
    MVector viewDir = fnCam.viewDirection(MSpace::kWorld);
    viewDir.normalize();
    is_ortho = fnCam.isOrtho();
    short x, y;
    event.getPosition(x, y);
    pressScreenPos = MPoint(x, y);
    MPoint raySource, rayEnd;
    view.viewToWorld(x, y, raySource, rayEnd);
    pressRay = rayEnd - raySource;
    pressRay.normalize();
    if (is_ortho)
        cout << "ortho:" << pressRay.x << "," << pressRay.y << "," << pressRay.z << endl;
    else
        cout << "persp:" << pressRay.x << "," << pressRay.y << "," << pressRay.z << endl;

    if (is_ortho && is_2D)
        camera_pos = raySource;
    int hitFace;
    bool hitted = false;
    float minHitParam = 100000;
    for (size_t i = 0; i < fnMesh.size(); i++) {
        MFloatPoint hitPoint;
        float hitParam;
        hitted = fnMesh[i]->closestIntersection(MFloatPoint(raySource), MFloatVector(rayEnd - raySource), NULL, NULL,
                                                false, MSpace::kWorld, 2, false, NULL, hitPoint, &hitParam, &hitFace,
                                                NULL, NULL, NULL);
        if (hitted && hitParam < minHitParam)
            x0 = hitPoint;
    }
    if (!hitted) {
        int conPointsCount = 0;
        MPoint conPointAverage(0, 0, 0);
        for (size_t i = 0; i < fnMesh.size(); i++) {
            MFloatPointArray allPoints;
            fnMesh[i]->getPoints(allPoints, MSpace::kObject);
            MFloatVectorArray allNormals;
            fnMesh[i]->getNormals(allNormals, MSpace::kObject);
            for (size_t j = 0; j < allPoints.length(); j++) {
                if (allNormals[i] * viewDir < 0)
                    continue;
                short temp_x, temp_y;
                view.worldToView(allPoints[j], temp_x, temp_y);
                MPoint screenPos(temp_x, temp_y);
                if (screenPos.distanceTo(pressScreenPos) <= brush_radius) {
                    conPointAverage = conPointAverage + allPoints[j];
                    conPointsCount++;
                }
            }

        }
        x0 = conPointAverage * (1.0 / conPointsCount);
    }
    drawBrushUI(x, y, drawMgr, context);
    return MS::kSuccess;
}

MStatus BrushMPxCtx::doRelease(MEvent &event, MHWRender::MUIDrawManager &drawMgr,
                               const MHWRender::MFrameContext &context) {
    MGlobal::displayInfo("doRelease!");
    // short x, y;
    // event.getPosition(x, y);
    //drawBrushUI(x,y,drawMgr,context);
    //view.refresh();
    return MS::kSuccess;
}

MStatus
BrushMPxCtx::doDrag(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context) {
    MGlobal::displayInfo("doDrag!");
    short x, y;
    event.getPosition(x, y);
    drawBrushUI(x, y, drawMgr, context);
    //lastScreenPos = currentScreenPos;
    currentScreenPos = MPoint(x, y);

    if (brushMode == mode_BRUSHSIZE) {
        if (currentScreenPos.x < pressScreenPos.x)
            brush_radius--;
        else
            brush_radius++;
        brush_radius = std::max(0.0, brush_radius);
        return MS::kSuccess;
    }

    //grab force
    MVector f = computeForce(pressScreenPos, currentScreenPos);
//    f.normalize();
    Eigen::Vector3f f_e(f.x, f.y, f.z);

    float f_length = f_e.norm();
    //twist force
    MPoint tp;
    MVector q;
    view.viewToWorld((int) pressScreenPos.x, (int) pressScreenPos.y, tp, q);
    q.normalize();
    q = q * f_length;
    if (currentScreenPos.x < pressScreenPos.x)
        q = -q;
    Eigen::Vector3f q_e(q.x, q.y, q.z);

    //scale direction
    float scaleDir = f_length;
    if (currentScreenPos.x < pressScreenPos.x)
        scaleDir = -f_length;

    //pinch matrix
    Eigen::Matrix3f pinchF;
    //    pinchF(0, 0) = 0;
    //    pinchF(1, 1) = 0;
    //    pinchF(2, 2) = 0;
    //    pinchF(0, 1) = f_e.z();
    //    pinchF(1, 0) = f_e.z();
    //    pinchF(0, 2) = f_e.y();
    //    pinchF(2, 0) = f_e.y();
    //    pinchF(1, 2) = f_e.x();
    //    pinchF(2, 1) = f_e.z();

    pinchF(0, 0) = f_e.x();
    pinchF(1, 1) = 0;
    pinchF(2, 2) = -f_e.x();
    pinchF(0, 1) = f_e.y();
    pinchF(1, 0) = f_e.y();
    pinchF(0, 2) = 0;
    pinchF(2, 0) = 0;
    pinchF(1, 2) = 0;
    pinchF(2, 1) = 0;
    Eigen::Matrix3f I;
    I = I.Identity();

    //symmtry
    Eigen::Vector3f plane_c(0., 0., 0.); //plane center
    Eigen::Vector3f plane_normal(1., 0., 0.);//plane nomal vector
    Eigen::Matrix3f M = Eigen::Matrix3f::Identity() - 2 * plane_normal * plane_normal.transpose();
    Eigen::Vector3f x0_e(x0.x, x0.y, x0.z);
    Eigen::Vector3f x1_e = plane_c + M * (x0_e - plane_c);
    MFloatPoint x1(x1_e.x(), x1_e.y(), x1_e.z());

    Eigen::Vector3f f_sym_e = M * f_e;
    f_sym_e = f_sym_e.transpose();
    f_e = f_e.transpose();

    Eigen::Matrix3f pinchF_sym = M * pinchF * M;

    for (size_t m = 0; m < fnMesh.size(); m++) {
        MFloatPointArray allPoints;
        fnMesh[m]->getPoints(allPoints, MSpace::kObject);
        MPointArray newPoints;
        for (size_t i = 0; i < allPoints.length(); i++) {
            short x, y;
            view.worldToView(allPoints[i], x, y);
            MPoint screenPos(x, y);
            //if (screenPos.distanceTo(lastScreenPos) <= brush_radius)
            {
                MVector vr;
                MVector vr_sym;
                if (!is_2D) {
                    vr = allPoints[i] - x0;
                    vr_sym = allPoints[i] - x1;
                } else {
                    if (!is_ortho) {
                        float projOnRay = MVector(allPoints[i] - camera_pos) * (pressRay);
                        vr = (allPoints[i] - camera_pos) - pressRay * projOnRay;
                    } else {
                        float projOnRay = MVector(allPoints[i] - camera_pos) * (pressRay);
                        vr = (allPoints[i] - camera_pos) - pressRay * projOnRay;
                    }
                }
                // float projOnRay = MVector(allPoints[i]-camera_pos)*(lastPointHitRay);
                // MVector vr = (allPoints[i]-camera_pos) - lastPointHitRay*projOnRay;
                Eigen::Vector3f r_e(vr.x, vr.y, vr.z);
                float r = vr.length();
                float re = sqrt(r * r + radius_scale_1 * radius_scale_1);
                Eigen::Vector3f r_sym_e(vr_sym.x, vr_sym.y, vr_sym.z);
                float r_sym = vr_sym.length();
                float re_sym = sqrt(r_sym * r_sym + radius_scale_1 * radius_scale_1);
                MVector dis;
                if (brushMode == mode_GRAB) {
                    Eigen::Matrix3f i = Eigen::Matrix3f::Identity();
                    Eigen::Matrix3f u = (a - b) / re * i + a * radius_scale_1 * radius_scale_1 / (2 * re * re * re) * i +
                                        b * r_e * r_e.transpose() / (re * re * re);

                    u *= c;
                    Eigen::Vector3f res = u * f_e * pressure;
                    dis.x = res.x();
                    dis.y = res.y();
                    dis.z = res.z();
                } else if (brushMode == mode_GRAB_SYMMETRY) {
                    Eigen::Matrix3f i = Eigen::Matrix3f::Identity();
                    Eigen::Matrix3f u = (a - b) / re * i + a * radius_scale_1 * radius_scale_1 / (2 * re * re * re) * i +
                                        b * r_e * r_e.transpose() / (re * re * re);
                    Eigen::Matrix3f u_sym = (a - b) / re_sym * i +
                                            a * radius_scale_1 * radius_scale_1 / (2 * re_sym * re_sym * re_sym) * i +
                                            b * r_sym_e * r_sym_e.transpose() / (re_sym * re_sym * re_sym);
                    u *= c;
                    u_sym *= c;
                    Eigen::Vector3f res = u * f_e * pressure;
                    Eigen::Vector3f res_sym = u_sym * f_sym_e * pressure;

                    dis.x = res.x() + res_sym.x();
                    dis.y = res.y() + res_sym.y();
                    dis.z = res.z() + res_sym.z();

                } else if (brushMode == mode_GRAB_BISCALE) {
                    Eigen::Matrix3f i = Eigen::Matrix3f::Identity();
                    float re_1 = re;
                    float re_2 = sqrt(r * r + radius_scale_2 * radius_scale_2);
                    Eigen::Matrix3f u_1 =
                            (a - b) / re_1 * i + a * radius_scale_1 * radius_scale_1 / (2 * re_1 * re_1 * re_1) * i +
                            b * r_e * r_e.transpose() / (re_1 * re_1 * re_1);
                    Eigen::Matrix3f u_2 =
                            (a - b) / re_2 * i + a * radius_scale_2 * radius_scale_2 / (2 * re_2 * re_2 * re_2) * i +
                            b * r_e * r_e.transpose() / (re_2 * re_2 * re_2);
                    Eigen::Matrix3f u = u_1 - u_2;

                    u *= (c / (1.0 / radius_scale_1 - 1 / radius_scale_2));
                    Eigen::Vector3f res = u * f_e * pressure;

                    dis.x = res.x();
                    dis.y = res.y();
                    dis.z = res.z();
                } else if (brushMode == mode_GRAB_TRISCALE) {
                    Eigen::Matrix3f i = Eigen::Matrix3f::Identity();
                    float re_1 = re;
                    float re_2 = sqrt(r * r + radius_scale_2 * radius_scale_2);
                    float re_3 = sqrt(r * r + radius_scale_3 * radius_scale_3);
                    Eigen::Matrix3f u_1 =
                            (a - b) * i / re_1 + a * radius_scale_1 * radius_scale_1 * i / (2 * re_1 * re_1 * re_1) +
                            b * r_e * r_e.transpose() / (re_1 * re_1 * re_1);
                    Eigen::Matrix3f u_2 =
                            (a - b) * i / re_2 + a * radius_scale_2 * radius_scale_2 * i / (2 * re_2 * re_2 * re_2) +
                            b * r_e * r_e.transpose() / (re_2 * re_2 * re_2);
                    Eigen::Matrix3f u_3 =
                            (a - b) * i / re_3 + a * radius_scale_3 * radius_scale_3 * i / (2 * re_3 * re_3 * re_3) +
                            b * r_e * r_e.transpose() / (re_3 * re_3 * re_3);
                    float w1 = 1;
                    float w2 = -(radius_scale_3 * radius_scale_3 - radius_scale_1 * radius_scale_1 + 0.0) /
                               (radius_scale_3 * radius_scale_3 - radius_scale_2 * radius_scale_2);
                    float w3 = (radius_scale_2 * radius_scale_2 - radius_scale_1 * radius_scale_1 + 0.0) /
                               (radius_scale_3 * radius_scale_3 - radius_scale_2 * radius_scale_2);
                    Eigen::Matrix3f u = u_1 * w1 + u_2 * w2 + u_3 * w3;
                    u *= (c / (w1 / radius_scale_1 + w2 / radius_scale_2 + w3 / radius_scale_3));
                    Eigen::Vector3f res = u * f_e * pressure;

                    dis.x = res.x();
                    dis.y = res.y();
                    dis.z = res.z();
                } else if (brushMode == mode_TWIST) {
                    float u = -a *
                              (1 / (re * re * re) + 3 * radius_scale_1 * radius_scale_1 / (2 * re * re * re * re * re));
//                    u *= c;
                    Eigen::Vector3f qr = q_e.cross(r_e);
                    Eigen::Vector3f dis_e = u * qr * pressure;
                    dis = MVector(dis_e.x(), dis_e.y(), dis_e.z());
                } else if (brushMode == mode_TWIST_SYMMETRY) {
                    float u = -a *
                              (1 / (re * re * re) + 3 * radius_scale_1 * radius_scale_1 / (2 * re * re * re * re * re));
                    float u_sym = -a * (1 / (re_sym * re_sym * re_sym) + 3 * radius_scale_1 * radius_scale_1 /
                                                                         (2 * re_sym * re_sym * re_sym * re_sym *
                                                                          re_sym));
                    Eigen::Vector3f qr = q_e.cross(r_e);
                    Eigen::Vector3f qr_e = -q_e.cross(r_sym_e);

                    Eigen::Vector3f dis_e = u * qr * pressure;
                    Eigen::Vector3f dis_sym = u_sym * qr_e * pressure;

                    dis = MVector(dis_e.x() + dis_sym.x(), dis_e.y() + dis_sym.y(), dis_e.z() + dis_sym.z());

                } else if (brushMode == mode_SCALE) {
                    float u = (2 * b - a) *
                              (1 / (re * re * re) + 3 * radius_scale_1 * radius_scale_1 / (2 * re * re * re * re * re));
                    u *= scaleDir;
                    dis = -vr * u;
                } else if (brushMode == mode_SCALE_SYMMETRY) {
                    float u = (2 * b - a) *
                              (1 / (re * re * re) + 3 * radius_scale_1 * radius_scale_1 / (2 * re * re * re * re * re));
                    float u_sym = (2 * b - a) * (1 / (re_sym * re_sym * re_sym) + 3 * radius_scale_1 * radius_scale_1 /
                                                                                  (2 * re_sym * re_sym * re_sym *
                                                                                   re_sym * re_sym));
                    u *= scaleDir;
                    u_sym *= scaleDir;
                    dis = -(vr * u + vr_sym * u_sym);

                } else if (brushMode == mode_PINCH) {
                    Eigen::Matrix3f i = Eigen::Matrix3f::Identity();
                    Eigen::Vector3f fr = pinchF * r_e;
                    Eigen::Matrix3f rfr = r_e.transpose() * pinchF * r_e * i;

                    Eigen::Vector3f res = (2 * b - a) / (re * re * re) * fr + 3 / (2 * re * re * re * re * re) *
                                                                              (2 * b * rfr +
                                                                               a * radius_scale_1 * radius_scale_1 *
                                                                               pinchF) * r_e;
                    res = -res * pressure;
                    dis = MVector(res.x(), res.y(), res.z());
                    dis *= c;
                } else if (brushMode == mode_PINCH_SYMMETRY) {
                    Eigen::Matrix3f i = Eigen::Matrix3f::Identity();
                    Eigen::Vector3f fr = pinchF * r_e;
                    Eigen::Vector3f fr_sym = pinchF_sym * r_sym_e;
                    Eigen::Matrix3f rfr = r_e.transpose() * pinchF * r_e * i;
                    Eigen::Matrix3f rfr_sym = r_sym_e.transpose() * pinchF_sym * r_sym_e * i;

                    Eigen::Vector3f res = (2 * b - a) / (re * re * re) * fr + 3 / (2 * re * re * re * re * re) *
                                                                              (2 * b * rfr +
                                                                               a * radius_scale_1 * radius_scale_1 *
                                                                               pinchF) * r_e;

                    Eigen::Vector3f res_sym = (2 * b - a) / (re_sym * re_sym * re_sym) * fr_sym +
                                              3 / (2 * re_sym * re_sym * re_sym * re_sym * re_sym) *
                                              (2 * b * rfr_sym + a * radius_scale_1 * radius_scale_1 * pinchF_sym) *
                                              r_sym_e;
                    res = -res * pressure;
                    res_sym = -res_sym * pressure;
                    dis = MVector(res.x() + res_sym.x(), res.y() + res_sym.y(), res.z() + res_sym.z());
                    dis *= c;
                }
                newPoints.append(allPoints[i] + dis);
            }
        }
        fnMesh[m]->setPoints(newPoints);
    }
    return MS::kSuccess;
}


MStatus
BrushMPxCtx::doEnterRegion(MEvent &event, MHWRender::MUIDrawManager &drawMgr, const MHWRender::MFrameContext &context) {
    MGlobal::displayInfo("doEnterRegion");
    view = M3dView::active3dView();
    return MS::kSuccess;
}

MStatus BrushMPxCtx::doEnterRegion(MEvent &event) {
    MGlobal::displayInfo("doEnterRegion");
    view = M3dView::active3dView();
    return MS::kSuccess;
}

MStatus BrushMPxCtx::doHold(MEvent &event) {
    MGlobal::displayInfo("doHold");
    return MS::kSuccess;
}

void BrushMPxCtx::getClassName(MString &name) const {
    name.set("KelvinBrushPaint");
}
