#include "BrushMPxCtxCmd.h"

MStatus BrushMPxCtxCmd::appendSyntax() {
    //和Mel脚本中变量映射
    MSyntax mySyntax = syntax();
    MStatus stat = mySyntax.addFlag("-rs1", "-radiusScale1", MSyntax::kDouble);
    stat = mySyntax.addFlag("-rs2", "-radiusScale2", MSyntax::kDouble);
    stat = mySyntax.addFlag("-rs3", "-radiusScale3", MSyntax::kDouble);
    stat = mySyntax.addFlag("-ps", "-pressure", MSyntax::kDouble);
    stat = mySyntax.addFlag("-esm", "-elasticShearModulus", MSyntax::kDouble);
    stat = mySyntax.addFlag("-pr", "-poissonRatio", MSyntax::kDouble);
    stat = mySyntax.addFlag("-bm", "-brushMode", MSyntax::kString);
    stat = mySyntax.addFlag("-2d", "-is2DBrush", MSyntax::kBoolean);
    stat = mySyntax.addFlag("-bs", "-brushSize", MSyntax::kDouble);
    stat = mySyntax.addFlag("-m", "-mesh", MSyntax::kString);
    //mySyntax.useSelectionAsDefault(true);
    //stat = mySyntax.setObjectType(MSyntax::kSelectionList, 1);
    return MS::kSuccess;
}


MStatus BrushMPxCtxCmd::doEditFlags() {
    MGlobal::displayInfo("doEditFlags!");
    MStatus status = MS::kSuccess;
    MArgParser argData = parser();
    //用来查询参数变量，并保存到自己的数据结构
    if (argData.isFlagSet("-rs1")) {
        double radiusScale;
        radiusScale = argData.flagArgumentDouble("-rs1", 0);
        m_KevinLetBrushCtx->setRadiusScale(radiusScale, 1);
    }

    if (argData.isFlagSet("-rs2")) {
        double radiusScale;
        radiusScale = argData.flagArgumentDouble("-rs2", 0);
        m_KevinLetBrushCtx->setRadiusScale(radiusScale, 2);
    }

    if (argData.isFlagSet("-rs3")) {
        double radiusScale;
        radiusScale = argData.flagArgumentDouble("-rs3", 0);
        m_KevinLetBrushCtx->setRadiusScale(radiusScale, 3);
    }

    if (argData.isFlagSet("-ps")) {
        double pressure;
        argData.getFlagArgument("-ps", 0, pressure);
        m_KevinLetBrushCtx->setPressure(pressure);
    }

    if (argData.isFlagSet("-bs")) {
        double brushSize;
        argData.getFlagArgument("-bs", 0, brushSize);
        m_KevinLetBrushCtx->setBrushRadius(brushSize);
    }

    if (argData.isFlagSet("-esm")) {
        double elasticShearModulus;
        argData.getFlagArgument("-esm", 0, elasticShearModulus);
        m_KevinLetBrushCtx->setElasticShearModulus(elasticShearModulus);
    }

    if (argData.isFlagSet("-pr")) {
        double poissonRatio;
        argData.getFlagArgument("-pr", 0, poissonRatio);
        m_KevinLetBrushCtx->setPoissonRatio(poissonRatio);
    }

    if (argData.isFlagSet("-bm")) {
        MString brushMode;
        argData.getFlagArgument("-bm", 0, brushMode);
        if (brushMode == "Grab")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_GRAB);
        else if (brushMode == "Grab_Symmetry")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_GRAB_SYMMETRY);
        else if (brushMode == "Grab_BiScale")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_GRAB_BISCALE);
        else if (brushMode == "Grab_TriScale")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_GRAB_TRISCALE);
        else if (brushMode == "Twist")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_TWIST);
        else if (brushMode == "Twist_Symmetry")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_TWIST_SYMMETRY);
        else if (brushMode == "Scale")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_SCALE);
        else if (brushMode == "Scale_Symmetry")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_SCALE_SYMMETRY);
        else if (brushMode == "Pinch")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_PINCH);
        else if (brushMode == "Pinch_Symmetry")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_PINCH_SYMMETRY);
        else if (brushMode == "Size_Adjust")
            m_KevinLetBrushCtx->setBrushMode(BrushMPxCtx::mode_BRUSHSIZE);

    }

    if (argData.isFlagSet("-2d")) {
        bool is2DBrush;
        argData.getFlagArgument("-2d", 0, is2DBrush);
        cout << "2d:" << is2DBrush << endl;
        m_KevinLetBrushCtx->setIs2DBrush(is2DBrush);
    }

    MSelectionList sList;
    MGlobal::getActiveSelectionList(sList);
    MStringArray allSelectedMesh;
    for (unsigned int kk = 0; kk < sList.length(); ++kk) {
        MDagPath dag;
        sList.getDagPath(kk, dag);
        if (dag.hasFn(MFn::kMesh)) {
            allSelectedMesh.append(dag.fullPathName());
        }
    }
    m_KevinLetBrushCtx->setMeshToManip(allSelectedMesh);

    return MS::kSuccess;
}

MStatus BrushMPxCtxCmd::doQueryFlags() {
    // 用来查询参数变量的值
    MArgParser argData = parser();
    if (argData.isFlagSet("-rs1")) {
        setResult(m_KevinLetBrushCtx->radius_scale_1);
    }
    if (argData.isFlagSet("-rs2")) {
        setResult(m_KevinLetBrushCtx->radius_scale_2);
    }
    if (argData.isFlagSet("-rs3")) {
        setResult(m_KevinLetBrushCtx->radius_scale_3);
    }

    if (argData.isFlagSet("-ps")) {
        setResult(m_KevinLetBrushCtx->pressure);
    }

    if (argData.isFlagSet("-esm")) {
        setResult(m_KevinLetBrushCtx->elastic_shear_modulus);
    }

    if (argData.isFlagSet("-pr")) {
        setResult(m_KevinLetBrushCtx->poisson_ratio);
    }

    if (argData.isFlagSet("-bs")) {
        if (m_KevinLetBrushCtx)
            setResult(m_KevinLetBrushCtx->brush_radius);
        else
            setResult(50.0);
    }

    if (argData.isFlagSet("-m")) {
        MString res;
        for (size_t i = 0; i < m_KevinLetBrushCtx->mesh_to_manip.length(); i++)
            res += (m_KevinLetBrushCtx->mesh_to_manip[i] + " , ");
        setResult(res);
    }

    if (argData.isFlagSet("-bm")) {
        if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_GRAB)
            setResult("Grab");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_GRAB_SYMMETRY)
            setResult("Grab_Symmetry");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_GRAB_BISCALE)
            setResult("Grab_BiScale");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_GRAB_TRISCALE)
            setResult("Grab_TriScale");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_TWIST)
            setResult("Twist");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_TWIST_SYMMETRY)
            setResult("Twist_Symmetry");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_SCALE)
            setResult("Scale");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_SCALE_SYMMETRY)
            setResult("Scale_Symmetry");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_PINCH)
            setResult("Pinch");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_PINCH_SYMMETRY)
            setResult("Pinch_Symmetry");
        else if (m_KevinLetBrushCtx->brushMode == BrushMPxCtx::mode_BRUSHSIZE)
            setResult("Size_Adjust");
    }

    if (argData.isFlagSet("-2d")) {
        setResult(m_KevinLetBrushCtx->is_2D);
    }

    return MS::kSuccess;
}

