#pragma once
#include "BrushMPxCtx.h"
#include <maya/MFnPlugin.h>

class BrushMPxCtxCmd :public MPxContextCommand
{
	// 继承Maya类MPxContextCommand
public:
	BrushMPxCtxCmd(){}
	virtual MStatus         doEditFlags();
	virtual MStatus         doQueryFlags();
	virtual MPxContext* 	makeObj()
	{
		MGlobal::displayInfo("makeobj!");
		m_KevinLetBrushCtx = new BrushMPxCtx();
		return m_KevinLetBrushCtx;
	}
	virtual MStatus         appendSyntax();
	static void*            creator()
	{
		return new BrushMPxCtxCmd;
	}

protected:
BrushMPxCtx* m_KevinLetBrushCtx;
};

MStatus initializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj, PLUGIN_COMPANY, "3.0", "Any");
	status = plugin.registerContextCommand("KevinLetBrushCtx",
			BrushMPxCtxCmd::creator,
		"KevinLetBrush",
		BrushMPxToolCmd::creator);
	if (!status) {
		status.perror("registerContextCommand");
		return status;
	}

	return status;
}

MStatus uninitializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj);
	status = plugin.deregisterContextCommand("KevinLetBrushCtx",
		"KevinLetBrush");
	if (!status) {
		status.perror("deregisterContextCommand");
		return status;
	}

	return status;
}


