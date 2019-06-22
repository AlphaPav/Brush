#pragma once

#include <maya/MPxToolCommand.h>
#include <maya/MPxContext.h>
#include <maya/MPxContextCommand.h>
#include <maya/MGlobal.h>
#include <maya/MArgList.h>

class BrushMPxToolCmd : public MPxToolCommand
{
	// 继承Maya类MPxToolCommand
public:
	BrushMPxToolCmd() { setCommandString("KevinLetBrush"); }
	virtual ~BrushMPxToolCmd() {}

	static void* creator() { return new BrushMPxToolCmd; }
	MStatus doIt(const MArgList& args) { return redoIt(); }
	MStatus redoIt();
	MStatus undoIt();
	bool isUndoable()const { return true; }

	MStatus finalize()
	{
		MArgList command;
		command.addArg(commandString());
		return MPxToolCommand::doFinalize(command);
	}

};







