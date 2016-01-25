#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "Block.hpp"
#include "Vector2D.hpp"
#include "Config.h"


class ToolPathTimeEstimator
{
public:
	ToolPathTimeEstimator();
	~ToolPathTimeEstimator();
	double calculate(FILE * out);

	void detail(FILE * out);
	void timeline(FILE * out, float resolution);
	void prepare();
	void addJump(double x, double y);
	void addBlock(double x, double y);
    Vector2D calcJunctionSpeed(Block* A, Block* B);
    double calcBlockBuildTime(Block* blockPtr);
    void backwardIterate(Block* blockPtr);
    void forwardIterate(Block* blockPtr);
	Vector2D calcAllowedSpeed(double acceleration, double target_speed, Vector2D length);
	double length;
	void updateNominalSpeed();
	Vector2D nominal_speed;
    Vector2D allow_speed;
	Vector2D acceleration;
	Vector2D jerk;
	double distance(double acceleration, double entry_speed, double exit_float);
	struct List{
		Block* blockHead;
		Block* blockTail;
	}blocks;
};