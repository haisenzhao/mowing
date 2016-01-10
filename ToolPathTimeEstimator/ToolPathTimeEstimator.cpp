#include "StdAfx.h"
#include "ToolPathTimeEstimator.hpp"
#define π 3.141592653

ToolPathTimeEstimator::ToolPathTimeEstimator() {
    this->nominal_speed = Vector2D(0, 0);
    this->blocks.blockHead = new Block(0.0, 0.0);
    this->blocks.blockTail = this->blocks.blockHead;
    this->acceleration = Vector2D(3000, 3000);
    this->nominal_speed = Vector2D(DEFAULT_NOMINALSPEED, DEFAULT_NOMINALSPEED);
	this->nominal_speed.direction = false;
	this->jerk = Vector2D(20, 20);
}

ToolPathTimeEstimator::~ToolPathTimeEstimator() {
}

void ToolPathTimeEstimator::prepare() {
	this->blocks.blockHead = this->blocks.blockHead->next;
}

void ToolPathTimeEstimator::addJump(double x, double y)
{
	if ((x == this->blocks.blockTail->x) && (y == this->blocks.blockTail->y)) {
		return;
	}
	Block* T = new Block(x, y);
	T->previous = this->blocks.blockTail;
	this->blocks.blockTail->next = T;
	T->length = Vector2D(x - T->previous->x, y - T->previous->y);
	T->nominal_speed = Vector2D(JUMPSPEED, JUMPSPEED);
	this->blocks.blockTail = T;
}

void ToolPathTimeEstimator::addBlock(double x, double y) {
	if ((x == this->blocks.blockTail->x) && (y == this->blocks.blockTail->y)) {
		return;
	}
    Block* T = new Block(x,y);
    T->previous = this->blocks.blockTail;
    this->blocks.blockTail->next = T;
    T->length = Vector2D(x-T->previous->x, y-T->previous->y);
	T->nominal_speed = this->nominal_speed;
	this->blocks.blockTail = T;
}

double ToolPathTimeEstimator::calculate(FILE * out) {
    Block* blockPtr;
    //Processing the 1..n-1 block
	blockPtr = this->blocks.blockHead->next;

    this->blocks.blockHead->exit_speed = Vector2D(0, 0);
    while (blockPtr->next != this->blocks.blockTail) {
        blockPtr->entry_speed = Vector2D::getVectorFromDirection(blockPtr->previous->exit_speed.value, blockPtr->length);
        blockPtr->exit_speed = calcJunctionSpeed(blockPtr, blockPtr->next);
        if (blockPtr->exit_speed.value < calcAllowedSpeed(-acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, blockPtr->length).value) {
            allow_speed = Vector2D::getVectorFromDirection(calcAllowedSpeed(acceleration.maxAlongDirection(blockPtr->length), blockPtr->exit_speed.value, blockPtr->length).value, blockPtr->length);
            if (allow_speed.value < blockPtr->entry_speed.value) {
                backwardIterate(blockPtr);
			}
			else {
				exit(1);
			}
		}
		blockPtr->build_time = calcBlockBuildTime(blockPtr);
		blockPtr = blockPtr->next;
	}

	//Processing the n block
	blockPtr->entry_speed = blockPtr->previous->exit_speed;
	allow_speed = Vector2D::getVectorFromDirection(sqrt(2 * acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value), blockPtr->length);
	blockPtr->exit_speed = Vector2D(0, 0);
	if (blockPtr->exit_speed.value < calcAllowedSpeed(-acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, blockPtr->length).value) {
		allow_speed = Vector2D::getVectorFromDirection(sqrt(2 * acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value + pow(blockPtr->exit_speed.value, 2)), blockPtr->length);
		if (allow_speed.value < blockPtr->entry_speed.value) {
			backwardIterate(blockPtr);
		}
		else {
			exit(1);
		}
	}
	blockPtr->build_time = calcBlockBuildTime(blockPtr);

	//Sum up the build time
	double buildTime = 0;
	blockPtr = this->blocks.blockHead->next;
	while (blockPtr->next != nullptr) {
		blockPtr->build_time = calcBlockBuildTime(blockPtr);
		buildTime += blockPtr->build_time;
		length += blockPtr->length.value;
		blockPtr = blockPtr->next;
		fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, blockPtr->entry_speed.value, blockPtr->exit_speed.value, blockPtr->length.value, blockPtr->build_time);
	}

	return buildTime;


}

double ToolPathTimeEstimator::calcBlockBuildTime(Block* blockPtr) {
	if (blockPtr->length.value <= ((2 * pow(blockPtr->nominal_speed.maxAlongDirection(blockPtr->length), 2) - pow(blockPtr->entry_speed.value, 2) - pow(blockPtr->exit_speed.value, 2)) / 2 / acceleration.maxAlongDirection(blockPtr->length))) {
		blockPtr->max_speed = Vector2D::getVectorFromDirection(sqrt((2 * acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value + pow(blockPtr->exit_speed.value, 2) + pow(blockPtr->entry_speed.value, 2)) / 2), blockPtr->length);
		return (2 * blockPtr->max_speed.value - blockPtr->entry_speed.value - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length);

	}
	else {
		return (2 * blockPtr->nominal_speed.maxAlongDirection(blockPtr->length) - blockPtr->entry_speed.value - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length) + (blockPtr->length.value - (2 * pow(blockPtr->nominal_speed.maxAlongDirection(blockPtr->length), 2) - pow(blockPtr->entry_speed.value, 2) - pow(blockPtr->exit_speed.value, 2)) / 2 / acceleration.maxAlongDirection(blockPtr->length)) / blockPtr->nominal_speed.maxAlongDirection(blockPtr->length);
	}
}

void ToolPathTimeEstimator::backwardIterate(Block *blockPtr) {
	if (calcJunctionSpeed(blockPtr->previous, blockPtr).value < calcAllowedSpeed(-acceleration.maxAlongDirection(blockPtr->length), blockPtr->exit_speed.value, blockPtr->length).value) {
		allow_speed = Vector2D::getVectorFromDirection(sqrt(2 * acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value + pow(calcJunctionSpeed(blockPtr->previous, blockPtr).value, 2)), blockPtr->length);
		if (allow_speed.value > Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length).value) {
			allow_speed = Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length);
		}
		if (blockPtr->exit_speed.value <= allow_speed.value) {
			return;
		}
		else {
			forwardIterate(blockPtr);
		}
	}
	else {
		blockPtr->entry_speed = allow_speed;
		blockPtr = blockPtr->previous;
		blockPtr->exit_speed = Vector2D::getVectorFromDirection(allow_speed.value, blockPtr->length);
		allow_speed = Vector2D::getVectorFromDirection(calcAllowedSpeed(acceleration.maxAlongDirection(blockPtr->length), blockPtr->exit_speed.value, blockPtr->length).value, blockPtr->length);
		if (blockPtr->entry_speed.value <= allow_speed.value) {
			return;
		}
		backwardIterate(blockPtr);
	}
}

void ToolPathTimeEstimator::forwardIterate(Block *blockPtr) {

	blockPtr->exit_speed = allow_speed;
	blockPtr = blockPtr->next;
	blockPtr->entry_speed = allow_speed;
	allow_speed = Vector2D::getVectorFromDirection(sqrt(2 * acceleration.maxAlongDirection(blockPtr->length)*blockPtr->length.value + pow(blockPtr->entry_speed.value, 2)), blockPtr->length);
	if (allow_speed.value > Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length).value) {
		allow_speed = Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(blockPtr->length), blockPtr->length);
	}
	if (blockPtr->exit_speed.value <= allow_speed.value) {
		return;
	}
	forwardIterate(blockPtr);

}

Vector2D ToolPathTimeEstimator::calcAllowedSpeed(double acceleration, double initialSpeed, Vector2D length)
{
	return Vector2D::getVectorFromDirection(sqrt(2 * acceleration*length.value + initialSpeed*initialSpeed), length);
}



Vector2D ToolPathTimeEstimator::calcJunctionSpeed(Block *A, Block *B) {
	if (A == nullptr) {
		A = new Block(0, 0);
	}
	Vector2D A_speed, junctionSpeed;
	double arc_cos, arc_cos_2;
	A_speed = calcAllowedSpeed(acceleration.maxAlongDirection(A->length), A->entry_speed.value, A->length);
	arc_cos = (A_speed*B->length) / (A_speed.value*B->length.value);
	if (arc_cos < -1) {
	arc_cos = -1;
}
	if (abs(arc_cos - 1) <= 1e-10) {
		if (A_speed.value >= nominal_speed.maxAlongDirection(A->length)) {
			return Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(A->length), A->length);
		}
		return A_speed;
	}
	arc_cos_2 = sqrt((1 + arc_cos) / 2);
	junctionSpeed = Vector2D::getVectorFromDirection(jerk.maxAlongDirection(B->length.normalize()-A->length.normalize()) / (2 * sqrt(1 - pow(arc_cos_2, 2))), A->length);
	if (junctionSpeed.value >= nominal_speed.maxAlongDirection(A->length)) {
		return Vector2D::getVectorFromDirection(nominal_speed.maxAlongDirection(A->length), A->length);
	}
	return junctionSpeed;
}





void ToolPathTimeEstimator::detail(FILE * out) {
	//Sum up the build time
	double buildTime = 0;
	double length = 0;
	Block * blockPtr;
	blockPtr = this->blocks.blockHead->next;
	while (blockPtr->next != nullptr) {
		blockPtr->build_time = calcBlockBuildTime(blockPtr);
		if (abs(blockPtr->entry_speed.value - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length) >= blockPtr->build_time + 1e-4) {
			//can't be small,only equal, this is the minimum time, add a epsilon in case the error
			fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, blockPtr->entry_speed.value, blockPtr->exit_speed.value, blockPtr->length.value, blockPtr->build_time);
		}
		else {
			if (((nominal_speed.maxAlongDirection(blockPtr->length) - blockPtr->entry_speed.value) + (nominal_speed.maxAlongDirection(blockPtr->length) - blockPtr->exit_speed.value)) /
				acceleration.maxAlongDirection(blockPtr->length) >= blockPtr->build_time + 1e-4) {
				//can not reach the nominal speed or just reach the nominal speed in a very short time
				if (blockPtr->entry_speed.value - blockPtr->exit_speed.value <= 1e-4) {
					double mid_x, mid_y, mid_speed;
					mid_x = (blockPtr->x + blockPtr->previous->x) / 2;
					mid_y = (blockPtr->y + blockPtr->previous->y) / 2;
					mid_speed = (blockPtr->entry_speed.value + blockPtr->exit_speed.value) / 2;
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid_x, mid_y, blockPtr->entry_speed.value, mid_speed, blockPtr->length.value / 2, blockPtr->build_time / 2);
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value, blockPtr->length.value / 2, blockPtr->build_time / 2);
				}
				else {
					if (blockPtr->entry_speed.value > blockPtr->exit_speed.value) {
						double mid_x, mid_y, mid_speed, mid_time;
						Vector2D travel;
						mid_speed = (blockPtr->build_time*acceleration.maxAlongDirection(blockPtr->length) + blockPtr->entry_speed.value + blockPtr->exit_speed.value) / 2;
						travel = Vector2D::getVectorFromDirection(distance(acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, mid_speed), blockPtr->length);
						mid_x = blockPtr->previous->x + travel.x;
						mid_y = blockPtr->previous->y + travel.y;
						mid_time = (mid_speed - blockPtr->entry_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid_x, mid_y, blockPtr->entry_speed.value, mid_speed, travel.value, mid_time);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value, blockPtr->length.value - travel.value, blockPtr->build_time - mid_time);
					}
					else {
						double mid_x, mid_y, mid_speed, mid_time;
						Vector2D travel;
						mid_speed = (blockPtr->build_time*acceleration.maxAlongDirection(blockPtr->length) + blockPtr->entry_speed.value + blockPtr->exit_speed.value) / 2;
						travel = Vector2D::getVectorFromDirection(distance(acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, mid_speed), blockPtr->length);
						mid_x = blockPtr->previous->x + travel.x;
						mid_y = blockPtr->previous->y + travel.y;
						mid_time = (mid_speed - blockPtr->entry_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid_x, mid_y, blockPtr->entry_speed.value, mid_speed, travel.value, mid_time);
						fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value, blockPtr->length.value - travel.value, blockPtr->build_time - mid_time);
					}
				}
			}
			else {
				//Must reach the maximum speed and move with this speed
				double mid1_x, mid1_y, mid_speed, mid1_time;
				double mid2_x, mid2_y, mid2_time_from_exit;
				Vector2D travel_acc, travel_curise;
				mid_speed = nominal_speed.maxAlongDirection(blockPtr->length);
				travel_acc = Vector2D::getVectorFromDirection(distance(acceleration.maxAlongDirection(blockPtr->length), blockPtr->entry_speed.value, mid_speed), blockPtr->length);
				mid1_x = blockPtr->previous->x + travel_acc.x;
				mid1_y = blockPtr->previous->y + travel_acc.y;
				mid1_time = (mid_speed - blockPtr->entry_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
				mid2_time_from_exit = (mid_speed - blockPtr->exit_speed.value) / acceleration.maxAlongDirection(blockPtr->length);
				travel_curise = Vector2D::getVectorFromDirection(mid_speed*(blockPtr->build_time - mid1_time - mid2_time_from_exit), blockPtr->length);
				mid2_x = mid1_x + travel_curise.x;
				mid2_y = mid1_y + travel_curise.y;
				if (blockPtr->entry_speed.value - mid_speed > 1e-5) {
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid1_x, mid1_y, blockPtr->entry_speed.value, mid_speed, travel_acc.value, mid1_time);
				}
				
				fprintf(out, "%f, %f, %f, %f, %f, %f\n", mid2_x, mid2_y, blockPtr->entry_speed.value, mid_speed, travel_curise.value,
					blockPtr->build_time - mid1_time - mid2_time_from_exit);
				if (mid_speed - blockPtr->exit_speed.value > 1e-5) {
					fprintf(out, "%f, %f, %f, %f, %f, %f\n", blockPtr->x, blockPtr->y, mid_speed, blockPtr->exit_speed.value,
						(blockPtr->exit_speed.value + mid_speed)*(blockPtr->build_time - mid1_time - mid2_time_from_exit) / 2, mid2_time_from_exit);
				}
			}
		}
			buildTime += blockPtr->build_time;
			length += blockPtr->length.value;
			blockPtr = blockPtr->next;
	}
}

double ToolPathTimeEstimator::distance(double acceleration, double small_speed, double large_speed) {
	return (pow(large_speed, 2) - pow(small_speed, 2)) / (2 * acceleration);
}


