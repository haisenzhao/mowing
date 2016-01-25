#include "stdafx.h"
#include "Vector2D.hpp"

class Block
{
public:
	Block(double x, double y);
	~Block();
    double x,y;
	float s1, s2, s3;
	bool jump;
    Vector2D max_speed;
	Vector2D nominal_speed;
	Vector2D entry_speed;
	Vector2D exit_speed;
	double build_time;
	Vector2D length;
    Vector2D direction;
    Block* next;
    Block* previous;
};