#include "stdafx.h"
#include "Block.hpp"


Block::Block(double x, double y) {
    this->next = nullptr;
    this->previous = nullptr;
    this->x = x;
    this->y = y;
    this->length = Vector2D(0, 0);
}