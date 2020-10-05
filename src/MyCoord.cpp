/***************************************************************************
 * file:        MyCoord.cc
 *
 * author:      Christian Frank
 *
 * copyright:   (C) 2004 Telecommunication Networks Group (TKN) at
 *              Technische Universitaet Berlin, Germany.
 *
 *              This program is free software; you can redistribute it
 *              and/or modify it under the terms of the GNU General Public
 *              License as published by the Free Software Foundation; either
 *              version 2 of the License, or (at your option) any later
 *              version.
 *              For further information see file COPYING
 *              in the top level directory
 ***************************************************************************
 * part of:     framework implementation developed by tkn
 **************************************************************************/

#include "assert.h"

#include "MyCoord.h"


const MyCoord MyCoord::NIL = MyCoord(NaN, NaN, NaN);
const MyCoord MyCoord::ZERO = MyCoord(0.0, 0.0, 0.0);
const MyCoord MyCoord::X_AXIS = MyCoord(1.0, 0.0, 0.0);
const MyCoord MyCoord::Y_AXIS = MyCoord(0.0, 1.0, 0.0);
const MyCoord MyCoord::Z_AXIS = MyCoord(0.0, 0.0, 1.0);

/**
 * On a torus the end and the begin of the axes are connected so you
 * get a circle. On a circle the distance between two points can't be greater
 * than half of the circumference.
 * If the normal distance between two points on one axis is bigger than
 * half of the size there must be a "shorter way" over the border on this axis
 */
static double dist(double coord1, double coord2, double size)
{
    double difference = fabs(coord1 - coord2);
    if (difference == 0)
        // NOTE: event if size is zero
        return 0;
    else {
        assert(size != 0);
        double dist = MyCoord::modulo(difference, size);
        return std::min(dist, size - dist);
    }
}

double MyCoord::sqrTorusDist(const MyCoord& b, const MyCoord& size) const
{
    double xDist = dist(x, b.x, size.x);
    double yDist = dist(y, b.y, size.y);
    double zDist = dist(z, b.z, size.z);
    return xDist * xDist + yDist * yDist + zDist * zDist;
}

