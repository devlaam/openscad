#pragma once

#include <memory>
#include "geometry/PolySet.h"
#include "geometry/Geometry.h"
//RUUD
#include "geometry/linalg.h"

std::unique_ptr<PolySet> applyHull(const Geometry::Geometries& children);
std::shared_ptr<const Geometry> applyMinkowski(const Geometry::Geometries& children);
//RUUD
std::unique_ptr<PolySet> applyBox(const Geometry::Geometries& children, const Vector3d& add);
