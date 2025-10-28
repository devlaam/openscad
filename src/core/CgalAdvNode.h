#pragma once

#include <string>

#include "core/node.h"
#include "core/ModuleInstantiation.h"
#include "geometry/linalg.h"

//RUUD
enum class CgalAdvType { MINKOWSKI, HULL, FILL, RESIZE, BOX };

class CgalAdvNode : public AbstractNode
{
public:
  VISITABLE();
  CgalAdvNode(const ModuleInstantiation *mi, CgalAdvType type) : AbstractNode(mi), type(type) {}
  std::string toString() const override;
  std::string name() const override;

  unsigned int convexity{1};
  Vector3d newsize;
  Eigen::Matrix<bool, 3, 1> autosize;
  CgalAdvType type;
  //RUUD
  Vector3d add;
  bool act;

};
