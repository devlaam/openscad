#pragma once
#undef DIFFERENCE  // #defined in winuser.h

//RUUD
enum class OpenSCADOperator { UNION, INTERSECTION, DIFFERENCE, MINKOWSKI, HULL, FILL, RESIZE, BOX };
