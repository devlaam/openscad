#include "geometry/boolean_utils.h"

#include <utility>
#include <memory>
#include <vector>

#ifdef ENABLE_CGAL
#include "geometry/cgal/CGALNefGeometry.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include "geometry/cgal/cgalutils.h"
#endif  // ENABLE_CGAL
#ifdef ENABLE_MANIFOLD
#include "geometry/manifold/ManifoldGeometry.h"
#include "geometry/manifold/manifoldutils.h"
#endif  // ENABLE_MANIFOLD

#include "glview/RenderSettings.h"
#include "geometry/PolySet.h"
#include "utils/printutils.h"

#include "geometry/Reindexer.h"
#include "geometry/GeometryUtils.h"
//RUUD
#include "geometry/linalg.h"

#ifdef ENABLE_CGAL
std::unique_ptr<PolySet> applyHull(const Geometry::Geometries& children)
{
  using Hull_kernel = CGAL::Epick;
  // Collect point cloud
  Reindexer<Hull_kernel::Point_3> reindexer;

  auto addCapacity = [&](const auto n) { reindexer.reserve(reindexer.size() + n); };

  auto addPoint = [&](const auto& v) { reindexer.lookup(v); };

  for (const auto& item : children) {
    auto& chgeom = item.second;
#ifdef ENABLE_CGAL
    if (const auto *N = dynamic_cast<const CGALNefGeometry *>(chgeom.get())) {
      if (!N->isEmpty()) {
        addCapacity(N->p3->number_of_vertices());
        for (auto it = N->p3->vertices_begin(); it != N->p3->vertices_end(); ++it) {
          addPoint(CGALUtils::vector_convert<Hull_kernel::Point_3>(it->point()));
        }
      }
#endif  // ENABLE_CGAL
#ifdef ENABLE_MANIFOLD
    } else if (const auto *mani = dynamic_cast<const ManifoldGeometry *>(chgeom.get())) {
      addCapacity(mani->numVertices());
      mani->foreachVertexUntilTrue([&](auto& p) {
        addPoint(CGALUtils::vector_convert<Hull_kernel::Point_3>(p));
        return false;
      });
#endif  // ENABLE_MANIFOLD
    } else if (const auto *ps = dynamic_cast<const PolySet *>(chgeom.get())) {
      addCapacity(ps->indices.size() * 3);
      for (const auto& p : ps->indices) {
        for (const auto& ind : p) {
          addPoint(CGALUtils::vector_convert<Hull_kernel::Point_3>(ps->vertices[ind]));
        }
      }
    }
  }

  const auto& points = reindexer.getArray();
  if (points.size() <= 3) return nullptr;

  // Apply hull
  if (points.size() >= 4) {
    try {
      CGAL::Polyhedron_3<Hull_kernel> r;
      CGAL::convex_hull_3(points.begin(), points.end(), r);
      PRINTDB("After hull vertices: %d", r.size_of_vertices());
      PRINTDB("After hull facets: %d", r.size_of_facets());
      PRINTDB("After hull closed: %d", r.is_closed());
      PRINTDB("After hull valid: %d", r.is_valid());
      // FIXME: Make sure PolySet is set to convex.
      // FIXME: Can we guarantee a manifold PolySet here?
      return CGALUtils::createPolySetFromPolyhedron(r);
    } catch (const CGAL::Failure_exception& e) {
      LOG(message_group::Error, "CGAL error in applyHull(): %1$s", e.what());
    }
  }
  return nullptr;
}

//RUUD(box-3D)
std::unique_ptr<PolySet> applyBox(const Geometry::Geometries& children, const Vector3d& add)
{
  using Box_kernel = CGAL::Epick;
  const double eps = 1e-9;

  // Track min/max; start with empty state
  bool hasPoint = false;
  bool resizefault = false;
  double minx = 0.0, miny = 0.0, minz = 0.0;
  double maxx = 0.0, maxy = 0.0, maxz = 0.0;

  auto updateBounds = [&](const Box_kernel::Point_3& p) {
    const double x = CGAL::to_double(p.x());
    const double y = CGAL::to_double(p.y());
    const double z = CGAL::to_double(p.z());
    if (!hasPoint) {
      minx = maxx = x; miny = maxy = y; minz = maxz = z;
      hasPoint = true;
    } else {
      if (x < minx) minx = x; if (x > maxx) maxx = x;
      if (y < miny) miny = y; if (y > maxy) maxy = y;
      if (z < minz) minz = z; if (z > maxz) maxz = z;
    }
  };

  // Collect bounds from all child geometry types
  for (const auto& item : children) {
    auto& chgeom = item.second;
#ifdef ENABLE_CGAL
    if (const auto *N = dynamic_cast<const CGALNefGeometry *>(chgeom.get())) {
      if (!N->isEmpty()) {
        for (auto it = N->p3->vertices_begin(); it != N->p3->vertices_end(); ++it) {
          updateBounds(CGALUtils::vector_convert<Box_kernel::Point_3>(it->point()));
        }
      }
#endif  // ENABLE_CGAL
#ifdef ENABLE_MANIFOLD
    } else if (const auto *mani = dynamic_cast<const ManifoldGeometry *>(chgeom.get())) {
      mani->foreachVertexUntilTrue([&](auto& v) {
        updateBounds(CGALUtils::vector_convert<Box_kernel::Point_3>(v));
        return false;
      });
#endif  // ENABLE_MANIFOLD
    } else if (const auto *ps = dynamic_cast<const PolySet *>(chgeom.get())) {
      // Iterate faces -> indices -> vertices, same as in applyHull()
      for (const auto& face : ps->indices) {
        for (const auto& idx : face) {
          updateBounds(CGALUtils::vector_convert<Box_kernel::Point_3>(ps->vertices[idx]));
        }
      }
    }
  }

  if (!hasPoint) return nullptr;

  // Build a CGAL hexahedron for the AABB and convert to PolySet
  auto make_box_polyset = [&](double ax, double ay, double az,
                              double bx, double by, double bz) -> std::unique_ptr<PolySet>
  {
    const Box_kernel::Point_3 p000(ax, ay, az);
    const Box_kernel::Point_3 p001(ax, ay, bz);
    const Box_kernel::Point_3 p010(ax, by, az);
    const Box_kernel::Point_3 p011(ax, by, bz);
    const Box_kernel::Point_3 p100(bx, ay, az);
    const Box_kernel::Point_3 p101(bx, ay, bz);
    const Box_kernel::Point_3 p110(bx, by, az);
    const Box_kernel::Point_3 p111(bx, by, bz);

    CGAL::Polyhedron_3<Box_kernel> poly;
    // For the order see: https://doc.cgal.org/latest/BGL/group__PkgBGLGeneratorFct.html#ga12fa3e202c24740dade5764e3ea80c41
    CGAL::make_hexahedron(p000, p100, p110, p010, p011, p001, p101, p111, poly);

    PRINTDB("Box closed: %d", poly.is_closed());
    PRINTDB("Box valid: %d", poly.is_valid());

    return CGALUtils::createPolySetFromPolyhedron(poly);
  };

  // See if we can expand/contract the bounding box. If not, ignore.
  if ( (maxx - minx + 2*add[0]) > eps ) { minx -= add[0]; maxx += add[0] } else { resizefault = true; }
  if ( (maxy - miny + 2*add[1]) > eps ) { miny -= add[1]; maxy += add[1] } else { resizefault = true; }
  if ( (maxz - minz + 2*add[2]) > eps ) { minz -= add[2]; maxz += add[2] } else { resizefault = true; }

  if (resizefault) { LOG("WARNING: Bounding box has one or more negative dimensions, ignoring offending size modifications."); }

  try {
    return make_box_polyset(minx, miny, minz, maxx, maxy, maxz);
  } catch (const CGAL::Failure_exception& e) {
    LOG(message_group::Warning, "CGAL error in applyBox() (degenerate Box?): %1$s", e.what());
  }

  return nullptr;
}

//END-RUUD


/*!
   children cannot contain nullptr objects

  FIXME: This shouldn't return const, but it does due to internal implementation details
 */
std::shared_ptr<const Geometry> applyMinkowski(const Geometry::Geometries& children)
{
#if ENABLE_MANIFOLD
  if (RenderSettings::inst()->backend3D == RenderBackend3D::ManifoldBackend) {
    return ManifoldUtils::applyMinkowski(children);
  }
#endif  // ENABLE_MANIFOLD
  return CGALUtils::applyMinkowski3D(children);
}
#else   // ENABLE_CGAL
std::unique_ptr<PolySet> applyHull(const Geometry::Geometries& children)
{
  return std::make_unique<PolySet>(3, true);
}

//RUUD
std::unique_ptr<PolySet> applyBox(const Geometry::Geometries& children, const Vector3d& add)
{
  return std::make_unique<PolySet>(3, true);
}

std::shared_ptr<const Geometry> applyMinkowski(const Geometry::Geometries& children)
{
  return std::make_shared<PolySet>(3);
}
#endif  // !ENABLE_CGAL
