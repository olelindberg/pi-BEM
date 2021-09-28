#include "../include/HeightRatioRefinement.h"
#include "my_utilities.h"
void
HeightRatioRefinement::refine(dealii::Triangulation<2, 3> &    tria,
                              const std::vector<TopoDS_Shape> &cad_surfaces)
{
  _pcout << "Height ratio refinement ...\n";

  for (int iter = 0; iter < itermax; ++iter)
  {
    bool isrefining = false;
    for (dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
         cell != tria.end();
         ++cell)
    {
      double height_ratio = 0.0;
      if (int(cell->material_id()) - 1 < (int)cad_surfaces.size())
      {
        // In the following lines, we try to come up with an estimation
        // of the cell normal. It is obtained from the average of the
        // normal to the 4 triangles in which the cell can be split using
        // the vertices and the center. The commented lines can be used
        // for checks in case something goes wrong.

        dealii::Point<3> t0 = cell->vertex(0) + (-1.0) * cell->center();
        dealii::Point<3> t1 = cell->vertex(1) + (-1.0) * cell->center();
        dealii::Point<3> t2 = cell->vertex(2) + (-1.0) * cell->center();
        dealii::Point<3> t3 = cell->vertex(3) + (-1.0) * cell->center();

        dealii::Point<3> nn0(t0(1) * t1(2) - t0(2) * t1(1),
                             t0(2) * t1(0) - t0(0) * t1(2),
                             t0(0) * t1(1) - t0(1) * t1(0));
        nn0 /= nn0.norm();
        dealii::Point<3> nn1(t1(1) * t3(2) - t1(2) * t3(1),
                             t1(2) * t3(0) - t1(0) * t3(2),
                             t1(0) * t3(1) - t1(1) * t3(0));
        nn1 /= nn1.norm();
        dealii::Point<3> nn2(t3(1) * t2(2) - t3(2) * t2(1),
                             t3(2) * t2(0) - t3(0) * t2(2),
                             t3(0) * t2(1) - t3(1) * t2(0));
        nn2 /= nn2.norm();
        dealii::Point<3> nn3(t2(1) * t0(2) - t2(2) * t0(1),
                             t2(2) * t0(0) - t2(0) * t0(2),
                             t2(0) * t0(1) - t2(1) * t0(0));
        nn3 /= nn3.norm();
        dealii::Point<3> n = (nn0 + nn1 + nn2 + nn3) / 4.0;
        n /= n.norm();

        // once the cell normal has beed created, we want to use it as
        // the direction of the projection onto the CAD surface first
        // though, let's check that we are using a CAD surface for the
        // refinement of the manifold_id associated with the present cell

        // if so, the cad_surface associated with the present
        // manifold_id is identified...
        TopoDS_Shape neededShape = cad_surfaces[int(cell->material_id()) - 1];
        // ...and used to set up a line intersection to project the
        // cell center on the CAD surface along the direction
        // specified by the previously computed cell normal
        dealii::Point<3> projection =
          dealii::my_line_intersection<3>(neededShape, cell->center(), n, tolerance);

        t0 = cell->vertex(0) + (-1.0) * projection;
        t1 = cell->vertex(1) + (-1.0) * projection;
        t2 = cell->vertex(2) + (-1.0) * projection;
        t3 = cell->vertex(3) + (-1.0) * projection;

        double proj0 = std::fabs(t0[0] * n[0] + t0[1] * n[1] + t0[2] * n[2]);
        double proj1 = std::fabs(t1[0] * n[0] + t1[1] * n[1] + t1[2] * n[2]);
        double proj2 = std::fabs(t2[0] * n[0] + t2[1] * n[1] + t2[2] * n[2]);
        double proj3 = std::fabs(t3[0] * n[0] + t3[1] * n[1] + t3[2] * n[2]);

        double dist_max = 0.0;
        dist_max        = std::max(proj0, dist_max);
        dist_max        = std::max(proj1, dist_max);
        dist_max        = std::max(proj2, dist_max);
        dist_max        = std::max(proj3, dist_max);

        double cell_extent = std::min(cell->extent_in_direction(0), cell->extent_in_direction(1));
        height_ratio       = dist_max / cell_extent;
      }

      if (height_ratio > height_ratio_max)
      {
        cell->set_refine_flag();
        isrefining = true;
      }
    }

    if (isrefining)
    {
      tria.prepare_coarsening_and_refinement();
      tria.execute_coarsening_and_refinement();
    }
    else
      break;
  }
}
