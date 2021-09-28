#include "../include/CurvatureRefinement.h"
#include "my_utilities.h"

void
CurvatureRefinement::refine(dealii::Triangulation<2, 3> &    tria,
                            const std::vector<TopoDS_Shape> &cad_surfaces)
{
  // the following refinement cycle is based upon the original CAD
  // geometry curvature. For this reason it can be activated not only
  // when the user requires it with the surface_curvature_refinement
  // option in the input file. Of course, this is possible only if
  // also the use_cad_surface_and_curves flag is set to thrue through
  // the input file. Only in this way in fact, CAD surfaces and curves
  // are prescribed for the triangulation refinements on some of
  // its manifold ids.
  if (use_cad_surface_and_curves && surface_curvature_refinement)
  {
    _pcout << "Refining based on surface curvature ...\n";
    const double tolerance          = cad_to_projectors_tolerance_ratio * max_tol;
    unsigned int refinedCellCounter = 1;
    unsigned int cycles_counter     = 0;
    // the refinement procedure is recursively repeated until no more cells
    // are flagged for refinement, or until the user specified maximum number
    // of curvature refinement cycles is reached
    while ((refinedCellCounter) && (cycles_counter < max_curvature_ref_cycles))
    {
      std::vector<double> cell_size_all;

      // the refined cells counter is zeroed at the start of each cycle
      refinedCellCounter = 0;
      // we loop on the all the triangulation active cells
      dealii::Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
      dealii::Triangulation<2, 3>::active_cell_iterator endc = tria.end();
      for (; cell != endc; ++cell)
      {
        double cell_size;
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
          dealii::Point<3> projection;
          {
            projection = dealii::my_line_intersection<3>(neededShape, cell->center(), n, tolerance);
          }
          double min_curvature = 0.0;
          double max_curvature = 0.0;
          {
            dealii::OpenCASCADE::surface_curvature(
              neededShape, projection, tolerance, min_curvature, max_curvature);
          }
          // among the differential point, we select the maximum
          // absolute curvature
          double max_abs_curv = fmax(fabs(min_curvature), fabs(max_curvature));

          // radius is computed from the maximum absolute curvatur
          double curvature_radius = 1.0 / fmax(max_abs_curv, tolerance);

          // the target cell size is selected so that it corresponds to
          // a cells_per_circle fraction of the circumference
          // corresponding to the minimum curvature radius
          cell_size = 2.0 * dealii::numbers::PI / cells_per_circle * curvature_radius;
        }
        else
        {
          // if the cell manifold_id is not associated to a CAD
          // surface, the target cell_size is set to and extremely high
          // value, so that the cell is never refined
          cell_size = 2 * dealii::numbers::PI / cells_per_circle / tolerance;
        }
        cell_size_all.push_back(cell_size);
        // the following line si for debug puropses and should be
        // uncommented if something is not working with the refinement

        // if the cell diameter is higher than the target cell size, the
        // refinement flag is set (unless the cell is already very small
        // ---which for us means 10xtolerance)
        if ((cell->diameter() > cell_size) && (cell->diameter() > 10 * tolerance))
        {
          cell->set_refine_flag();
          refinedCellCounter++;
        }
      }
      double min_cell_size = *std::min_element(cell_size_all.begin(), cell_size_all.end());
      double max_cell_size = *std::max_element(cell_size_all.begin(), cell_size_all.end());

      _pcout << "min_cell_size: " << min_cell_size << std::endl;
      _pcout << "max_cell_size: " << max_cell_size << std::endl;

      // the number of cells to be refined in this cycle is reported, the
      // refinement is carried out and the make_edges_conformal function is
      // called to check no edge presents non comformities
      _pcout << "Curvature Based Local Refinement Cycle: " << cycles_counter << " ("
             << refinedCellCounter << ")" << std::endl;
      tria.execute_coarsening_and_refinement();
      //      make_edges_conformal(_withDoubleNodes);
      cycles_counter++;
    }
  }
  _pcout << "We have a tria of " << tria.n_active_cells() << " cells." << std::endl;
}