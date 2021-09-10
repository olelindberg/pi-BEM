#include <deal.II/dofs/dof_handler.h>

class dof_handler_util
{
public:
  template <int dim>
  static int
  number_of_materials(const dealii::DoFHandler<dim - 1, dim> &dh)
  {
    int numMat = 0;
    for (auto cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      numMat = std::max(numMat, (int)(cell->material_id() + 1));

      for (int j = 0; j < 4; ++j)
      {
        auto line = cell->line(j);
        numMat    = std::max(numMat, (int)(line->manifold_id() + 1));
      }
    }
    return numMat;
  }
};
