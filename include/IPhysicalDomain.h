#ifndef IPHYSICAL_DOMAIN_H
#define IPHYSICAL_DOMAIN_H

#include <deal.II/grid/tria.h>
#include <memory>
#include <string>

template <int dim>
class IPhysicalDomain
{
public:
  virtual const dealii::Triangulation<dim - 1, dim> &getTria() const                          = 0;
  virtual dealii::Triangulation<dim - 1, dim> &      getTria()                                = 0;
  virtual std::vector<unsigned int>                  get_dirichlet_boundary_ids()             = 0;
  virtual std::vector<unsigned int>                  get_neumann_boundary_ids()               = 0;
  virtual void                                       read_domain(std::string input_path = "") = 0;
  virtual void refine_and_resize(std::string input_path = "")                                 = 0;
  virtual void update_triangulation()                                                         = 0;

  virtual void update_domain(double, double, double){};
};

#endif // IPHYSICAL_DOMAIN_H
