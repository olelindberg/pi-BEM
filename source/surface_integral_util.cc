

#include "../include/surface_integral_util.h"

#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/error_estimator.h>

#include <iomanip>
#include <iostream>
#include <limits>

#include "../include/dof_handler_util.h"
#include "../include/laplace_kernel.h"
#include "Teuchos_TimeMonitor.hpp"

template <>
surface_integral_util<3>::surface_integral_util(ComputationalDomain<3> &comp_dom, MPI_Comm comm)
  : comp_dom(comp_dom)
  , parsed_fe("Scalar FE", "FE_Q(1)")
  , parsed_gradient_fe("Vector FE", "FESystem[FE_Q(1)^3]", "u,u,u", 3)
  , dh(comp_dom.getTria())
  , gradient_dh(comp_dom.getTria())
  , _mpi_communicator(comm)
{}

template <int dim>
void surface_integral_util<dim>::reinit()
{
  fe          = parsed_fe();
  gradient_fe = parsed_gradient_fe();

  dh.distribute_dofs(*fe);
  gradient_dh.distribute_dofs(*gradient_fe);

  DoFRenumbering::component_wise(dh);
  DoFRenumbering::component_wise(gradient_dh);

  DoFRenumbering::subdomain_wise(dh);
  DoFRenumbering::subdomain_wise(gradient_dh);

  if (mapping_type == "FE")
  {
    map_vector.reinit(gradient_dh.n_dofs());
    VectorTools::get_position_vector(gradient_dh, map_vector);
  }

  if (!mapping)
  {
    if (mapping_type == "FE")
      mapping = std::make_shared<MappingFEField<dim - 1, dim>>(gradient_dh, map_vector);
    else
      mapping = std::make_shared<MappingQ<dim - 1, dim>>(mapping_degree);
  }

  const types::global_dof_index n_dofs = dh.n_dofs();

  std::vector<types::subdomain_id> dofs_domain_association(n_dofs);

  DoFTools::get_subdomain_association(dh, dofs_domain_association);
  std::vector<types::subdomain_id> vector_dofs_domain_association(gradient_dh.n_dofs());

  DoFTools::get_subdomain_association(gradient_dh, vector_dofs_domain_association);

  this_cpu_set.clear();
  this_cpu_set.set_size(n_dofs);

  for (types::global_dof_index i = 0; i < n_dofs; ++i)
    if (dofs_domain_association[i] == Utilities::MPI::this_mpi_process(_mpi_communicator))
    {
      this_cpu_set.add_index(i);
      for (unsigned int idim = 0; idim < dim; ++idim)
      {}
    }

  this_cpu_set.compress();

  IndexSet trial_index_set;
  trial_index_set.clear();
  trial_index_set = DoFTools::dof_indices_with_subdomain_association(
    gradient_dh, Utilities::MPI::this_mpi_process(_mpi_communicator));
}


template <int dim>
void surface_integral_util<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Quadrature rules");
  {
    prm.declare_entry("Quadrature type",
                      "gauss",
                      Patterns::Selection(QuadratureSelector<(dim - 1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
  }
  prm.leave_subsection();

  prm.declare_entry("Mapping Type", "FE", Patterns::Selection("FE|Q"));

  prm.declare_entry("Mapping Q Degree", "1", Patterns::Integer());
}

template <int dim>
void surface_integral_util<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Quadrature rules");
  {
    quadrature = std::shared_ptr<Quadrature<dim - 1>>(
      new QuadratureSelector<dim - 1>(prm.get("Quadrature type"),
                                      prm.get_integer("Quadrature order")));
    _quadrature1d = std::shared_ptr<Quadrature<dim - 2>>(
      new QuadratureSelector<dim - 2>(prm.get("Quadrature type"),
                                      prm.get_integer("Quadrature order")));
  }
  prm.leave_subsection();

  mapping_type   = prm.get("Mapping Type");
  mapping_degree = prm.get_integer("Mapping Q Degree");
}



template <int dim>
double surface_integral_util<dim>::ssurffint(const Body &body, int dimId)
{
  FEValues<dim - 1, dim> fe_v(*mapping,
                              dh.get_fe(),
                              *quadrature,
                              update_values | update_normal_vectors | update_quadrature_points |
                                update_JxW_values);

  FEFaceValues<dim - 1, dim> fe_face_values(*mapping,
                                            dh.get_fe(),
                                            *_quadrature1d,
                                            update_values | update_normal_vectors |
                                              update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = fe_v.n_quadrature_points;

  const types::global_dof_index n_dofs = dh.n_dofs();
  std::vector<Point<dim>>       support_points(n_dofs);
  DoFTools::map_dofs_to_support_points<dim - 1, dim>(*mapping, dh, support_points);
  std::vector<types::global_dof_index> local_dof_indices(fe->dofs_per_cell);

  double area = 0.0;
  for (const auto &cell : dh.active_cell_iterators())
  {
    if (cell->subdomain_id() == Utilities::MPI::this_mpi_process(_mpi_communicator))
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      for (int j = 0; j < 4; ++j)
      {
        auto line = cell->line(j);
        if (line->at_boundary() && body.isWaterline(line->manifold_id()))
        {
          fe_face_values.reinit(cell, j);

          auto qpoints  = fe_face_values.get_quadrature_points();
          auto qnormals = fe_face_values.get_normal_vectors();

          for (unsigned int q = 0; q < fe_face_values.n_quadrature_points; ++q)
          {
            for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i)
            {
              //              area += qpoints[q][dimId] * fe_face_values.shape_value(i, q) *
              //              fe_face_values.JxW(q);
            } // for i
          }   // for q
        }     // if line at bnd
      }       // for j faces
    }         // if this cpu
  }           // for cell in active cells

  area = Utilities::MPI::sum(area, _mpi_communicator);
  return area;
}



template class surface_integral_util<3>;
