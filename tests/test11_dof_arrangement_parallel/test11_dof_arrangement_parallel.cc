#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_c1.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <boost/archive/text_oarchive.hpp>

#include <fstream>
#include <iostream>


int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, dealii::numbers::invalid_unsigned_int);
  auto np  = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto pid = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::cout << pid << "/" << np << std::endl;

  int degree = 3;
  if (argc > 1)
    degree = std::stoi(argv[1]);

  int subdivisions = 2;
  if (argc > 2)
    subdivisions = std::stoi(argv[2]);


  const int                     dim = 3;
  dealii::Triangulation<2, dim> tria;

  dealii::GridGenerator::subdivided_hyper_cube(tria, subdivisions);

  // std::string grid_filename =
  //   "/home/ole/dev/projects/pi-BEM/docs/data/DTC_MASHCON_2019/C3/mesh.inp";
  // std::ifstream in;
  // in.open(grid_filename);
  // dealii::GridIn<dim - 1, dim> gridIn;
  // gridIn.attach_triangulation(tria);
  // gridIn.read_ucd(in, true);


  dealii::GridTools::partition_triangulation(np, tria);

  std::ofstream   out("/home/ole/dev/temp/grid.inp");
  dealii::GridOut grid_out;
  grid_out.write_ucd(tria, out);

  dealii::FE_Q<2, dim>     fe_scalar(degree);
  dealii::FESystem<2, dim> fe_vector(fe_scalar, 3);

  dealii::DoFHandler dh_scalar(tria);
  dh_scalar.distribute_dofs(fe_scalar);

  dealii::DoFHandler dh_vector(tria);
  dh_vector.distribute_dofs(fe_vector);

  std::vector<dealii::types::subdomain_id> dofs_domain_association(dh_scalar.n_dofs());
  dealii::DoFTools::get_subdomain_association(dh_scalar, dofs_domain_association);
  dealii::IndexSet this_cpu_set(dh_scalar.n_dofs());
  for (dealii::types::global_dof_index i = 0; i < dh_scalar.n_dofs(); ++i)
    if (dofs_domain_association[i] == pid)
      this_cpu_set.add_index(i);
  this_cpu_set.compress();

  std::vector<dealii::types::subdomain_id> dofs_domain_association_vector(dh_vector.n_dofs());
  dealii::DoFTools::get_subdomain_association(dh_vector, dofs_domain_association_vector);
  dealii::IndexSet this_cpu_set_vector(dh_vector.n_dofs());
  for (dealii::types::global_dof_index i = 0; i < dh_vector.n_dofs(); ++i)
    if (dofs_domain_association_vector[i] == pid)
      this_cpu_set_vector.add_index(i);
  this_cpu_set_vector.compress();


  std::vector<dealii::types::global_dof_index> local_dofs_scalar(fe_scalar.dofs_per_cell);
  std::vector<dealii::types::global_dof_index> local_dofs_vector(fe_vector.dofs_per_cell);

  dealii::TrilinosWrappers::MPI::Vector scalar_solution(this_cpu_set, MPI_COMM_WORLD);
  dealii::TrilinosWrappers::MPI::Vector vector_solution(this_cpu_set_vector, MPI_COMM_WORLD);

  int  mapping_degree = 3;
  auto mapping        = std::make_shared<dealii::MappingQGeneric<dim - 1, dim>>(mapping_degree);
  std::vector<dealii::Point<dim>> support_points(dh_scalar.n_dofs());
  dealii::DoFTools::map_dofs_to_support_points(*mapping, dh_scalar, support_points);

  auto cell_vector = dh_vector.begin_active();
  for (const auto &cell_scalar : dh_scalar.active_cell_iterators())
  {
    if (cell_scalar->subdomain_id() == pid)
    {
      cell_scalar->get_dof_indices(local_dofs_scalar);
      cell_vector->get_dof_indices(local_dofs_vector);

      for (std::size_t i = 0; i < local_dofs_scalar.size(); ++i)
      {
        auto pnt = support_points[local_dofs_scalar[i]];

        for (std::size_t j = 0; j < 3; ++j)
          vector_solution[local_dofs_vector[i * 3 + j]] = pnt[j];

        scalar_solution[local_dofs_scalar[i]] = pnt[0];
      }
    }
    ++cell_vector;
  }

  const dealii::Vector<double> scalar_solution_local(scalar_solution);

  if (pid == 0)
  {
    dealii::DataOut<dim - 1, dealii::DoFHandler<dim - 1, dim>> dataout_scalar;
    dataout_scalar.attach_dof_handler(dh_scalar);

    dataout_scalar.add_data_vector(
      scalar_solution_local,
      "scalar",
      dealii::DataOut<dim - 1, dealii::DoFHandler<dim - 1, dim>>::type_dof_data);

    dataout_scalar.build_patches(
      *mapping,
      mapping_degree,
      dealii::DataOut<dim - 1, dealii::DoFHandler<dim - 1, dim>>::curved_inner_cells);

    std::string   filename_scalar = "/home/ole/dev/temp/scalar.vtu";
    std::ofstream file_scalar(filename_scalar);
    dataout_scalar.write_vtu(file_scalar);
  }

  const dealii::Vector<double> vector_solution_local(vector_solution);
  if (pid == 0)
  {
    dealii::DataOut<dim - 1, dealii::DoFHandler<dim - 1, dim>> dataout_vector;
    dataout_vector.attach_dof_handler(dh_vector);

    std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, dealii::DataComponentInterpretation::component_is_part_of_vector);
    dataout_vector.add_data_vector(
      vector_solution_local,
      std::vector<std::string>(dim, "vector"),
      dealii::DataOut<dim - 1, dealii::DoFHandler<dim - 1, dim>>::type_dof_data,
      data_component_interpretation);

    dataout_vector.build_patches(
      *mapping,
      mapping_degree,
      dealii::DataOut<dim - 1, dealii::DoFHandler<dim - 1, dim>>::curved_inner_cells);

    std::string   filename_vector = "/home/ole/dev/temp/vector.vtu";
    std::ofstream file_vector(filename_vector);
    dataout_vector.write_vtu(file_vector);
  }
}