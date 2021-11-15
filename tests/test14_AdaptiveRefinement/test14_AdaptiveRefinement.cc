//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_c1.h>
#include <deal.II/grid/grid_tools.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "./../tests.h"
#include "AdaptiveRefinement.h"
#include "bem_problem.h"
#include "boundary_conditions.h"
#include "computational_domain.h"



void
flowAroundSphere(double  x0,
                 double  y0,
                 double  z0,
                 double  x,
                 double  y,
                 double  z,
                 double  U,
                 double  r,
                 double &pot,
                 double &u,
                 double &v,
                 double &w)
{
  x = x - x0;
  y = y - y0;
  z = z - z0;

  double r3 = r * r * r;
  double x2 = x * x;
  double y2 = y * y;
  double z2 = z * z;

  double tmp1 = std::pow(x2 + y2 + z2, -1.5);
  double tmp2 = std::pow(x2 + y2 + z2, -2.5);
  pot         = 0.5 * U * r3 * x * tmp1;
  u           = -1.5 * U * r3 * x2 * tmp2 + 0.5 * U * r3 * tmp1;
  v           = -1.5 * U * r3 * x * y * tmp2;
  w           = -1.5 * U * r3 * x * z * tmp2;
}

class Options
{
public:
  Options(int argc, char *argv[])
  {
    _desc = std::make_shared<boost::program_options::options_description>(
      "test08_center_of_pressure options");
    _desc->add_options()("help", "produce help message");
    _desc->add_options()("input-file", boost::program_options::value<std::string>(), "Input file");


    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, *_desc),
                                  _vm);
    boost::program_options::notify(_vm);

    if (_vm.count("help"))
      std::cout << *_desc << "\n";
  }

  std::string
  get_input_file()
  {
    return _vm["input-file"].as<std::string>();
  };

private:
  std::shared_ptr<boost::program_options::options_description> _desc;
  boost::program_options::variables_map                        _vm;
};


int
main(int argc, char **argv)
{
  unsigned int np;
  if (argc == 1)
    np = numbers::invalid_unsigned_int;
  else
    np = atoi(argv[1]);
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, np);
  auto                             pid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  std::cout << "MPI initialized with " << np << " np" << std::endl;

  int numRefinements = 1;
  if (argc > 1)
    numRefinements = atoi(argv[2]);

  dealii::ConditionalOStream pcout(std::cout);
  double                     errorEstimatorMax = 1.0;
  double                     aspectRatioMax    = 2.5;

  //---------------------------------------------------------------------------
  // Initialize and read grid:
  //---------------------------------------------------------------------------
  dealii::Triangulation<2, 3> tria;
  std::string                 filename = "mesh.inp";
  std::ifstream               in;
  in.open(filename);
  dealii::GridIn<2, 3> gridIn;
  if (in.is_open())
  {
    gridIn.attach_triangulation(tria);
    gridIn.read_ucd(in, true);
  }
  else
  {
    std::cout << "Unable to open file\n";
  }


  for (int j = 0; j < numRefinements; ++j)
  {
    //---------------------------------------------------------------------------
    // Initialize scalar stuff:
    //---------------------------------------------------------------------------
    std::cout << pid << " Initialize scalar stuff 1\n";
    dealii::FE_Q<2, 3> fe_scalar(1);
    std::cout << pid << " Initialize scalar stuff 2\n";
    dealii::DoFHandler dh_scalar(tria);
    std::cout << pid << " Initialize scalar stuff 3\n";
    dh_scalar.distribute_dofs(fe_scalar);
    DoFRenumbering::component_wise(dh_scalar);
    std::cout << pid << " Initialize scalar stuff 4\n";
    std::vector<dealii::types::subdomain_id> dofs_domain_association(dh_scalar.n_dofs());
    std::cout << pid << " Initialize scalar stuff 5\n";
    dealii::DoFTools::get_subdomain_association(dh_scalar, dofs_domain_association);
    std::cout << pid << " Initialize scalar stuff 6\n";
    dealii::IndexSet this_cpu_set;
    this_cpu_set.set_size(dh_scalar.n_dofs());
    std::cout << pid << " Initialize scalar stuff 7\n";
    for (dealii::types::global_dof_index i = 0; i < dh_scalar.n_dofs(); ++i)
      if (dofs_domain_association[i] == pid)
        this_cpu_set.add_index(i);
    std::cout << pid << " Initialize scalar stuff 8\n";
    this_cpu_set.compress();
    std::cout << pid << " Initialize scalar stuff 9\n";
    TrilinosWrappers::MPI::Vector sol_scalar;
    std::cout << pid << " Initialize scalar stuff 10\n";
    sol_scalar.reinit(this_cpu_set, MPI_COMM_WORLD);
    std::cout << pid << " Initialize scalar stuff 11\n";
    MPI_Barrier(MPI_COMM_WORLD);

    //---------------------------------------------------------------------------
    // Initialize vector stuff:
    //---------------------------------------------------------------------------
    std::cout << pid << " Initialize vector stuff\n";
    dealii::FESystem<2, 3> fe_vector(fe_scalar, 3);
    dealii::DoFHandler     dh_vector(tria);
    dh_vector.distribute_dofs(fe_vector);
    //    DoFRenumbering::component_wise(dh_vector);
    std::vector<dealii::types::subdomain_id> dofs_domain_association_vector(dh_vector.n_dofs());
    dealii::DoFTools::get_subdomain_association(dh_vector, dofs_domain_association_vector);
    dealii::IndexSet this_cpu_set_vector(dh_vector.n_dofs());
    for (dealii::types::global_dof_index i = 0; i < dh_vector.n_dofs(); ++i)
      if (dofs_domain_association_vector[i] == pid)
        this_cpu_set_vector.add_index(i);
    this_cpu_set_vector.compress();
    TrilinosWrappers::MPI::Vector sol_vector;
    sol_vector.reinit(this_cpu_set_vector, MPI_COMM_WORLD);

    std::cout << pid << " Calculate support points\n";
    dealii::MappingQ<2, 3> mapping(1);

    std::vector<dealii::Point<3>> support_points(dh_scalar.n_dofs());
    dealii::DoFTools::map_dofs_to_support_points(mapping, dh_scalar, support_points);

    //---------------------------------------------------------------------------
    // Assign values to scalar and vector stuff:
    //---------------------------------------------------------------------------
    std::cout << pid << " Assign values to scalar and vector \n";
    std::vector<dealii::types::global_dof_index> local_dof_indices(fe_scalar.dofs_per_cell);
    std::vector<types::global_dof_index>         local_dof_indices_v(fe_vector.dofs_per_cell);
    auto                                         cell_v = dh_vector.begin_active();
    for (auto cell = dh_scalar.begin_active(); cell != dh_scalar.end(); ++cell)
    {
      if (cell->subdomain_id() == pid)
      {
        cell->get_dof_indices(local_dof_indices);
        cell_v->get_dof_indices(local_dof_indices_v);

        for (int i = 0; i < (int)local_dof_indices.size(); ++i)
        {
          auto   pos = support_points[local_dof_indices[i]];
          double pot, u, v, w;

          flowAroundSphere(0, 0, 0, pos[0], pos[1], pos[2], 1.0, 1.0, pot, u, v, w);
          sol_scalar[local_dof_indices[i]]           = pot;
          sol_vector[local_dof_indices_v[i * 3 + 0]] = u;
          sol_vector[local_dof_indices_v[i * 3 + 1]] = v;
          sol_vector[local_dof_indices_v[i * 3 + 2]] = w;
        }
      }
      ++cell_v;
    }

    const Vector<double> localized_sol_scalar(sol_scalar);
    const Vector<double> localized_sol_vector(sol_vector);
    if (pid == 0)
    {
      std::cout << pid << " Writing output ...\n";

      std::string                          filename_scalar = "sol_scalar_results.vtu";
      dealii::DataOut<2, DoFHandler<2, 3>> dataout_scalar;
      dataout_scalar.attach_dof_handler(dh_scalar);
      dataout_scalar.add_data_vector(localized_sol_scalar,
                                     "sol_scalar",
                                     dealii::DataOut<2, dealii::DoFHandler<2, 3>>::type_dof_data);
      dataout_scalar.build_patches(
        mapping, 1, dealii::DataOut<2, dealii::DoFHandler<2, 3>>::curved_inner_cells);
      std::ofstream file_scalar(filename_scalar.c_str());
      dataout_scalar.write_vtu(file_scalar);



      std::string filename_vector = "sol_vector_results.vtu";
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
          3, dealii::DataComponentInterpretation::component_is_part_of_vector);
      dealii::DataOut<2, DoFHandler<2, 3>> dataout_vector;
      dataout_vector.attach_dof_handler(dh_vector);
      dataout_vector.add_data_vector(localized_sol_vector,
                                     std::vector<std::string>(3, "sol_vector"),
                                     dealii::DataOut<2, dealii::DoFHandler<2, 3>>::type_dof_data,
                                     data_component_interpretation);
      dataout_vector.build_patches(
        mapping, 1, dealii::DataOut<2, dealii::DoFHandler<2, 3>>::curved_inner_cells);
      std::ofstream file_vector(filename_vector.c_str());
      dataout_vector.write_vtu(file_vector);
    }
    std::cout << pid << " Writing output, done\n";

    std::cout << pid << " Refining mesh adaptively ...\n";
    AdaptiveRefinement adaptiveRefinement(pcout, MPI_COMM_WORLD, errorEstimatorMax, aspectRatioMax);
    adaptiveRefinement.refine(
      pid, fe_scalar, fe_vector, dh_scalar, dh_vector, sol_scalar, sol_vector, tria);
    std::cout << pid << " Refining mesh adaptively, done\n";

    if (pid == 0)
    {
      std::cout << pid << " Writing new mesh ...\n";
      std::string   filename0 = ("meshResult.inp");
      std::ofstream logfile0(filename0.c_str());
      GridOut       grid_out0;
      grid_out0.write_ucd(tria, logfile0);
    }
  }
  std::cout << pid << " done\n";
}