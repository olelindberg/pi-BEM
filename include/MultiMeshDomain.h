#ifndef MULTI_MESH_DOMAIN_H
#define MULTI_MESH_DOMAIN_H

#include <IPhysicalDomain.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/grid/manifold_lib.h>

template <int dim>
class MultiMeshDomain : public IPhysicalDomain<dim>
{
public:
  MultiMeshDomain(MPI_Comm comm = MPI_COMM_WORLD);
  ~MultiMeshDomain();

  virtual void read_domain(std::string input_path = "") override;

  virtual void refine_and_resize(std::string input_path = "") override;

  virtual void update_triangulation() override;

  virtual const dealii::Triangulation<dim - 1, dim> &getTria() const override
  {
    return *_mesh;
  }

  virtual dealii::Triangulation<dim - 1, dim> &getTria() override
  {
    return *_mesh;
  }

  virtual std::vector<unsigned int> get_dirichlet_boundary_ids() override
  {
    return dirichlet_boundary_ids;
  };

  virtual std::vector<unsigned int> get_neumann_boundary_ids() override
  {
    return neumann_boundary_ids;
  };


private:
  MPI_Comm     mpi_communicator;
  unsigned int n_mpi_processes;
  unsigned int this_mpi_process;

  std::shared_ptr<dealii::Triangulation<2, 3>> _mesh;
  std::vector<unsigned int>                    dirichlet_boundary_ids;
  std::vector<unsigned int>                    neumann_boundary_ids;

  std::shared_ptr<dealii::TransfiniteInterpolationManifold<2, 3>> inner_manifold;

  dealii::ConditionalOStream pcout;
};

#endif // MULTI_MESH_DOMAIN_H1
