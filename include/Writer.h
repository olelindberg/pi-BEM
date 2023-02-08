#ifndef WRITER_H
#define WRITER_H

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

class Writer
{
public:
  Writer(){};
  virtual ~Writer()
  {}


  bool saveScalarFields(std::string                     filename,
                        const dealii::DoFHandler<2, 3> &dh,
                        const dealii::Mapping<2, 3> &   mapping,
                        unsigned int                    mapping_degree)
  {
    dealii::DataOut<2, 3> dataout;
    dataout.attach_dof_handler(dh);

    auto nameIt = _scalarFieldNames.begin();
    for (auto &field : _scalarFields)
    {
      dataout.add_data_vector(field, *nameIt, dealii::DataOut<2, 3>::type_cell_data);
      ++nameIt;
    }

    dataout.build_patches(mapping, mapping_degree, dealii::DataOut<2, 3>::curved_inner_cells);
    dataout.build_patches();

    std::ofstream file(filename);
    if (file.is_open())
      dataout.write_vtu(file);
    else
    {
      std::cout << "Error: Failed saving " << filename << std::endl;
      return false;
    }
    return true;
  }


  void addScalarField(std::string name, const dealii::Vector<double> &scalarField)
  {
    _scalarFieldNames.push_back(name);
    _scalarFields.push_back(scalarField);
  }

  bool save(std::string filename, const dealii::Triangulation<2, 3> &trimesh)
  {
    dealii::DataOut<2, 3> dataout;
    dataout.attach_triangulation(trimesh);
    dataout.build_patches();
    std::ofstream file(filename);
    if (file.is_open())
      dataout.write_vtu(file);
    else
    {
      std::cout << "Error: Failed saving " << filename << std::endl;
      return false;
    }
    return true;
  }

private:
  std::vector<std::string>            _scalarFieldNames;
  std::vector<dealii::Vector<double>> _scalarFields;
};

#endif // WRITER_H
