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


  bool
  saveScalarFields(std::string                                   filename,
                   const dealii::DoFHandler<2, 3> &              dh,
                   const std::shared_ptr<dealii::Mapping<2, 3>> &mapping,
                   unsigned int                                  mapping_degree)
  {
    std::cout << "hello1\n";
    dealii::DataOut<2, dealii::DoFHandler<2, 3>> dataout_scalar;
    std::cout << "hello2\n";
    dataout_scalar.attach_dof_handler(dh);
    std::cout << "hello3\n";

    auto nameIt = _scalarFieldNames.begin();
    for (auto &field : _scalarFields)
    {
      std::cout << "hello4\n";
      dataout_scalar.add_data_vector(field,
                                     *nameIt,
                                     dealii::DataOut<2, dealii::DoFHandler<2, 3>>::type_cell_data);
      std::cout << "hello5\n";
      ++nameIt;
      std::cout << "hello6\n";
    }

    std::cout << "hello7\n";

    dataout_scalar.build_patches(*mapping,
                                 mapping_degree,
                                 dealii::DataOut<2, dealii::DoFHandler<2, 3>>::curved_inner_cells);
    std::cout << "hello8\n";

    std::cout << "hello9\n";
    std::ofstream file(filename);
    std::cout << "hello10\n";
    if (file.is_open())
      dataout_scalar.write_vtu(file);
    else
    {
      std::cout << "Error: Failed saving " << filename << std::endl;
      return false;
    }
    std::cout << "hello11\n";
    return true;
  }


  void
  addScalarField(std::string name, const dealii::Vector<double> &scalarField)
  {
    _scalarFieldNames.push_back(name);
    _scalarFields.push_back(scalarField);
  }

private:
  std::vector<std::string>            _scalarFieldNames;
  std::vector<dealii::Vector<double>> _scalarFields;
};

#endif // WRITER_H
