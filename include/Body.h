#ifndef BODY_H
#define BODY_H

#include <deal.II/base/tensor.h>

#include <vector>


class Body
{
public:
  Body()
  {
    _centerOfGravity = 0.0;
  }
  virtual ~Body()
  {}

  void
  setMaterialIndices(const std::vector<int> &materialIndices)
  {
    _materialIndices = materialIndices;
  }

  void
  setWaterlineIndices(const std::vector<int> &indices)
  {
    _waterlineIndices = indices;
  }

  const std::vector<int> &
  getMaterialIndices() const
  {
    return _materialIndices;
  }

  bool
  hasMaterial(int index) const
  {
    if (std::find(_materialIndices.begin(), _materialIndices.end(), index) !=
        _materialIndices.end())
      return true;
    else
      return false;
  }

  bool
  isWaterline(int index) const
  {
    if (std::find(_waterlineIndices.begin(), _waterlineIndices.end(), index) !=
        _waterlineIndices.end())
      return true;
    else
      return false;
  }

  void
  setDraft(double draft)
  {
    _draft = draft;
  }

  double
  getDraft() const
  {
    return _draft;
  }

  void
  setName(const std::string &name)
  {
    _name = name;
  }

  const std::string &
  getName() const
  {
    return _name;
  }

  void
  setCenterOfGravity(const dealii::Tensor<1, 3> &centerOfGravity)
  {
    _centerOfGravity = centerOfGravity;
  }

  const dealii::Tensor<1, 3> &
  getCenterOfGravity() const
  {
    return _centerOfGravity;
  }


  void
  print()
  {
    std::cout << "Printing body parameters ..." << std::endl;
    std::cout << "name             : " << _name << std::endl;
    std::cout << "draft            : Tm = " << _draft << std::endl;

    std::cout << "material indices : ";
    for (const auto &id : _materialIndices)
      std::cout << id << ", ";
    std::cout << std::endl;

    std::cout << "waterline indices: ";
    for (const auto &id : _waterlineIndices)
      std::cout << id << ", ";
    std::cout << std::endl;

    std::cout << "center of gravity: COG = " << _centerOfGravity << std::endl;
  }

private:
  std::string          _name;
  double               _draft = 0.0;
  std::vector<int>     _materialIndices;
  std::vector<int>     _waterlineIndices;
  dealii::Tensor<1, 3> _centerOfGravity;
};

#endif // BODY_H