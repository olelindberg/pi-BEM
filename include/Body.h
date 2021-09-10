#ifndef BODY_H
#define BODY_H

#include <vector>


class Body
{
public:
  Body()
  {}
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


private:
  std::vector<int> _materialIndices;
  std::vector<int> _waterlineIndices;
};

#endif // BODY_H