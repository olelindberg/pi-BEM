#ifndef BODY_SETTINGS_READER_JSON_H
#define BODY_SETTINGS_READER_JSON_H
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <exception>
#include <iostream>
#include <vector>

#include "Body.h"



class BodySettingsReaderJSON
{
public:
  static bool
  read(const std::string &filename, Body &body)
  {
    namespace pt = boost::property_tree;
    pt::ptree root;
    try
    {
      pt::read_json(filename, root);

      auto name = root.get<std::string>("name");
      body.setName(name);

      auto draft = root.get<double>("draft");
      body.setDraft(draft);

      std::vector<int> materialIndices;
      for (pt::ptree::value_type &id : root.get_child("materialIndices"))
        materialIndices.push_back(id.second.get_value<int>());
      body.setMaterialIndices(materialIndices);

      std::vector<int> waterlineIndices;
      for (pt::ptree::value_type &id : root.get_child("waterlineIndices"))
        waterlineIndices.push_back(id.second.get_value<int>());
      body.setWaterlineIndices(waterlineIndices);

      std::vector<double> centerOfGravity;
      for (pt::ptree::value_type &id : root.get_child("centerOfGravity"))
        centerOfGravity.push_back(id.second.get_value<double>());
      Tensor<1, 3> tmp;
      for (int i = 0; i < 3; ++i)
        tmp[i] = centerOfGravity[i];
      body.setCenterOfGravity(tmp);



      return true;
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      return false;
    }
  }

private:
};

#endif // BODY_SETTINGS_READER_JSON_H
