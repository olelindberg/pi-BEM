#ifndef JSON_BODYREADER_H
#define JSON_BODYREADER_H
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <exception>
#include <iostream>
#include <vector>

#include "Body.h"



class JSON_BodyReader
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

      std::vector<int> materialIndices;
      for (pt::ptree::value_type &id : root.get_child("materialIndices"))
        materialIndices.push_back(id.second.get_value<int>());
      body.setMaterialIndices(materialIndices);

      std::vector<int> waterlineIndices;
      for (pt::ptree::value_type &id : root.get_child("waterlineIndices"))
        waterlineIndices.push_back(id.second.get_value<int>());
      body.setWaterlineIndices(waterlineIndices);

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

#endif // JSON_BODYREADER_H
