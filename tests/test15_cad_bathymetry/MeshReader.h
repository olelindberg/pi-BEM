#ifndef MESH_READER_H
#define MESH_READER_H


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <fstream>

class MeshReader
{
public:
  static std::shared_ptr<dealii::Triangulation<2, 3>>
  read(const std::filesystem::path &root, const std::vector<GeometryInput> &inputs)
  {
    //---------------------------------------------------------------------------
    // Read meshes:
    //---------------------------------------------------------------------------
    std::vector<std::shared_ptr<dealii::Triangulation<2, 3>>> meshes;
    for (auto &input : inputs)
    {
      auto filename = std::filesystem::path(root).append(input.get_filename());
      std::cout << "Reading inp mesh file: " << filename << std::endl;
      std::ifstream file;
      file.open(filename);
      auto mesh = std::make_shared<dealii::Triangulation<2, 3>>();
      if (file.is_open())
      {
        dealii::GridIn<2, 3> grid_in;
        grid_in.attach_triangulation(*mesh);
        grid_in.read_ucd(file, true);
        file.close();
      }
      else
      {
        std::cout << filename << std::endl;
        std::cout << "Error" << std::endl;
        // Error handling
      }

      dealii::Tensor<1, 3> pre_translate;
      pre_translate[0] = input.get_pre_translate_x();
      pre_translate[1] = input.get_pre_translate_y();
      pre_translate[2] = input.get_pre_translate_z();
      dealii::GridTools::shift(pre_translate, *mesh);

      dealii::GridTools::scale(input.get_scale(), *mesh);

      dealii::GridTools::rotate(input.get_rotate_x(), 0, *mesh);
      dealii::GridTools::rotate(input.get_rotate_y(), 1, *mesh);
      dealii::GridTools::rotate(input.get_rotate_z(), 2, *mesh);

      dealii::Tensor<1, 3> post_translate;
      post_translate[0] = input.get_post_translate_x();
      post_translate[1] = input.get_post_translate_y();
      post_translate[2] = input.get_post_translate_z();
      dealii::GridTools::shift(post_translate, *mesh);

      meshes.push_back(mesh);
    }

    auto mesh = std::make_shared<dealii::Triangulation<2, 3>>();
    for (auto &tmp : meshes)
      dealii::GridGenerator::merge_triangulations(
        *tmp, *mesh, *mesh, std::numeric_limits<double>::epsilon(), true);

    return mesh;
  }
};

#endif // MESH_READER_H