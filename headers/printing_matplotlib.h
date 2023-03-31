// #include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
// #include <cppcolormap.h> // https://github.com/tdegeus/cppcolormap


// function to simplify imshow plotting
void Imshow_matrix(vector<float> &M, int rows, int col, string &Title, string &File_name, bool show){

  const int colors = 1;
  const float* zptr = &(M[0]); // pointer to put into imshow. It works weird on C++
  // auto cmap = cppcolormap::jet();
  plt::title(Title);
  PyObject* mat;
  plt::imshow( zptr,rows,col,colors, {{"cmap", "jet"}}, &mat);
  plt::colorbar(mat);
  if(show){
    plt::show();
  }
  //plt::show();
  plt::save(File_name);
  plt::close();
  Py_DECREF(mat);

}
