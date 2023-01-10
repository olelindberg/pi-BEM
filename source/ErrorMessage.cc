#include "ErrorMessage.h"

#include <iostream>
#include <sstream>

void
ErrorMessage::print(std::string file, int line, std::string msg)
{
  std::cout << message(file, line, msg);
}

std::string
ErrorMessage::message(std::string file, int line, std::string msg)
{
  std::stringstream sstrm;
  sstrm << "Error in " << file << "(" << line << "):\n";
  sstrm << msg << "\n";
  return sstrm.str();
}