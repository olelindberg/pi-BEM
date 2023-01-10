#ifndef ERROR_MESSAGE_H
#define ERROR_MESSAGE_H

#include <string>

class ErrorMessage
{
public:
  static void
  print(std::string file, int line, std::string msg);
  static std::string
  message(std::string file, int line, std::string msg);
};

#endif // ERROR_MESSAGE_H