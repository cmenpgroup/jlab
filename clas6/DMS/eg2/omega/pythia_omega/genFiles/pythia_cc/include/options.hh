#ifndef options_hh
#define options_hh 1

#include "cpp_headers.hh"

bool parse_file(const char*);
string get_line_data(const char*, const char*);
string get_next_line_data(const char* tag, const char* key);
//void print_state (const std::ios& stream);

template <class T>
bool get_opt(const char* tag, const char* key, vector< T >& values)
{
  string data = get_line_data(tag,key);
  //cout << "data (1) is " << data << ", length is " << data.length()  << endl;
  if(data.length() == 0) return false;
  stringstream in(data.c_str());  
  while(1){
    T tmp;
    in>>tmp;
    if(!in) break;
    values.push_back(tmp);
  }
  return true;
}


template <class T>
bool get_opt(const char* tag, const char* key, T& value)
{
  string data = get_line_data(tag,key);
  //cout << "data (2) is " << data << ", length is " << data.length()  << endl;
  stringstream in(data.c_str());  
  //print_state(in);
  //if(data.length() == 0) return false;
  in>>value;
  //cout << "value (2) is "<<value<<endl;
  //print_state(in);

  if(in.fail())
   return false;
  else
   return true;

/*
  if(in.good())
  {
    cout << "in is good"<<endl;
   return true;
  }
  else
  {
    cout << "in is not good"<<endl;
   return false;
  }
*/
}

template <class T>
bool get_all_opt(const char* tag, const char* key, T& value)
{
  string data = get_next_line_data(tag, key);
  if(data.empty()) return false;
  
  stringstream in(data.c_str());  
  in>>value;
  if(in) {
    return true;
  }
  else {
    return false;  
  }
}

/*
void print_state (const std::ios& stream) {
  std::cout << " good()=" << stream.good();
  std::cout << " eof()=" << stream.eof();
  std::cout << " fail()=" << stream.fail();
  std::cout << " bad()=" << stream.bad()<<endl;
}
*/

#endif
