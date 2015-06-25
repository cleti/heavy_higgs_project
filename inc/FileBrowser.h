

#ifndef FILEBROWSER_H
#define FILEBROWSER_H 

#include <iostream>
#include <iomanip>
#include <sstream>
#include "boost/filesystem.hpp" 

class FileBrowser
{
private:
  boost::filesystem::path               cdir_path;
  boost::filesystem::directory_iterator cdir_itr;
  boost::filesystem::directory_iterator end_itr;
  std::vector<boost::filesystem::directory_entry> cdir_content;
  std::string ext_filter;
  
  int ls_cdir();
  std::string get_entry_path(int entry_num);
  int select_entry(int entry_num);
  int apply_filter(const boost::filesystem::path &path);
  
public:
  FileBrowser(const std::string &base_path = ".");
  ~FileBrowser() {}

  std::string browse();
  int set_cdir(const boost::filesystem::path &path);
  
};



#endif
