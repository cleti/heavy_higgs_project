

#include "../inc/FileBrowser.h"







int FileBrowser::ls_cdir()
{
  cdir_content.clear();
  int cnt = 0;
  for ( cdir_itr = boost::filesystem::directory_iterator(cdir_path);
	cdir_itr != end_itr;
	++cdir_itr )
    {
      if(is_regular_file(*cdir_itr) || is_directory(*cdir_itr) )
	{
	  if (apply_filter(cdir_itr->path())) continue;
	  std::cout << std::setw(5) << cnt << std::setw(20) << cdir_itr->path().filename() << std::endl;
	  cdir_content.push_back(*cdir_itr);
	  ++cnt;
	}
    }
  
  ext_filter.clear();
  return cnt;
}


int FileBrowser::select_entry(int entry_num)
{
  cdir_itr = boost::filesystem::directory_iterator(cdir_path);
  if (entry_num<0) return 0;
  for (int i=0;i<entry_num;++i)
    {
      ++cdir_itr;
      if (cdir_itr == end_itr) return 0;
    }
  return 1;
}

int FileBrowser::apply_filter(const boost::filesystem::path &path)
{
  std::string str = path.filename().string();
  if (str.front() == '.') return 1;
  if (str.back()  == '~') return 1;
  std::string ext = path.extension().string();
  if (!ext_filter.empty())
    {
      if (!(ext_filter == ext)) return 1;
    }
  return 0;
}


FileBrowser::FileBrowser(const std::string &base_path)
{
  if (!set_cdir(base_path))
    {
      std::cout << std::endl << " FileBrowser:  Path '" << base_path << "' not found." << std::endl;
    }
}

int FileBrowser::set_cdir(const boost::filesystem::path &path)
{
  if(boost::filesystem::exists(path))
    {
      cdir_path = path;
      return 1;
    }
  else
    {
      cdir_path = boost::filesystem::current_path();
      return 0;
    }
}

std::string FileBrowser::browse()
{
  int nentries = 0;
  int entry_num = -1;
  std::string user_input;
  std::stringstream user_input_stream;
  
  boost::filesystem::directory_iterator entry;
  while(1)
    {
      user_input_stream.clear();
      user_input_stream.flush();
      std::cout << std::endl << " Content of " << cdir_path << ":";
      if (!ext_filter.empty())
	{
	  std::cout << "  [filtering " << ext_filter << " files ]";
	}
      std::cout <<  std::endl;

      // list entries in the current directory
      nentries = ls_cdir();

      std::cout << std::endl << " Choose entry number [u:up, q:quit, .ext:set filter]: " << std::endl; 
      std::cout << " >> ";
      std::cin >> user_input;
      user_input_stream.str(user_input);

      char tc = ' ';
      // inspect the first character of user input (ignore spaces)
      while (tc==' ')
	{
	  user_input_stream.get(tc);
	}
      user_input_stream.unget();
      if (tc == '.')
	{
	  // set extension filter
	  user_input_stream >> ext_filter;
	  if (ext_filter == ".") ext_filter.clear();
	  continue;
	}
      else if (tc == 'u' || tc == 'U')
	{
	  // go to parent dir.
	  cdir_path /= "..";// cdir_path.parent_path();
	  continue;
	}
      else if (tc == 'q' || tc == 'Q')
	{
	  // quit
	  std::cout << std::endl << " FileBrowser:  Exit." << std::endl;
	  return "";
	}
      else if (tc >= '0' && tc <= '9')
	{
	  // extract entry number
	  user_input_stream >> entry_num;
	}
      else
	{
	  std::cout << std::endl << " Invalid input " << user_input << "." << std::endl;
	  continue;
	}
      
      if (entry_num>=0 && entry_num<nentries)
	{
	  // regular file selected
	  if (boost::filesystem::is_regular_file(cdir_content[entry_num]))
	    {
	      // the path to the regular file will be returned
	      break;
	    }
	  // directory selected
	  else if (boost::filesystem::is_directory(cdir_content[entry_num]))
	    {
	      if (!set_cdir(cdir_content[entry_num].path()))
		{
		  std::cout << std::endl << " Could not select new cdir: " << cdir_content[entry_num].path() << std::endl;
		}
	      continue;
	    }
	}
      else
	{
	  std::cout << std::endl << " Could not select entry " << entry_num << "." << std::endl;
	  entry_num = -1;
	}
    }
  return cdir_content[entry_num].path().string();
}

