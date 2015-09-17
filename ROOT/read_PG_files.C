

int read_PG_files()
{

  ifstream file;
  string filename;
  while (!file.is_open())
    {
      std::cout << std::endl << " Enter path+filename: " << std::endl;
      std::cout << " >> "; std::cin >> filename;
      file.open(filename);
    }
  std::cout << std::endl << " Reading data from file '" << filename << "' ..." << std::endl;


  const int colWidth1 = 16;
  const int colWidth2 = 20;


  std::string sbuff;
  double bin_lo;
  double lo_phi[3];
  double nlo_phi[3];
  double lo_qcd[3];
  double nlo_qcd[3];
  double tmp;
      while(1)
	{
	  if (file.good())
	    {
	      std::getline(file,sbuff);
	      // remove leading whitespace
	      while (sbuff.front() == ' ' || sbuff.front() == '\t')
		{
		  sbuff.erase(std::begin(sbuff));
		}
	      // inspect line
	      if (sbuff.size()>0) 
		{
		  if (sbuff.front() == '#')
		    { // print comment
		      std::cout << sbuff << std::endl;
		    }
		  else
		    { // data
		      std::stringstream iss(sbuff);
		      ///////////////////////////////
		      // lower bin boundary
		      iss >> bin_lo;
		      // LO PHI full, without int.
		      iss >> tmp;
		      iss >> tmp;
		      iss >> tmp;
		      // LO PHI full, with int.
		      iss >> lo_phi[0];
		      iss >> lo_phi[1];
		      iss >> lo_phi[2];
		      // LO PHI eff., without int.
		      iss >> tmp;
		      iss >> tmp;
		      iss >> tmp;
		      // LO PHI eff., with int.
		      iss >> tmp;
		      iss >> tmp;
		      iss >> tmp;
		      // NLO PHI, without int.
		      iss >> tmp;
		      iss >> tmp;
		      iss >> tmp;
		      // NLO PHI, with int.
		      iss >> nlo_phi[0];
		      iss >> nlo_phi[1];
		      iss >> nlo_phi[2];
		      // LO QCD
		      iss >> lo_qcd[0];
		      iss >> lo_qcd[1];
		      iss >> lo_qcd[2];
		      // NLO QCD
		      iss >> nlo_qcd[0];
		      iss >> nlo_qcd[1];
		      iss >> nlo_qcd[2];
		      ///////////////////////////////
		      // print relevant data
		      ost << std::setw(colWidth2) << std::setprecision(5)  << bin_lo;
		      ost << std::setw(colWidth2) << std::setprecision(10) << lo_qcd[0];
		      ost << std::setw(colWidth2) << std::setprecision(10) << lo_qcd[0]+lo_phi[0];
		      ost << std::setw(colWidth2) << std::setprecision(10) << nlo_qcd[0];
		      ost << std::setw(colWidth2) << std::setprecision(10) << nlo_qcd[0]+nlo_phi[0] << std::endl;
		    }
		}
	    }
	  else
	    {
	      if (file.eof())
		{
		  std::cout << "\n End-of-File reached on input operation. " << std::endl;
		  return 1;
		}
	      if (file.fail())
		{
		  std::cout << std::endl << " Error! ";
		  std::cout << "\n Logical error on i/o operation. " << std::endl;
		  return 0;
		}
	      if (file.bad())
		{
		  std::cout << std::endl << " Error! ";
		  std::cout << "\n Read/writing error on i/o operation. " << std::endl;
		  return 0;
		};
	    }
	}









  return 1;
}
