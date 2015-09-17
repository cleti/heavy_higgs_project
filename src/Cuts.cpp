

#include "../inc/Cuts.h"



Cut::Cut(OBSFnc obs, const std::initializer_list<double>& il): d_obs(obs), d_intervals(il)
{
  if(!check_intervals())
    {
       WARNING("object not properly constructed, use set_intervals() to reset");
    }
}

bool Cut::check_intervals()
{
  double temp = -DBL_MAX;
  for (auto it=d_intervals.begin(); it!=d_intervals.end();it++)
    {
      if (*it<=temp)
	{
	  WARNING("need intervals in ascending order");
	  return false;
	}
	temp = *it;
    }
  return true;
}

bool Cut::set_intervals(const std::initializer_list<double>& il)
{
  d_intervals = il;
  return check_intervals();
}

bool Cut::operator()(const PS_2* ps_lab, const PS_2* ps_tt) const
{
  double obs_val = d_obs(ps_lab,ps_tt);
  
// #ifdef DEBUG  
//   PRINT(obs_val);
// #endif
  
  int len = d_intervals.size();
  for (int i=0; i<len;i+=2)
    {
// #ifdef DEBUG
//       std::cout << "interval : [" << d_intervals[i] << ", ";
//       if (len-i == 1)
// 	{
// 	  std::cout <<  "+infty ]\n";
// 	}
//       else
// 	{
// 	  std::cout << d_intervals[i+1] << " ]\n";
// 	}
// #endif
      // obs_val lies above lower bound?
      if (obs_val > d_intervals[i])
	{
	  // if this is the last value
	  // obs_val lies within the unbounded interval [d_intervals[i],+infty]
	  if (len-i == 1) return true;
	  else
	    {
	      // else check also the upper limit
	      if (obs_val < d_intervals[i+1]) return true;
	    }
	}
    }
  return false;
}








bool EvalCuts(const CutVec* cuts,const PS_2* ps_lab, const PS_2* ps_tt)
{
  if (cuts)
    {
      for (auto const& cut: *cuts)
	{
	  if (!(*cut)(ps_lab,ps_tt)) return false;
	}
    }
  return true;
}
