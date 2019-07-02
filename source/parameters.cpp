#include "parameters.h"

Parameters::Parameters(std::string file_name, std::string ctcf,
                       std::string prob) {
  this->read_parm(file_name);
  this->read_ctcf(ctcf);
  this->read_coupling_prob(prob);
}

Parameters::~Parameters() {}

bool Parameters::read_parm(std::string file_name) {
  // read file
  std::ifstream parm_file;
  parm_file.open(file_name, std::fstream::in);

  // return error if read file fail
  if (parm_file.fail()) {
    throw "ERROR: Impossible to open parameters file " + file_name;
    return 1;
  }

  std::string line;
  while (std::getline(parm_file, line)) {
    std::istringstream is_line(line);
    std::string name;
    if (std::getline(is_line, name, '=')) {
      std::string par;
      if (std::getline(is_line, par)) {
        if (name == "nstep")
          set_parm(m_nstep, std::stoi(par));
        else if (name == "print")
          set_parm(m_print, std::stoi(par));
        else if (name == "timestep")
          set_parm(m_timestep, std::stof(par));
        else if (name == "temp")
          set_parm(m_temp, std::stof(par));
        else if (name == "init")
          set_parm(m_init, par);
        else if (name == "nsphere")
          set_parm(m_psphere, std::stoi(par));
        else if (name == "distance")
          set_parm(m_distance, std::stof(par));
        else if (name == "spring_k")
          set_parm(m_spring, std::stof(par));
        else if (name == "rinit")
          set_parm(m_rinit, std::stof(par));
        else if (name == "epsilon")
          set_parm(m_epsilon, std::stof(par));
        else if (name == "rmin")
          set_parm(m_rmin, std::stof(par));
        else if (name == "rcut")
          set_parm(m_rcut, std::stof(par));
        else if (name == "gamma")
          set_parm(m_gamma, std::stof(par));
        else if (name == "k_on")
          set_parm(m_k_on, std::stof(par));
        else if (name == "k_off")
          set_parm(m_k_off, std::stof(par));
        else if (name == "rate_fwl")
          set_parm(m_rate_fwl, std::stof(par));
        else if (name == "rate_fwr")
          set_parm(m_rate_fwr, std::stof(par));
        else if (name == "rate_bwl")
          set_parm(m_rate_bwl, std::stof(par));
        else if (name == "rate_bwr")
          set_parm(m_rate_bwr, std::stof(par));
        else if (name == "n_max_extr")
          set_parm(m_n_max_extr, std::stoi(par));
        else if (name == "permeability_ctcf")
          set_parm(m_perm_ctcf, std::stof(par));
      }
    }
  }
  parm_file.close();
  return 0;
}

bool Parameters::read_ctcf(std::string file_name) {
  int x;

  // read file
  std::ifstream ctcf_file;
  ctcf_file.open(file_name, std::fstream::in);

  // return error if read file fail
  if (ctcf_file.fail()) {
    throw "ERROR: Impossible to open parameters file " + file_name;
    return 1;
  }

  std::string line;
  while (std::getline(ctcf_file, line)) {
    std::istringstream iss(line);
    iss >> x;
    m_ctcf.push_back(x);
  }
  ctcf_file.close();
  return 0;
}

bool Parameters::read_coupling_prob(std::string file_name) {
  double p;

  // read file
  std::ifstream coupling_prob_file;
  coupling_prob_file.open(file_name, std::fstream::in);

  // return error if read file fail
  if (coupling_prob_file.fail()) {
    throw "ERROR: Impossible to open parameters file " + file_name;
    return 1;
  }

  std::string line;
  while (std::getline(coupling_prob_file, line)) {
    std::istringstream iss(line);
    iss >> p;
    m_coupling_prob.push_back(p);
  }
  coupling_prob_file.close();
  return 0;
}
