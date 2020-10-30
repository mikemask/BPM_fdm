class Bond
{
public:
  Bond() :
    sigmac_(0.),
    tauc_(0.),
    rad_(0.),
    kn_(0.),
    ks_(0.),
    fn_(0.),
    fs_(0.),
    partid_{0, 0},
    inertia_{0., 0.},
    r_old_{0., 0., 0.},
    r_old2_{0., 0., 0.},
    vc_old_{0., 0., 0.},
    force_n_{0., 0., 0.},
    force_s_{0., 0., 0.},
    torque_n_{0., 0., 0.},
    torque_s_{0., 0., 0.},
    torque_t_{0., 0., 0.},
    torque_tr_{0., 0., 0.}
  {}

  Bond(double sigma, double tau, double radius, double kn, double ks, double fn, double fs) :
    sigmac_(sigma),
    tauc_(tau),
    rad_(radius),
    kn_(kn),
    ks_(ks),
    fn_(fn),
    fs_(fs)
  {}

  Bond(int i, int j, double I, double J) :
    partid_{i, j},
    inertia_{I, J}
  {}

  Bond(double r_oldx, double r_oldy, double r_oldz, double r_old2x, double r_old2y, double r_old2z, double vc_oldx, double vc_oldy, double vc_oldz, double fnx, double fny, double fnz, double fsx, double fsy, double fsz, double tnx, double tny, double tnz, double tsx, double tsy, double tsz, double ttx, double tty, double ttz, double ttrx, double ttry, double ttrz) :
    r_old_{r_oldx, r_oldy, r_oldz},
    r_old2_{r_old2x, r_old2y, r_old2z},
    vc_old_{vc_oldx, vc_oldy, vc_oldz},
    force_n_{fnx, fny, fnz},
    force_s_{fsx, fsy, fsz},
    torque_n_{tnx, tny, tnz},
    torque_s_{tsx, tsy, tsz},
    torque_t_{ttx, tty, ttz},
    torque_tr_{ttrx, ttry, ttrz}
  {}


  ~Bond() {}

  double getSigma()
  { return sigmac_; }

  double getTau()
  { return tauc_; }

  double getRad()
  { return rad_; }

  double getStiffNorm()
  { return kn_; }

  double getStiffShear()
  { return ks_; }

  double getFn()
  {return fn_; }

  double getFs()
  {return fs_; }

  int* getIds()
  { return partid_; }

  double* getInertia()
  { return inertia_; }

  double* getRold()
  { return r_old_; }

  double* getRold2()
  { return r_old2_; }

  double* getVcold()
  { return vc_old_; }

  double* getForce_n()
  { return force_n_; }

  double* getForce_s()
  { return force_s_; }

  double* getTorque_n()
  { return torque_n_; }

  double* getTorque_s()
  { return torque_s_; }

  double* getTorque_t()
  { return torque_t_; }

  double* getTorque_tr()
  { return torque_tr_; }


  void setSigma(double sigma)
  { sigmac_ = sigma; }

  void setTau(double tau)
  { tauc_ = tau; }

  void setRad(double radius)
  { rad_ = radius; }

  void setStiffNorm(double kn)
  { kn_ = kn; }

  void setStiffShear(double ks)
  { ks_ = ks; }

  void setFn(double fn)
  { fn_ = fn; }

  void setFs(double fs)
  { fs_ = fs; }

  void setIds(int* ids)
  {
    partid_[0] = ids[0];
    partid_[1] = ids[1];
  }

  void setInertia(double* in)
  {
    inertia_[0] = in[0];
    inertia_[1] = in[1];
  }

  void setRold(double* r_old)
  {
    r_old_[0] = r_old[0];
    r_old_[1] = r_old[1];
    r_old_[2] = r_old[2];
  }

  void setRold2(double* r_old2)
  {
    r_old2_[0] = r_old2[0];
    r_old2_[1] = r_old2[1];
    r_old2_[2] = r_old2[2];
  }

  void setVcold(double* vc_old)
  {
    vc_old_[0] = vc_old[0];
    vc_old_[1] = vc_old[1];
    vc_old_[2] = vc_old[2];
  }

  void setForce_n(double* fn)
  {
    force_n_[0] = fn[0];
    force_n_[1] = fn[1];
    force_n_[2] = fn[2];
  }

  void setForce_s(double* fs)
  {
    force_s_[0] = fs[0];
    force_s_[1] = fs[1];
    force_s_[2] = fs[2];
  }

  void setTorque_n(double* tn)
  {
    torque_n_[0] = tn[0];
    torque_n_[1] = tn[1];
    torque_n_[2] = tn[2];
  }

  void setTorque_s(double* ts)
  {
    torque_s_[0] = ts[0];
    torque_s_[1] = ts[1];
    torque_s_[2] = ts[2];
  }

  void setTorque_t(double* tt)
  {
    torque_t_[0] = tt[0];
    torque_t_[1] = tt[1];
    torque_t_[2] = tt[2];
  }

  void setTorque_tr(double* ttr)
  {
    torque_tr_[0] = ttr[0];
    torque_tr_[1] = ttr[1];
    torque_tr_[2] = ttr[2];
  }


private:
  double sigmac_, tauc_, rad_, kn_, ks_, fn_, fs_, inertia_[2];
  int partid_[2];
  double r_old_[3], r_old2_[3], vc_old_[3], force_n_[3], force_s_[3], torque_n_[3], torque_s_[3], torque_t_[3], torque_tr_[3];

};




