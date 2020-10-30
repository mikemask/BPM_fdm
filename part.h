class Part
{
public:
    Part() :
        r_{0., 0., 0.},
        r_in_{0., 0., 0.},
        r_old_{0., 0., 0.},
        U_{0., 0., 0.},
        U_old_{0., 0., 0.},
        om_{0., 0., 0.},
        f_{0., 0., 0.},
        f_old_{0., 0., 0.},
        t_{0., 0., 0.},
        t_old_{0., 0., 0.},
        m_(0.),
        d_(0.),
        id_(-1)
    {}

    Part(double x, double y, double z, double x_in, double y_in, double z_in, double x_old, double y_old, double z_old, double u, double v, double w, double omx, double omy, double omz, double fx, double fy, double fz, double tx, double ty, double tz) :
        r_{x, y, z},
        r_in_{x_in, y_in, z_in},
        r_old_{x_old, y_old, z_old},
        U_{u, v, w},
        om_{omx, omy, omz},
        f_{fx, fy, fz},
        t_{tx, ty, tz}
    {}

    Part(double u_old, double v_old, double w_old, double fx_old, double fy_old, double fz_old, double tx_old, double ty_old, double tz_old) :
        U_old_{u_old, v_old, w_old},
        f_old_{fx_old, fy_old, fz_old},
        t_old_{tx_old, ty_old, tz_old}
    {}

    Part(double m, double d, int id) :
        m_(m),
        d_(d),
        id_(id)
    {}

    ~Part() {}

    double* getPos()
    { return r_; }

    double* getPosIn()
    { return r_in_; }

    double* getPosOld()
    { return r_old_; }

    double* getVel()
    { return U_; }

    double* getVelOld()
    { return U_old_; }

    double* getRot()
    { return om_; }

    double* getForce()
    { return f_; }

    double* getForceOld()
    { return f_old_; }

    double* getTorque()
    { return t_; }

    double* getTorqueOld()
    { return t_old_; }

    double getMass()
    { return m_; }

    double getDiameter()
    { return d_; }

    double getId()
    { return id_; }


    void setPos(double* r)
    {
        r_[0] = r[0];
        r_[1] = r[1];
        r_[2] = r[2];
    }

    void setPosIn(double* r_in)
    {
        r_in_[0] = r_in[0];
        r_in_[1] = r_in[1];
        r_in_[2] = r_in[2];
    }

    void setPosOld(double* r_old)
    {
        r_old_[0] = r_old[0];
        r_old_[1] = r_old[1];
        r_old_[2] = r_old[2];
    }

    void setVel(double* U)
    {
        U_[0] = U[0];
        U_[1] = U[1];
        U_[2] = U[2];
    }

    void setVelOld(double* U_old)
    {
        U_old_[0] = U_old[0];
        U_old_[1] = U_old[1];
        U_old_[2] = U_old[2];
    }

    void setRot(double* om)
    {
        om_[0] = om[0];
        om_[1] = om[1];
        om_[2] = om[2];
    }

    void setForce(double* f)
    {
        f_[0] = f[0];
        f_[1] = f[1];
        f_[2] = f[2];
    }

    void setForceOld(double* f_old)
    {
        f_old_[0] = f_old[0];
        f_old_[1] = f_old[1];
        f_old_[2] = f_old[2];
    }

    void setTorque(double* t)
    {
        t_[0] = t[0];
        t_[1] = t[1];
        t_[2] = t[2];
    }

    void setTorqueOld(double* t_old)
    {
        t_old_[0] = t_old[0];
        t_old_[1] = t_old[1];
        t_old_[2] = t_old[2];
    }

    void setMass(double m)
    { m_ = m; }

    void setDiameter(double d)
    { d_ = d; }

    void setId(int id)
    { id_ = id; }


private:
    double r_[3], r_in_[3], r_old_[3], U_[3], U_old_[3], om_[3], f_[3], f_old_[3], t_[3], t_old_[3], m_, d_;
    int id_;
};
