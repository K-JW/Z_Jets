/*
 * FileName: particle.h
 * Version: 1.0
 * License: MIT License
 *
 * Author: KANG Jin-Wen
 * E-Mail: kangjinwen@vip.qq.com
 * Date: 2019-10-15 16:51:18
 *
 * LastEditors: KANG Jin-Wen
 * LastEditTime: 2019-11-03 10:21:43
 * Description: Define a class that include some information of particle
 *              Refs PseduoJet.hh
 */

#ifndef __IHEPTOOLS_PARTICLE_H__
#define __IHEPTOOLS_PARTICLE_H__

#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <valarray>

using std::cout;

namespace iHepTools {

static const double PI = acos(-1);
static const double TwoPI = 2.0 * PI;
static constexpr double EPSILON = 1.0E-20;
static constexpr int DefaultPDGID = 0;
static constexpr double NaNValue = std::numeric_limits<double>::quiet_NaN();
static constexpr double INFVALUE = std::numeric_limits<double>::infinity();

class Particle {
public:
    // ------------------- Constructors -------------------------------
    Particle()
        : _E(0.),
          _px(0.),
          _py(0.),
          _pz(0.),
          _pdg_id(DefaultPDGID),
          _mass(0.),
          _phi(NaNValue),
          _rap(NaNValue),
          _pt2(0.),
          _user_index(-1) {}
    Particle(const int &pdg_id, const double &E, const double &px,
             const double &py, const double &pz, const double &mass)
        : _pdg_id(pdg_id), _E(E), _px(px), _py(py), _pz(pz), _mass(mass) {
        _finish_init();
        _user_index = -1;
    }
    Particle(const double &E, const double &px, const double &py,
             const double &pz, const double &mass)
        : _pdg_id(DefaultPDGID), _E(E), _px(px), _py(py), _pz(pz), _mass(mass) {
        _finish_init();
        _user_index = -1;
    }
    Particle(const int &pdg_id, const double &E, const double &px,
             const double &py, const double &pz)
        : _pdg_id(pdg_id), _E(E), _px(px), _py(py), _pz(pz) {
        _set_mass();
        _finish_init();
        _user_index = -1;
    }
    Particle(const double &E, const double &px, const double &py,
             const double &pz)
        : _pdg_id(DefaultPDGID), _E(E), _px(px), _py(py), _pz(pz) {
        _set_mass();
        _finish_init();
        _user_index = -1;
    }
    // ------------------- END OF CONSTRUCTORS -------------------------

    // ------------------- PUBLIC: Member Variables --------------------
    inline int pdg_id() const { return _pdg_id; }
    inline void reset_pdg_id(int pdg_id) { _pdg_id = pdg_id; }
    inline double E() const { return _E; }
    inline double e() const { return _E; }
    inline double px() const { return _px; }
    inline double py() const { return _py; }
    inline double pz() const { return _pz; }
    inline double m() const { return _mass; }
    inline double mass() const { return _mass; }
    inline double m2() const { return _mass * _mass; }
    // return the squared mass from calculate
    inline double m2_calc() const { return _get_mass_square(); }
    // return the mass from calculte
    inline double m_calc() const {
        return _get_mass_square() > 0 ? sqrt(_get_mass_square())
                                      : -sqrt(-_get_mass_square());
    }
    // return the squared transverse mass = pt^2 + m^2
    inline double mt2() const { return (_E + _pz) * (_E - _pz); }
    // return the transverse mass = sqrt(pt^2 + m^2)
    inline double mt() const { return sqrt(std::fabs(mt2())); }
    // return the squared transverse momentum
    inline double pt2() const { return _pt2; }
    // return the scalar transverse momentum
    inline double pt() const { return sqrt(_pt2); }

    // returns phi (in the range 0 to 2*PI)
    inline double phi() const {
        _ensure_valid_rap_phi();
        return _phi;
    }
    // returns phi in the range -PI to PI
    inline double phi_std() const {
        return this->phi() > PI ? this->phi() - TwoPI : this->phi();
    }
    // return rapidity
    inline double rap() const {
        _ensure_valid_rap_phi();
        return _rap;
    }
    inline double rapidity() const { return this->rap(); }
    // return the psedo-rapidity
    double pseudorapidity() const;
    inline double eta() const { return pseudorapidity(); }
    // return the squared 3-vector modulus = px^2 + py^2 + pz^2
    inline double mod_p_2() const { return _pt2 + _pz * _pz; }
    // return the 3-vector modulus = sqrt(px^2 + py^2 + pz^2)
    inline double mod_p() const { return sqrt(_pt2 + _pz * _pz); }
    // return the transverse energy
    inline double Et() const {
        return (_pt2 <= EPSILON) ? 0.0 : _E / sqrt(1.0 + _pz * _pz / _pt2);
    }
    // return the transverse energy squared
    inline double Et2() const {
        return (_pt2 <= EPSILON) ? 0.0 : _E * _E / (1.0 + _pz * _pz / _pt2);
    }
    // cos of the polar angle
    inline double cos_theta() const {
        return std::min(1.0, std::max(-1.0, _pz / sqrt(mod_p_2())));
    }
    // polar angle
    inline double theta() const { return acos(cos_theta()); }
    // get and set user index
    inline int user_index() const { return _user_index; }
    inline void set_user_index(int i) { _user_index = i; }

    // return component i, where E==0, X==1, Y==2, Z==3
    double operator()(int i) const;
    inline double operator[](int i) const { return (*this)(i); }

    // return pt-distance between this particle and another one
    double pt_distance_with(const Particle &other) const;
    // return squared cylinder (rap-phi) distance between this particle and
    // another
    double squared_distance(const Particle &other) const;
    // return the cylinder (rap-phi) distance between this particle and another
    // Delta_R = sqrt(Delta y^2 + Delta phi^2)
    inline double delta_R_with(const Particle &other) const {
        return sqrt(squared_distance(other));
    }
    // returns other.phi() - this.phi(), i.e. the phi distance to
    // other, constrained to be in range -pi .. pi
    double delta_phi_to(const Particle &other) const;
    // returns other.phi() - this.phi(), i.e. the phi distance to
    // other, constrained to be in range 0 .. pi
    double delta_phi_with(const Particle &other) const;

    // return a valarray containing the four-momentum (components 1-3
    // are 3-mom, component 0 is energy).
    std::valarray<double> four_mom() const;

    // transform this particle (given in the rest frame of prest) into a
    // particle in the lab frame
    Particle &boost(const Particle &prest);
    // transform this particle (given in lab) into a particle in the rest
    // frame of prest
    Particle &inverse_boost(const Particle &prest);

    // overload operator
    Particle &operator*=(double);
    Particle &operator/=(double);
    Particle &operator+=(const Particle &);
    Particle &operator-=(const Particle &);

    // reset the 4-momentum according th the supplied components
    inline void reset(double E, double px, double py, double pz) {
        _E = E;
        _px = px;
        _py = py;
        _pz = pz;
        _set_mass();
        _finish_init();
    }
    inline void reset(double E, double px, double py, double pz, double m) {
        _E = E;
        _px = px;
        _py = py;
        _pz = pz;
        _mass = m;
        _finish_init();
    }
    inline void reset(const Particle &mParticle) { (*this) = mParticle; }

    // return the 4-vector dot product with another particle
    inline double dot_product_with(const Particle &P) {
        return _E * P.e() - _px * P.px() - _py * P.py() - _pz * P.pz();
    }

    // print particle Info
    void print_particle_info();
    // ------------------- END OF MEMBER VARIABLES ---------------------

private:
    // ------------- Member Variables ----------------------------------
    double _E, _px, _py, _pz;
    double _mass;
    int _pdg_id;
    mutable double _phi, _rap;
    double _pt2;
    int _user_index;
    // ------------- END OF MEMBER VARIABLES ---------------------------

    // ------------- Member functions ----------------------------------
    // calculate phi, rap, pt2 based on the 4-momentum components
    inline void _finish_init();
    // set
    void _set_rap_phi() const;
    // set mass
    inline double _get_mass_square() const {
        return pow(_E, 2) - pow(_px, 2) - pow(_py, 2) - pow(_pz, 2);
    }
    inline void _set_mass() {
        double mass_square = _get_mass_square();
        _mass = mass_square < 0.0 ? 0.0 : sqrt(mass_square);
    }
    inline void _set_pt2() { _pt2 = pow(_px, 2) + pow(_py, 2); }
    inline void _ensure_valid_rap_phi() const {
        if (std::isnan(_phi)) _set_rap_phi();
    }
    // ------------- END OF MEMBER VARIABLES ---------------------------

};  // class Particle

/************************************************************************/
/*                Routines for basic binary operations                  */
/************************************************************************/
Particle operator+(const Particle &, const Particle &);
Particle operator-(const Particle &, const Particle &);
Particle operator*(double, const Particle &);
Particle operator*(const Particle &, double);
Particle operator/(const Particle &, double);
// -----------------------------------------------------------------------

double Particle::pseudorapidity() const {
    if (_px <= EPSILON && _py <= EPSILON) return INFVALUE;
    if (_pz <= EPSILON) return 0.0;
    double theta = atan(this->pt() / this->_pz);
    if (theta < 0) theta += PI;
    return -log(tan(theta / 2.0));
}

double Particle::operator()(int i) const {
    switch (i) {
        case 0:
            return e();
            break;
        case 1:
            return px();
            break;
        case 2:
            return py();
            break;
        case 3:
            return pz();
            break;

        default:
            std::string err = "Particle subscripting: bad index (" +
                              std::to_string(i) + ") don't belong to 0 to 3.";
            throw std::out_of_range(err);
            break;
    }
}

// return pt-distance between this particle and another one
double Particle::pt_distance_with(const Particle &other) const {
    double distance = std::min(this->_pt2, other._pt2);
    double dphi = fabs(this->phi() - other.phi());
    if (dphi > PI) dphi = TwoPI - dphi;
    double drap = this->rap() - other.rap();
    distance = distance * (dphi * dphi + drap * drap);
    return distance;
}

// return squared cylinder (rap-phi) distance between this particle and another
double Particle::squared_distance(const Particle &other) const {
    double dphi = fabs(this->phi() - other.phi());
    if (dphi > PI) dphi = TwoPI - dphi;
    double drap = this->rap() - other.rap();
    return (dphi * dphi + drap * drap);
}

// returns other.phi() - this.phi(), i.e. the phi distance to
// other, constrained to be in range -pi .. pi
double Particle::delta_phi_to(const Particle &other) const {
    double dphi = other.phi() - this->phi();
    if (dphi > PI) dphi -= TwoPI;
    if (dphi < -PI) dphi += TwoPI;
    return dphi;
}
// returns other.phi() - this.phi(), i.e. the phi distance to
// other, constrained to be in range 0 .. pi
double Particle::delta_phi_with(const Particle &other) const {
    double dphi = delta_phi_to(other);
    return (dphi >= 0.0 ? dphi : dphi + PI);
}

// return a valarray containing the four-momentum (components 1-3
// are 3-mom, component 0 is energy).
std::valarray<double> Particle::four_mom() const {
    std::valarray<double> mom(4);
    mom[0] = _E;
    mom[1] = _px;
    mom[2] = _py;
    mom[3] = _pz;
    return mom;
}

// transform this particle (given in the rest frame of prest) into a
// particle in the lab frame
Particle &Particle::boost(const Particle &prest) {
    if (fabs(prest.px()) <= EPSILON && fabs(prest.py()) <= EPSILON &&
        fabs(prest.pz()) <= EPSILON)
        return *this;

    double m_local = prest.m();
    assert(m_local > EPSILON);

    double pf4 = (px() * prest.px() + py() * prest.py() + pz() * prest.pz() +
                  e() * prest.e()) /
                 m_local;
    double fn = (pf4 + e()) / (prest.e() + m_local);
    _px += fn * prest.px();
    _py += fn * prest.py();
    _pz += fn * prest.pz();
    _E = pf4;

    _finish_init();
    return *this;
}

// transform this particle (given in lab) into a particle in the rest
// frame of prest
Particle &Particle::inverse_boost(const Particle &prest) {
    if (fabs(prest.px()) <= EPSILON && fabs(prest.py()) <= EPSILON &&
        fabs(prest.pz()) <= EPSILON)
        return *this;

    double m_local = prest.m();
    assert(m_local > EPSILON);

    double pf4 = (-px() * prest.px() - py() * prest.py() - pz() * prest.pz() +
                  e() * prest.e()) /
                 m_local;
    double fn = (pf4 + e()) / (prest.e() + m_local);
    _px -= fn * prest.px();
    _py -= fn * prest.py();
    _pz -= fn * prest.pz();
    _E = pf4;

    _finish_init();
    return *this;
}

// multiply the particle's momentum by the coefficient
Particle &Particle::operator*=(double coeff) {
    _ensure_valid_rap_phi();
    _px *= coeff;
    _py *= coeff;
    _pz *= coeff;
    _E *= coeff;
    _pt2 *= coeff * coeff;
    _mass *= coeff;
    return *this;
}

// divide the particle's momentum by the coefficient
Particle &Particle::operator/=(double coeff) {
    (*this) *= 1.0 / coeff;
    return *this;
}

// add the other particle's momentum to this particle
Particle &Particle::operator+=(const Particle &other) {
    _px += other._px;
    _py += other._py;
    _pz += other._pz;
    _E += other._E;
    _set_mass();
    _finish_init();
    return *this;
}

// subtract the other particle's momentum from this particle
Particle &Particle::operator-=(const Particle &other) {
    _px -= other._px;
    _py -= other._py;
    _pz -= other._pz;
    _E -= other._E;
    _set_mass();
    _finish_init();
    return *this;
}

// print particle Info
void Particle::print_particle_info() {
    cout << std::setiosflags(std::ios::left)
         << "--Particle \033[1;32mINFO:\033[0m "
         << "PDG ID: " << pdg_id() << '\n'
         << "                 "
         << "E:    " << std::setw(8) << e() << "      "
         << "px:   " << std::setw(8) << px() << "      "
         << "py:   " << std::setw(8) << py() << "      "
         << "pz:   " << std::setw(8) << pz() << '\n'
         << "                 "
         << "mass：" << std::setw(8) << m() << "      "
         << "m_c:  " << std::setw(8) << m_calc() << "      "
         << "m2:   " << std::setw(8) << m2() << "      "
         << "m2_c: " << std::setw(8) << m2_calc() << '\n'
         << "                 "
         << "mt:   " << std::setw(8) << mt() << "      "
         << "pt:   " << std::setw(8) << pt() << "      "
         << "phi:  " << std::setw(8) << phi() << "      "
         << "phi_s:" << std::setw(8) << phi_std() << '\n'
         << "                 "
         << "rap:  " << std::setw(8) << rap() << "      "
         << "prap: " << std::setw(8) << pseudorapidity() << "      "
         << "modp2:" << std::setw(8) << mod_p_2() << "      "
         << "Et2:   " << std::setw(8) << Et2() << '\n'
         << "                 " << '\n';
}

/*************** PRIVATE MEMBER FUNCTIONS DEFINE *******************/
// Define of Member functions
inline void Particle::_finish_init() {
    _set_pt2();
    _phi = NaNValue;
    _rap = NaNValue;
}

void Particle::_set_rap_phi() const {
    if (_pt2 < EPSILON)
        _phi = 0.0;
    else
        _phi = atan2(_py, _px);
    if (_phi < 0.0) _phi += TwoPI;
    if (_phi >= TwoPI) _phi -= TwoPI;
    if (fabs(_E - fabs(_pz)) < EPSILON && _pt2 < EPSILON) {
        if (_pz >= 0.0)
            _rap = INFINITY;
        else
            _rap = -INFINITY;
    } else {
        double effective_m2 = std::max(0.0, _get_mass_square());
        double E_plus_pz = _E + fabs(_pz);
        _rap = 0.5 * log((_pt2 + effective_m2) / (E_plus_pz * E_plus_pz));
        if (_pz > 0.0) _rap = -_rap;
    }
}

/************************************************************************/
/*                Routines for basic binary operations                  */
/************************************************************************/
// return sum of two particle
Particle operator+(const Particle &P1, const Particle &P2) {
    return Particle(P1.e() + P2.e(), P1.px() + P2.px(), P1.py() + P2.py(),
                    P1.pz() + P2.pz());
}
// return difference of two particle
Particle operator-(const Particle &P1, const Particle &P2) {
    return Particle(P1.e() - P2.e(), P1.px() - P2.px(), P1.py() - P2.py(),
                    P1.pz() - P2.pz());
}
Particle operator*(double coeff, const Particle &P) {
    Particle coeff_times_particle = P;
    coeff_times_particle *= coeff;
    return coeff_times_particle;
}
Particle operator*(const Particle &P, double coeff) { return coeff * P; }
Particle operator/(const Particle &P, double coeff) {
    return (1.0 / coeff) * P;
}
// -----------------------------------------------------------------------

}  // namespace iHepTools

#endif  // __IHEPTOOLS_PARTICLE_H__