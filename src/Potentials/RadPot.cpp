#include "RadPot.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp" //for AMM corrections (need core of wf)
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/FGRadPot.hpp"
#include "Wavefunction/Wavefunction.hpp" //for AMM corrections (need core of wf)
#include "qip/Vector.hpp"
#include <algorithm>
#include <cstdlib> //for std::abort
#include <vector>

namespace QED {

//==============================================================================
RadPot::Xl::Xl() : m_x{1.0} {}
RadPot::Xl::Xl(std::vector<double> x) : m_x{std::move(x)} {}

double RadPot::Xl::operator()(int l) const {
  if (m_x.empty())
    return 1.0;
  if (l >= (int)m_x.size())
    return m_x.back();
  return m_x.at(std::size_t(l));
}

//==============================================================================
RadPot::RadPot()
  : m_Z(0),
    m_rN(0.0),
    m_rcut(0.0),
    m_f({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}),
    m_xl{},
    print(true) {}

//==============================================================================
RadPot::RadPot(const std::vector<double> &r, double Z, double rN, double rcut,
               Scale f, Xl xl, bool tprint, bool do_readwrite,
               const std::string &label)
  : m_Z(Z), m_rN(rN), m_rcut(rcut), m_f(f), m_xl(xl), print(tprint) {

  bool read_ok = false;
  if (do_readwrite)
    read_ok = read_write(r, IO::FRW::RoW::read, label);
  if (!read_ok) {
    form_potentials(r);
    if (do_readwrite)
      read_write(r, IO::FRW::RoW::write, label);
  }
}

//==============================================================================
bool RadPot::read_write(const std::vector<double> &r, IO::FRW::RoW rw,
                        const std::string &label_x) {

  std::string label = "_";
  if (m_f.u != 0.0)
    label += 'u';
  if (m_f.h != 0.0)
    label += 'h';
  if (m_f.l != 0.0)
    label += 'l';
  if (m_f.m != 0.0)
    label += 'm';
  if (m_f.dm != 0.0)
    label += "DM";
  if (m_f.mV2 != 0.0)
    label += "mV2";
  if (m_f.wk != 0.0)
    label += "w";
  if (m_rN == 0.0)
    label += "_pt";
  const auto fname = AtomData::atomicSymbol(int(m_Z + 0.001)) + label +
                     (label_x == "" ? "" : "_" + label_x) + ".qed.abf";

  const auto readQ = rw == IO::FRW::read;

  if (readQ && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = !readQ ? "Writing to " : "Reading from ";
  std::cout << rw_str << "QED rad. pot. file: " << fname << "\n" << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // read/write info + Grid:
  std::vector<double> t_r;
  if (!readQ) {
    t_r = r; // copy; cannot write const ref to file?
  }
  rw_binary(iofs, rw, m_rN, m_rcut, t_r);
  if (readQ) {
    std::cout << "(rcut = " << m_rcut << "au, rN=" << m_rN * PhysConst::aB_fm
              << " fm)\n";
  }
  // read/write potentials::
  rw_binary(iofs, rw, mVu, mVh, mVl, mHm, mVwk);

  const auto grid_same =
    (r.size() == t_r.size() && (qip::compare_eps(r, t_r).first < 1.0e-3));

  if (readQ && !grid_same) {
    std::cout << "Interpolating QED rad-pot onto current grid.\n";

    if (!mVu.empty())
      mVu = Interpolator::interpolate(t_r, mVu, r);
    if (!mVh.empty())
      mVh = Interpolator::interpolate(t_r, mVh, r);
    if (!mVl.empty())
      mVl = Interpolator::interpolate(t_r, mVl, r);
    if (!mHm.empty())
      mHm = Interpolator::interpolate(t_r, mHm, r);
    if (!mDeltaH.empty())
      mDeltaH = Interpolator::interpolate(t_r, mDeltaH, r);
    if (!mHmV2.empty())
      mHmV2 = Interpolator::interpolate(t_r, mHmV2, r);
    if (!mVwk.empty())
      mVwk = Interpolator::interpolate(t_r, mVwk, r);
  }

  return true;
}

//==============================================================================
void RadPot::form_potentials(const std::vector<double> &r) {
  if (print) {
    std::cout << "Forming QED radiative potential ";
    if (m_rcut > 0.0) {
      std::cout << "(rcut = " << m_rcut << ") ";
    }
    std::cout << "with Rn = " << m_rN * PhysConst::aB_fm << "fm\n";
  }

  if (m_f.u != 0.0) {
    if (print)
      std::cout << "Uehling; scale=" << m_f.u << "\n";
    mVu = fill(FGRP::V_Uehling, r);
  }
  if (m_f.h != 0.0) {
    if (print)
      std::cout << "Self-energy (high freq); scale=" << m_f.h << "\n";
    const auto stride = std::max(r.size() / 400, 2ul);
    // Stride of 5 with 2000 points
    mVh = fill(FGRP::V_SEh, r, stride);
  }
  if (m_f.l != 0.0) {
    if (print)
      std::cout << "Self-energy (low freq); scale=" << m_f.l << "\n";
    mVl = fill(FGRP::V_SEl, r);
  }
  if (m_f.m != 0.0) {
    if (print)
      std::cout << "Self-energy (magnetic); scale=" << m_f.m << "\n";
    mHm = fill(FGRP::H_Magnetic, r);
  }
  if (m_f.dm != 0.0) {
    if (print)
      std::cout << "Self-energy (magnetic) DeltaH (Overleaf 5.2); scale="
                << m_f.dm << "\n";
  }
  if (m_f.mV2 != 0.0) {
    if (print)
      std::cout << "Self-energy (magnetic) V2 (Overleaf 5.1); scale=" << m_f.mV2
                << "\n";
  }
  if (m_f.wk != 0.0) {
    if (print)
      std::cout << "Wickman-Kroll; scale=" << m_f.wk << "\n";
    mVwk = fill(FGRP::V_WK, r);
  }
}

//==============================================================================
std::vector<double> RadPot::Vel(int l) const {
  using namespace qip::overloads;
  const auto a = FGRP::Fit::A(m_Z, l);
  const auto b = FGRP::Fit::B(m_Z, l);
  const auto xl = m_xl(l);
  return xl * ((m_f.u * mVu) + (a * m_f.h * mVh) + (b * m_f.l * mVl) +
               (m_f.wk * mVwk));
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Hmag(int l) const {
  using namespace qip::overloads;
  const auto xl = m_xl(l);
  return xl * m_f.m * mHm;
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Hmag_AMM(int l) const {
  using namespace qip::overloads;
  const auto xl = m_xl(l);

  if (m_f.mV2 != 0.0) {
    return xl * m_f.mV2 * mHmV2;
  }

  return xl * (m_f.m * mHm + m_f.dm * mDeltaH);
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Vu(int l) const {
  using namespace qip::overloads;
  const auto xl = m_xl(l);
  return xl * m_f.u * mVu;
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Vl(int l) const {
  using namespace qip::overloads;
  const auto b = FGRP::Fit::B(m_Z, l);
  const auto xl = m_xl(l);
  return xl * b * m_f.l * mVl;
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Vh(int l) const {
  using namespace qip::overloads;
  const auto a = FGRP::Fit::A(m_Z, l);
  const auto xl = m_xl(l);
  return xl * a * m_f.h * mVh;
}

void RadPot::update_core_corrections(const Wavefunction &wf) {

  // Build Ncore_r by integrating each core state's f and g components.
  // Build mDeltaH from Ncore_r.
  // Build nuclear enclosed charge Z_r.
  // Build mHmV2 from Z_r - Ncore_r.

  const double scale_dm = m_f.dm;
  const double scale_mV2 = m_f.mV2;

  if ((scale_dm != 0.0) && (scale_mV2 != 0.0)) {
    std::cout << "dm and mV2 can't be non-zero at the same time!" << "\n";
  } else if (scale_dm != 0.0) {
    mDeltaH = form_DeltaH(wf, scale_dm);
  } else if (scale_mV2 != 0.0) {
    const std::vector<double> DeltaH =
      form_DeltaH(wf, 1.0); // scale both at same time at end
    const std::vector<double> Zr_r2 = form_Zr_r2(wf);

    // check if both vectors same length
    // if not then error
    if (DeltaH.size() != Zr_r2.size()) {
      std::cout << "QED AMM correction vectors have different sizes" << "\n";
      std::abort();
    }
    //else calculate DeltaH + Zr_r2
    mHmV2.resize(DeltaH.size());
    for (std::size_t i = 0; i < DeltaH.size(); ++i) {
      mHmV2[i] = scale_mV2 * (DeltaH[i] + Zr_r2[i]); //probably overload qip::overloads
    }
  }
}
//==============================================================================
std::vector<double> RadPot::form_DeltaH(const Wavefunction &wf,
                                        const double scale) const {
  const Grid &gr = wf.grid();
  std::vector<double> Ncore_r(gr.num_points(), 0.0);
  for (const auto &Fc : wf.core()) {

    const std::vector<double> fac_Ncore(gr.num_points(), Fc.twojp1());

    NumCalc::additivePIntegral<NumCalc::zero_to_r>(Ncore_r, fac_Ncore, Fc.f(),
                                                   Fc.f(), gr);
    NumCalc::additivePIntegral<NumCalc::zero_to_r>(Ncore_r, fac_Ncore, Fc.g(),
                                                   Fc.g(), gr);
  } // end for loop over core states

  //double a_e = PhysConst::alpha / (2.0 * M_PI);
  //double fac = a_e * PhysConst::alpha / 2.0;
  double global_fac = PhysConst::alpha2 / (4.0 * M_PI);
  std::vector<double> DeltaH(gr.num_points(), 0.0);
  for (std::size_t i = 0; i < gr.num_points(); ++i) {
    DeltaH[i] = -scale * global_fac * Ncore_r[i] / (gr.r(i) * gr.r(i)); //unsure if numerically stable
    // ask Ben probably
  }
  return DeltaH;
}
//==============================================================================
std::vector<double> RadPot::form_Zr_r2(const Wavefunction &wf) const {

  //need to get Z, \rho_nuc(r), and core f(r)/g(r)/j
  const Grid &gr = wf.grid();
  const std::vector<double> rho_nuc = nuclear_density(wf);

  std::vector<double> Zr(gr.num_points(), 0.0); // enclosed charge
  const std::vector<double> fac_Z(gr.num_points(),
                                  4.0 * M_PI * wf.nucleus().z());
  std::vector<double> r2(gr.num_points(), 0.0);
  for (std::size_t i = 0; i < gr.num_points(); ++i) {
    r2[i] = gr.r(i) * gr.r(i);
  }

  NumCalc::additivePIntegral<NumCalc::zero_to_r>(Zr, fac_Z, r2, rho_nuc, gr);

  //return full result, is bc void NumCalc::additivePIntegral adds to 1st argument
  double global_fac = PhysConst::alpha2 / (4.0 * M_PI);
  std::vector<double> Zr_r2(gr.num_points(), 0.0);
  for (std::size_t i = 0; i < gr.num_points(); ++i) {
    Zr_r2[i] = global_fac * Zr[i] / (gr.r(i) * gr.r(i));
  }

  return Zr_r2;
}
//==============================================================================
std::vector<double> RadPot::nuclear_density(const Wavefunction &wf) const {

  const auto &gr = wf.grid();
  const auto &nucleus = wf.nucleus();

  const double r_rms = nucleus.r_rms() / PhysConst::aB_fm; // convert to au
  const auto nucleus_type =
    (r_rms > 0.0 || nucleus.type() == Nuclear::ChargeDistro::custom) ?
      nucleus.type() :
      Nuclear::ChargeDistro::point;

  switch (nucleus_type) { // Fermi, spherical, point, Gaussian, custom, Error

  case Nuclear::ChargeDistro::Fermi: {
    return Nuclear::fermiNuclearDensity_tcN(nucleus.t(), nucleus.c(), 1.0, gr);
  }
  case Nuclear::ChargeDistro::spherical: {
    const double rN = sqrt(5.0 / 3.0) * r_rms;
    const double rho_inside = 3.0 / (4.0 * M_PI * pow(rN, 3.0));

    std::vector<double> rho(gr.num_points(), 0.0);
    for (std::size_t i = 0; i < gr.num_points(); ++i) {
      if (gr.r(i) < rN) {
        rho[i] = rho_inside;
      }
    }
    return rho;
  }
  case Nuclear::ChargeDistro::Gaussian: {
    std::vector<double> rho(gr.num_points(), 0.0);

    const double factor = std::pow(3.0 / (2.0 * M_PI * r_rms * r_rms), 1.5);

    for (std::size_t i = 0; i < gr.num_points(); ++i) {
      const double r = gr.r(i);
      rho[i] = factor * std::exp((-3.0) * r * r / (2.0 * r_rms * r_rms));
    }
    return rho;
  }
  case Nuclear::ChargeDistro::point: {
    std::cout << "\nCannot do AMM correction for case "
              << Nuclear::parseType(nucleus_type) << "\n";
    std::abort();
  }
  case Nuclear::ChargeDistro::custom: {
    std::cout << "\nCannot do AMM correction for case "
              << Nuclear::parseType(nucleus_type) << "\n";
    std::abort();
  }
  case Nuclear::ChargeDistro::Error: {
    std::cout << "\nCannot do AMM correction for case "
              << Nuclear::parseType(nucleus_type) << "\n";
    std::abort();
  }
  }
}

//==============================================================================
RadPot ConstructRadPot(const std::vector<double> &r, double Z_eff, double rN_au,
                       const IO::InputBlock &input, bool print,
                       bool do_readwrite, const std::string &label) {

  input.check(
    {{"", "QED Radiative potential will be included if this block is present"},
     {"", "The following 5 are all doubles. Scale to include * potential; "
          "usually either 0.0 or 1.0, but can take any value:"},
     {"Ueh", "  Uehling (vacuum polarisation). [1.0]"},
     {"SE", "  Self-energy. [1.0]"},
     {"", "  The following are for indevidual components of the self-energy "
          "tern; will over-write the above:"},
     {"SE_h", "    self-energy high-freq electric."},
     {"SE_l", "    self-energy low-freq electric."},
     {"SE_m", "    self-energy magnetic."},
     {"SE_dm", "    self-energy magnetic DeltaH (Overleaf 5.2)."},
     {"SE_mV2", "    self-energy magnetic V2 (Overleaf 5.1)."},
     {"", ""},
     {"WK",
      "  Wickman-Kroll (approximate form of higher-order vac. pol.). [0.0]"},
     {"rcut", "Maximum radius (au) to calculate Rad Pot for [5.0]"},
     {"scale_rN", "Scale factor for Nuclear size. 0 for pointlike, 1 for "
                  "typical [1.0]"},
     {"scale_l", "List of doubles. Extra scaling factor for each l e.g., "
                 "1,0,1 => include for s and d, but not for p [1.0]"}});
  if (input.has_option("help")) {
    return {};
  }

  const auto x_Ueh = input.get("Ueh", 1.0);
  const auto x_SE = input.get("SE", 1.0);
  const auto x_SEe_h = input.get("SE_h", x_SE);
  const auto x_SEe_l = input.get("SE_l", x_SE);
  const auto x_SEm = input.get("SE_m", x_SE);
  const auto x_SEdm = input.get("SE_dm", x_SE);
  const auto x_SEmV2 = input.get("SE_mV2", x_SE);
  const auto x_wk = input.get("WK", 0.0);
  const auto rcut = input.get("rcut", 5.0);
  const auto scale_rN = input.get("scale_rN", 1.0);
  const auto x_spd = input.get("scale_l", std::vector{1.0});

  return RadPot(r, Z_eff, scale_rN * rN_au, rcut,
                {x_Ueh, x_SEe_h, x_SEe_l, x_SEm, x_SEdm, x_SEmV2, x_wk}, x_spd,
                print, do_readwrite, label);
}

} // namespace QED
