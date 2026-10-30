#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "jls.hpp"

namespace DiracOperator {

//! anomalous moment in the magnetic field operator: <a||h1||b>
/*! @details
\f[ <a||h1||b> = a_e R \f]
\f[R = \int ( 2<k_a||S||k_b> f_af_b - 2<-k_a||S||-k_b> g_ag_b ) \, dr\f]
Units of |mu_B|.
*/

class h1 final : public TensorOperator {
public:
  h1(const Grid &gr)
    : TensorOperator(1, Parity::even, +0.0, {}, Realness::real, false),
      m_r(gr.r()) {}
  h1 &operator=(const h1 &) = delete;
  h1(const h1 &) = default;
  ~h1() = default;
  std::string name()  const override final { return std::string("h1"); }
  std::string units() const override final { return std::string("|mu_B|"); }

  double angularF(int, int) const override final {
	// free electron a_e for now (one-loop (Schwinger) value)
	return PhysConst::alpha / (2.0 * M_PI);
  }
  double angularCff(const int ka, const int kb) const override final {
	//if (ka != kb) {
      return 2.0 * Angular::S_kk(ka, kb);
    //} else {
    //  return 0.0;
    //}
  }
  double angularCgg(const int ka, const int kb) const override final {
	//if (ka != kb) {
      return (-2.0) * Angular::S_kk(-ka, -kb);
    //} else {
    //  return 2.0 * (Angular::S_kk(ka, ka) + Angular::S_kk(-ka, -ka));
    //}
  }
  double angularCfg(int, int) const override final { return 0.0; }
  double angularCgf(int, int) const override final { return 0.0; }

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<h1>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"", "No input options"}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<h1>(wf.grid());
  }

private:
  const std::vector<double> m_r; // store radial vector (copy)
};

} // namespace DiracOperator