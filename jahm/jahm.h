
#ifndef JAHM_H
#define JAHM_H

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>

namespace ja {
const int32_t ROUND_VALUE = 1000;
const int32_t  H_BIAS = 40000; //bias magnetic field

class JAHM
{
    /**
     * @brief optimal initial values for the secant iteration method
     */
    const double   PROP_0 = 0.15;
    const double   PROP_1 = 0.21;
    /**
     * @brief error values eps1 - 0.000001, eps2 - 0.00000001
     */
    const double   ERR_1  = 1E-06;
    const double   ERR_2  = 1E-08;

    int32_t sat_mag;       // M_{s} (saturation magnetization) A/m
    int32_t anhyst_mag;    // a (hape factor) A/m
    int32_t aver_enrg;     // k (inked to hysteresis losses) A/m
    double   alpha;        // alpha (inked to inter-domain coupling)
    double   revers_coeff; // c (eversibility coefficient)

    int32_t init_hyst_value = H_BIAS; // H_{m}
    std::vector<double> h_t; // vector of the initial steps values for calculate array of a magnetic strength values
    std::vector<double> h_spe;
    std::vector<double> h0_points;
    std::vector<double> h_intrm;

    std::vector<double> h_points; // vector of the magnetic field strength values
    std::vector<double> m1_vec;

    /**
     * @brief m_an - The vector of a anhysteretic magnetization values
     *               is determined by the modified Langevin formula
     */
    std::vector<double> m_an;
    /**
     * @brief dman_dh - The differential equation describes anhysteretic magnetization (M_{an})
     *                  as a function of magnetic field (H)
     */
    std::vector<double> dman_dh;
    /**
     * @brief mag_tot - The magnetic moment of a unit volume M
     *                  is determined weighted sum of anhysteretic magnetization M_anl
     *                  and irreversible magnetization M_irr component
     */
    std::vector<double> mag_tot;
    std::vector<double> field_tot;

    double m_bias = 0.;
    bool mag_bias_done = false;

    void jaInitialStepVector(const double &begin, const double &end, const double &step);
    template <typename T> bool jaIsMember(const T& a, const std::vector<T>& b_vec) const;
    template <typename T> int32_t jaSignDelta(const T& a) const;
    long double jaCoth(long double a) const;
    void jaSecantRoutine();
    void jaRK4Routine();

public:
    JAHM() = delete;
    JAHM(uint32_t m_s, uint32_t a, uint32_t k, double alpha, double c);
    JAHM(uint32_t m_s, uint32_t a, uint32_t k, double alpha, double c, double begin, double end, double step);
    JAHM(uint32_t m_s, uint32_t a, uint32_t k, double alpha, double c, std::vector<double>&& h_values);
    ~JAHM();
    void jaSetInitialMagneticFieldValue(uint32_t a);
    void jaSetMagneticFieldStrengthVector(std::vector<double>&& other);

    void jaTotalMagnetizeCalc(double begin, double end, double step);
    void jaTotalMagnetizeCalc();

    std::vector<double> jaGetH();
    std::vector<double> jaGetM();
};
}
#endif // JAHM_H
