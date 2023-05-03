#include "jahm.h"

/**
 * @brief ja::JAHM::jaInitialStepVector initialize vector with values of steps
 * @param begin - start value
 * @param end - finish value
 * @param step - step value
 */
void ja::JAHM::jaInitialStepVector(const double &begin, const double &end, const double &step)
{
    double k = begin;
    int32_t size = (end-begin)/step+1;
    int32_t count = 0;
    h_t.resize(size, 0);
    while(k < end)
    {
        h_t.at(count) = std::round(k * ROUND_VALUE) / ROUND_VALUE;
        k += step;
        ++count;
    }
}

/**
 * @brief ja::JAHM::jaCoth - 1/th(x) - hyperbolic cotangent
 * @param a - input value
 * @return hyperbolic cotangent value
 */
inline long double ja::JAHM::jaCoth(long double a) const
{
    return 1/std::tanh(a);
}

/**
 * @brief ja::JAHM::jaIsMember - check contain the current value in the array
 * @param a - test value
 * @param b_vec - array
 * @return boolean variable
 */
template<typename T = double>
inline bool ja::JAHM::jaIsMember(const T &a, const std::vector<T> &b_vec) const
{
    if(std::find(std::begin(b_vec), std::end(b_vec), a) != std::end(b_vec)){ return true; }
    return false;
}

/**
 * @brief ja::JAHM::jaSignDelta - check of sign the variable
 * @param a - input value
 * @return "sign" 1 or -1
 */
template<typename T = double>
inline int32_t ja::JAHM::jaSignDelta(const T &a) const
{
    return (T(0) < a) - (a < T(0));
}

/**
 * @brief ja::JAHM::JAHM - constructor
 * @param m_s - Saturation magnetization
 * @param a - Anhysteretic shape coefficient
 * @param k - Irreversible loss coefficient
 * @param alpha - Internal coupling parameter
 * @param c - Reversible coefficient
 */
ja::JAHM::JAHM(uint32_t m_s,
               uint32_t a,
               uint32_t k,
               double alpha,
               double c)
    :sat_mag(m_s)
    ,anhyst_mag(a)
    ,aver_enrg(k)
    ,alpha(alpha)
    ,revers_coeff(c)
{}

/**
 * @brief ja::JAHM::JAHM - constructor
 * @param m_s - Saturation magnetization
 * @param a - Anhysteretic shape coefficient
 * @param k - Irreversible loss coefficient
 * @param alpha - Internal coupling parameter
 * @param c - Reversible coefficient
 * @param begin - vector of steps, start value
 * @param end - finish value
 * @param step - step value
 */
ja::JAHM::JAHM(uint32_t m_s,
               uint32_t a,
               uint32_t k,
               double alpha,
               double c,
               double begin,
               double end,
               double step)
    :sat_mag(m_s)
    ,anhyst_mag(a)
    ,aver_enrg(k)
    ,alpha(alpha)
    ,revers_coeff(c)
{
    jaInitialStepVector(begin, end, step);
}

/**
 * @brief ja::JAHM::JAHM - constructor
 * @param m_s - Saturation magnetization
 * @param a - Anhysteretic shape coefficient
 * @param k - Irreversible loss coefficient
 * @param alpha - Internal coupling parameter
 * @param c - Reversible coefficient
 * @param h_values - using input vector with magnetic field values
 */
ja::JAHM::JAHM(uint32_t m_s,
               uint32_t a,
               uint32_t k,
               double alpha,
               double c,
               std::vector<double> &&h_values)
    :sat_mag(m_s)
    ,anhyst_mag(a)
    ,aver_enrg(k)
    ,alpha(alpha)
    ,revers_coeff(c),
    h_points(std::move(h_values))

{}

ja::JAHM::~JAHM()
{}

/**
 * @brief ja::JAHM::jaSetInitialMagneticFieldValue
 * @param a - initial hysteresis
 */
void ja::JAHM::jaSetInitialMagneticFieldValue(uint32_t a)
{
    init_hyst_value = a;
}

/**
 * @brief ja::JAHM::jaSetMagneticFieldStrengthVector - swap with other h_points vector
 * @param other - vector with magnetic field values
 */
void ja::JAHM::jaSetMagneticFieldStrengthVector(std::vector<double> &&other)
{
    h_points.resize(other.size(), 0);
    std::swap(h_points, other);
}

/**
 * @brief ja::JAHM::jaTotalMagnetizeCalc
 * @param begin - vector of steps, start value
 * @param end - finish value
 * @param step - step value
 */
void ja::JAHM::jaTotalMagnetizeCalc(double begin, double end, double step)
{

    jaInitialStepVector(begin, end, step);

    for(uint32_t k=0; k<h_t.size(); ++k)
    {
        auto tmp = h_t[k] * ROUND_VALUE;
        auto var = std::floor(tmp) / static_cast<double>(ROUND_VALUE);
        auto result = (std::abs(std::remainder(var, 0.5))) == 0;
        if(result){ h_spe.push_back(k); }
    }

    for(uint32_t i=1; i<h_spe.size(); i+=2)
    {
        h0_points.push_back(h_spe[i-1]);
    }

    h_points.resize(h_t.size(), 0);
    for(uint32_t i=0; i<h_t.size(); ++i)
    {
        h_points[i] = (init_hyst_value * std::sin(M_PI * h_t[i]));
    }

    jaSecantRoutine();
    jaRK4Routine();
    mag_bias_done = true;
    h_intrm.resize(h_t.size(), 0);

    for(uint32_t i=0; i<h_t.size(); ++i)
    {
        h_points[i] = init_hyst_value * std::sin(M_PI * h_t[i]);
        h_intrm[i] = h_points[i];
        h_points[i] += H_BIAS;
    }

    jaSecantRoutine();
    jaRK4Routine();
}

/**
 * @brief ja::JAHM::jaGetH
 * @return vector of the magnetic field values
 */
std::vector<double> ja::JAHM::jaGetH()
{
    for(uint32_t i=0; i<field_tot.size(); ++i)
    {
        field_tot[i] = std::floor(field_tot[i]);
    }
    return field_tot;
}

/**
 * @brief ja::JAHM::jaGetM
 * @return vector of the total magnetization
 */
std::vector<double> ja::JAHM::jaGetM()
{
    for(uint32_t i=0; i<field_tot.size(); ++i)
    {
        mag_tot[i] = std::floor(mag_tot[i]);
    }
    return mag_tot;
}

/**
 * @brief ja::JAHM::jaSecantRoutine - estimate anhysteretic magnetization vector(M_{an})
 *        and defivative anhysteretic magnetization to magnetic field strength(dM_{an}/dH)
 *        The tangent method or secant method has the advantages of fast convergence.
 *        Considering the derivative information of M_an was not given,
 *        the tangent method was not suitable, and the secant method was employed
 *        to reduce iterative numbers.
 *        See for equation descript Xue G. - Numerical Solving Method for Jiles-Atherton Model and Influence Analysis of the Initial Magnetic Field on Hysteresis
 */
void ja::JAHM::jaSecantRoutine()
{
    m_an.resize(h_points.size(), 0);
    dman_dh.resize(h_points.size(), 0);

    auto check_ = [&](bool _mag_bias_done, double var, uint32_t index) -> bool
    {
        if(_mag_bias_done){
            return var <= ERR_2;
        }else{
            return jaIsMember<double>(static_cast<double>(index), h0_points);
        }
    };

    for(uint32_t i=0; i<h_points.size(); ++i)
    {
        auto h_p = h_points[i];
        bool tmp = false;
        tmp = check_(mag_bias_done, h_p, i);

        if(tmp)
        {
            m_an[i] = 0;
            dman_dh[i] = 1./(((3. * anhyst_mag)/sat_mag) - alpha); // Eq.(11)
        }
        else
        {
            auto m_an_sec0 = sat_mag / anhyst_mag * h_p * PROP_0; //Eq.(8)
            auto m_an_sec1 = sat_mag / anhyst_mag * h_p * PROP_1; //Eq.(8)
            auto fm_an_sec0 = m_an_sec0 - sat_mag * (jaCoth((h_p + alpha * m_an_sec0) / anhyst_mag) - (anhyst_mag / (h_p + alpha * m_an_sec0))); //Eq.(6)
            auto fm_an_sec1 = m_an_sec1 - sat_mag * (jaCoth((h_p + alpha * m_an_sec1) / anhyst_mag) - (anhyst_mag / (h_p + alpha * m_an_sec1))); //Eq.(6)
            auto i_sec = 1;
            while(std::abs(fm_an_sec1) >= ERR_1 * std::abs(m_an_sec1))
            {
                /* secant method formula
                 *
                 * x_{i+1} = x_{i} - (f(x_{i}) * (x_{i-1} - x_{i}))/(f(x_{i-1}) - f(x_{i}))
                 *
                 * */
                auto m_an_sec2 = m_an_sec1 - (fm_an_sec1 * (m_an_sec1 - m_an_sec0) / (fm_an_sec1 - fm_an_sec0));
                m_an_sec0 = m_an_sec1;
                fm_an_sec0 = fm_an_sec1;
                m_an_sec1 = m_an_sec2;
                fm_an_sec1 = m_an_sec1 - sat_mag * (jaCoth((h_p + alpha * m_an_sec1) / anhyst_mag) - (anhyst_mag / (h_p + alpha * m_an_sec1))); // Eq. (5)
                i_sec++;
            }
            m_an[i] = m_an_sec1; // M_{an}
            auto t = ((h_p + alpha * m_an_sec1) / anhyst_mag);
            auto t_sin2 = std::pow(std::sinh(t), 2);
            auto t_2 = std::pow(t, 2);
            dman_dh[i] = (t_sin2 - t_2) / (anhyst_mag * t_sin2 * t_2 / sat_mag - alpha * (t_sin2 - t_2)); // Appendix (A1) Equation - (dM_{an})/(dH)
        }
    }
    m_an[h_points.size()-1] = m_an[0];
    dman_dh[h_points.size()-1] = dman_dh[0];
}

/**
 * @brief ja::JAHM::jaRK4Routine - solve the total magnetization
 *        To solve the nonlinear ordinary differential equation
 *        the 4th order Runge-Kutta method could be employed. Imposing
 *        phi(H,M) = (delta*k*c*(dM_{an}/dH) + delta_M(1 - c)(M_{an} - M))
 *                                            /
 *                   (delta*k - alpha*delta_M(1 - c)(M_{an} - M))
 *        the computing approach is
 *        M_{n+1} = M_n + (h/6)(k-1 + 2k_2 + 2k_3 + k_4)
 *        k_1 = phi(H_m, M_n)
 *        k_2 = phi(H_m + (h/2), M_n + (hk_1/2))
 *        k_3 = phi(H_m + (h/2), M_n + (hk_2/2))
 *        k_4 = phi(H_m + h, M_n + hk_3)
 *        See for more Xue G. - Modification and Numerical Method fo the Jiles-Atherton Hysteresis Model
 *
 */
void ja::JAHM::jaRK4Routine()
{
    uint32_t p_num = h_t.size();
    m1_vec.resize(p_num, 0);
    std::vector<double> h_vec = h_points;

    auto max_ = [&](const std::vector<double>& arr) -> int32_t
    {
        std::vector<double> tmp;
        std::copy(std::begin(arr), std::end(arr), std::back_inserter(tmp));
        std::sort(std::begin(tmp), std::end(tmp));
        return tmp[tmp.size()-1];
    };

    if(mag_bias_done){
        for(uint32_t i = 0; i < p_num; ++i)
        {
            m1_vec[i] = 0;
            m1_vec[i] += m_bias;
        }
    }

    for(uint32_t i=1; i<p_num; ++i)
    {
        auto delta = jaSignDelta<>(h_vec[i] - h_vec[i-1]);
        auto h = (h_vec[i] - h_vec[i-1]);
        auto man_ = m_an[i];
        auto dman_dh_ = dman_dh[i-1];
        auto M1 = m1_vec[i-1];
        auto delM = jaSignDelta<>(jaSignDelta<>(delta * (man_ - M1)) + 1);
        auto phi1 = (1 - revers_coeff) * delM * (man_ - M1);

        double k1 = ((delta * aver_enrg * revers_coeff * dman_dh_) + phi1)/((delta * aver_enrg) - (alpha * phi1));

        M1 = (m1_vec[i-1] + h * k1 / 2);
        delM = jaSignDelta<>(jaSignDelta<>(delta * (man_ - M1)) + 1);
        phi1 = (1 - revers_coeff) * delM * (man_ - M1);

        double k2 = ((delta * aver_enrg * revers_coeff * dman_dh_) + phi1)/((delta * aver_enrg) - (alpha * phi1));

        M1 = (m1_vec[i-1] + h * k2 / 2);
        delM = jaSignDelta<>(jaSignDelta<>(delta * (man_ - M1)) + 1);
        phi1 = (1 - revers_coeff) * delM * (man_ - M1);

        double k3 = ((delta * aver_enrg * revers_coeff * dman_dh_) + phi1)/((delta * aver_enrg) - (alpha * phi1));

        M1 = (m1_vec[i-1] + h * k3);
        delM = jaSignDelta<>(jaSignDelta<>(delta * (man_ - M1)) + 1);
        phi1 = (1 - revers_coeff) * delM * (man_ - M1);

        double k4 = ((delta * aver_enrg * revers_coeff * dman_dh_) + phi1)/((delta * aver_enrg) - (alpha * phi1));
        m1_vec[i] = (m1_vec[i-1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4));
    }

    if(mag_bias_done)
    {
        std::size_t size = h_points.size() - std::trunc(h_points.size() * 0.9);
        mag_tot.resize(size, 0);
        field_tot.resize(size, 0);
        for(uint32_t i=std::trunc(h_points.size()*0.9)-1, k=0; i<h_points.size() && k<size; ++i, ++k)
        {
            mag_tot[k] = m1_vec[i];
            field_tot[k] = h_intrm[i];
        }

        for(uint32_t k=0; k< mag_tot.size(); ++k)
        {
            mag_tot[k] -= m_bias;
        }
    }
    else
    {
        m_bias = max_(m1_vec);
    }
}
