#include "jahm2.h"

double ja::JAHM2::coth(double x)
{
    return 1.0 / std::tanh(x);
}

std::vector<size_t> ja::JAHM2::find_special_indices(const std::vector<double> &tH)
{
    std::vector<size_t> Hspe;
    for (size_t i = 0; i < tH.size(); ++i) {
        if (std::abs(std::remainder(tH[i], 0.5)) < 1e-9) {
            Hspe.push_back(i);
        }
    }
    return Hspe;
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
void ja::JAHM2::secant_method(const std::vector<double> &Hpvec, std::vector<double> &Man, std::vector<double> &dMandHvec, const std::vector<size_t> &H0point)
{
    for (size_t iH = 0; iH < Hpvec.size(); ++iH) {
        double Hp = Hpvec[iH];
        if (std::find(H0point.begin(), H0point.end(), iH) != H0point.end()) {
            Man[iH] = 0.0;
            dMandHvec[iH] = 1.0 / ((3.0 * a / Ms) - alpha); // Eq.(11)
        } else {
            double Man_sec0 = Ms / a * Hp * prop0; //Eq.(8)
            double Man_sec1 = Ms / a * Hp * prop1; //Eq.(8)
            double fMan_sec0 = Man_sec0 - Ms * (coth((Hp + alpha * Man_sec0) / a) - a / (Hp + alpha * Man_sec0)); //Eq.(6)
            double fMan_sec1 = Man_sec1 - Ms * (coth((Hp + alpha * Man_sec1) / a) - a / (Hp + alpha * Man_sec1)); //Eq.(6)

            int isec = 1;
            int max_iter = 100;
            while (std::abs(fMan_sec1) >= eps * std::abs(Man_sec1) && isec < max_iter) {

                /* secant method formula
                 *
                 * x_{i+1} = x_{i} - (f(x_{i}) * (x_{i-1} - x_{i}))/(f(x_{i-1}) - f(x_{i}))
                 *
                 * */
                double Man_sec2 = Man_sec1 - fMan_sec1 * (Man_sec1 - Man_sec0) / (fMan_sec1 - fMan_sec0);
                Man_sec0 = Man_sec1;
                fMan_sec0 = fMan_sec1;
                Man_sec1 = Man_sec2;
                fMan_sec1 = Man_sec1 - Ms * (coth((Hp + alpha * Man_sec1) / a) - a / (Hp + alpha * Man_sec1)); // Eq. (5)
                isec++;
            }
            if (isec == max_iter) {
                std::cerr << "Secant method did not converge!\n";
            }
            Man[iH] = Man_sec1; // M_{an}
            double T = (Hp + alpha * Man_sec1) / a;
            double Tsin2 = std::pow(std::sinh(T), 2);
            double T2 = std::pow(T, 2);
            double denominator = (a * Tsin2 * T2 / Ms) - (alpha * (Tsin2 - T2));

            // Check for division by zero
            if (std::abs(denominator) < 1e-9) {
                denominator = 1e-9;  // or another smal value
            }

            dMandHvec[iH] = (Tsin2 - T2) / denominator;
        }
    }
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
void ja::JAHM2::rk4_method(const std::vector<double> &Hvec, const std::vector<double> &Manvec, const std::vector<double> &dMandHvec, std::vector<double> &M1vec, double Mbias) const
{
    for (size_t n = 1; n < Hvec.size(); ++n) {
                double delta = (Hvec[n] > Hvec[n-1]) ? 1.0 : -1.0;
                double h = Hvec[n] - Hvec[n-1];
                double Man = Manvec[n];
                double dMandH = dMandHvec[n-1];
                double M1 = M1vec[n-1];
                double delM = ((delta * (Man - M1) > 0) ? 1.0 : -1.0);
                double phi1 = (1.0 - c) * delM * (Man - M1);

                // Runge-Kutta methods coefficients
                double k1 = (delta * k * c * dMandH + phi1) / (delta * k - alpha * phi1);

                M1 = M1vec[n-1] + h * k1 / 2.0;
                delM = ((delta * (Man - M1) > 0) ? 1.0 : -1.0);
                phi1 = (1.0 - c) * delM * (Man - M1);
                double k2 = (delta * k * c * dMandH + phi1) / (delta * k - alpha * phi1);

                M1 = M1vec[n-1] + h * k2 / 2.0;
                delM = ((delta * (Man - M1) > 0) ? 1.0 : -1.0);
                phi1 = (1.0 - c) * delM * (Man - M1);
                double k3 = (delta * k * c * dMandH + phi1) / (delta * k - alpha * phi1);

                M1 = M1vec[n-1] + h * k3;
                delM = ((delta * (Man - M1) > 0) ? 1.0 : -1.0);
                phi1 = (1.0 - c) * delM * (Man - M1);
                double k4 = (delta * k * c * dMandH + phi1) / (delta * k - alpha * phi1);

                M1vec[n] = M1vec[n-1] + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
}

void ja::JAHM2::calculate()
{
    // Initialization of tH
    for (double t = 0.0; t <= 40.0; t += tin) {
        tH.push_back(t);
    }

    // Hpvec initialization
    Hpvec.resize(tH.size());
    for (size_t i = 0; i < tH.size(); ++i) {
        Hpvec[i] = Hbias * std::sin(M_PI * tH[i]);
    }

    // Search for special indexes
    Hspe = find_special_indices(tH);
    for (size_t i = 0; i < Hspe.size(); i += 2) {
        H0point.push_back(Hspe[i]);
    }

    // Initializing Man and dMandHvec
    Man.resize(Hpvec.size(), 0.0);
    dMandHvec.resize(Hpvec.size(), 0.0);

    // Secant method
    secant_method(Hpvec, Man, dMandHvec, H0point);

    // M1vec initialization
    M1vec.resize(Hpvec.size(), 0.0);

    // Runge-Kutta method
    rk4_method(Hpvec, Man, dMandHvec, M1vec, 0.0);

    // Search Mbias
    Mbias = *std::max_element(M1vec.begin(), M1vec.end());

    // Main loop for each Hm
    for (size_t iHmax = 0; iHmax < Hmvec.size(); ++iHmax) {
        double Hm = Hmvec[iHmax];
        HACvec.resize(tH.size());
        for (size_t i = 0; i < tH.size(); ++i) {
            HACvec[i] = Hm * std::sin(M_PI * tH[i]);
        }
        std::vector<double> Hpvec_new(tH.size());
        for (size_t i = 0; i < tH.size(); ++i) {
            Hpvec_new[i] = HACvec[i] + Hbias;
        }

        // Recalculating Man and dMandHvec
        std::vector<double> Man_new(Hpvec_new.size(), 0.0);
        std::vector<double> dMandHvec_new(Hpvec_new.size(), 0.0);
        secant_method(Hpvec_new, Man_new, dMandHvec_new, H0point);

        // Recalculating M1vec
        std::vector<double> M1vec_new(Hpvec_new.size(), Mbias);
        rk4_method(Hpvec_new, Man_new, dMandHvec_new, M1vec_new, Mbias);

        // Extraction of a stationary section
        size_t numst_start = static_cast<size_t>(std::floor(Hpvec_new.size() * 0.9));
        std::vector<double> tst(tH.begin() + numst_start, tH.end());
        Hst.assign(Hpvec_new.begin() + numst_start, Hpvec_new.end());
        std::vector<double> HACst(HACvec.begin() + numst_start, HACvec.end());
        Mst.assign(M1vec_new.begin() + numst_start, M1vec_new.end());
        for (auto& val : Mst) {
            val -= Mbias;
        }

        // Calculation B
        Bst.resize(Hst.size());
        for (size_t i = 0; i < Hst.size(); ++i) {
            Bst[i] = mu0 * (Hst[i] + Mst[i]);
        }
    }
}
