#ifndef JAHM2_H
#define JAHM2_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <fstream>

namespace ja {
// Auxiliary constants
constexpr double prop0 = 0.15;
constexpr double prop1 = 0.21;
// Magnetic constant, H/m
constexpr double mu0 = 4.0 * M_PI * 1e-7;
// Accuracy for the secant method
const double eps = 1e-6;

class JAHM2
{
private:
    // Model parameters
    double alpha; // alpha (inked to inter-domain coupling)
    double a; // a (hape factor) A/m
    double Ms; // M_{s} (saturation magnetization) A/m
    double k; // k (inked to hysteresis losses) A/m
    double c; // c (eversibility coefficient)
    double tin;
    double Hbias; // H_{m}
    std::vector<double> Hmvec;

    // Vectors for storing intermediate and final values
    std::vector<double> tH;
    std::vector<double> Hpvec;
    std::vector<size_t> Hspe;
    std::vector<size_t> H0point;

    /**
     * @brief m_an - The vector of a anhysteretic magnetization values
     *               is determined by the modified Langevin formula
     */
    std::vector<double> Man;
    /**
     * @brief dman_dh - The differential equation describes anhysteretic magnetization (M_{an})
     *                  as a function of magnetic field (H)
     */
    std::vector<double> dMandHvec;
    std::vector<double> M1vec;
    std::vector<double> HACvec;
    std::vector<double> Hst;
    std::vector<double> Mst;
    std::vector<double> Bst;
    double Mbias;

    // Auxiliary function for hyperbolic cotangent
    inline double coth(double x);

    // Function to find indices where remainder(tH, 0.5) == 0
    std::vector<size_t> find_special_indices(const std::vector<double>& tH);

    // Secant method for calculating Man and dMandH
    void secant_method(const std::vector<double>& Hpvec, std::vector<double>& Man, std::vector<double>& dMandHvec, const std::vector<size_t>& H0point);

    // Implements the 4th order Runge-Kutta method for integrating the differential equation of the Jiles-Atherton model.
    void rk4_method(const std::vector<double>& Hvec, const std::vector<double>& Manvec, const std::vector<double>& dMandHvec, std::vector<double>& M1vec, double Mbias) const;

public:
    JAHM2(double alpha_ = -0.03,
          double a_ = 12000.0,
          double Ms_ = 389898,//800000, 389898
          double k_ = 4200,//3000, 21000
          double c_ = 0.2,
          double tin_ = 0.01,
          double Hbias_ = 0.0,//40000,
          const std::vector<double>& Hmvec_ = {10000.0, 40000.0, 80000.0})
        :alpha(alpha_)
        , a(a_)
        , Ms(Ms_)
        , k(k_)
        , c(c_)
        , tin(tin_)
        , Hbias(Hbias_)
        , Hmvec(Hmvec_) {}

    void calculate();

    // Getters for result output
    void viewB()
    {
        for (int i = 0; i < Bst.size(); ++i) {
            std::cout << "B: " << Bst[i] << '\n';
        }
    }

    void viewH()
    {
        for (int i = 0; i < Hst.size(); ++i) {
            std::cout << "H: " << Hst[i] << '\n';
        }
    }

    void viewM()
    {
        for (int i = 0; i < Mst.size(); ++i) {
            std::cout << "M: " << Mst[i] << '\n';
        }
    }

    const std::vector<double>& getB() const
    {
        return Bst;
    }
    const std::vector<double>& getH() const
    {
        return Hst;
    }
    const std::vector<double>& getM() const
    {
        return Mst;
    }
    const std::vector<double>& getTime() const
    {
        return tH;
    }

    // Method for save B and H values in to file
    void saveBHToFile(const std::string& filename, size_t step = 10) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error of opening file: " << filename << std::endl;
            return;
        }

        file << "H,B\n";  // Header
        for (size_t i = 0; i < Bst.size(); i += step) {
            file << Hst[i] << "," << Bst[i] << "\n";
        }
        file.close();
    }

};
}
#endif // JAHM2_H
