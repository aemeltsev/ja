#ifndef JAHM3_H
#define JAHM3_H
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>
#include <QFile>
#include <QTextStream>
#include <QDebug>

namespace ja {

// Mathematical and physical constants
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps = 1e-4;

enum class SolverType
{
    Euler,
    RK4
};

// Common structure of parameters for all components of model
struct HysteresisParams
{
    double a;      // Form parameter (A/m)
    double k;      // Hysteresis parameter (A/m)
    double c;      // Reversibility parameter (0...1)
    double Ms;     // Saturation magnetization (A/m)
    double alpha;  // Interaction parameter
    double Kan;    // Anisotropy constant
    double psi;    // Anisotropy angle
};

struct MaterialMetrics
{
    double Hc; // Coercive force (A/m)
    double Br; // Residual induction (T)
    double Bs; // Saturation induction (T)
};

struct CoilGeometry
{
    double N;  // Turns
    double S;  // Square (m^2)
    double Le; // Line length (m)
};

class JAHM3
{
private:
    HysteresisParams p;
    /* Time integration step (discretization).
     * The smaller it is, the higher the accuracy
     * and the less drift.
     */
    double tin;
    /* Сonstant component of
     * the field strength (loop displacement).
     */
    double Hbias;
    /*
     * Vector of amplitudes of the alternating field \(H_{m}\).
     * Allows one to calculate a family of loops.
     */
    std::vector<double> Hmvec;

    /*
     * The resulting vectors of one stationary period
     * (field, magnetization, induction, and time).
     * The suffix st means stationary.
     * */
    std::vector<double> Hst, Mst, Bst, tH;
    // Calculated average magnetization value for loop centering.
    double currentMbias;

    /*
     * The Langevin function \(L(x) = \coth(x) - 1/x\).
     * Describes the ideal (hysteresis-free) state.
     * */
    inline double langevin(double x) const;

    /*
     * The derivative of the Langevin function.
     * Necessary for calculating differential permeability.
     * */
    inline double dlangevin(double x) const;

     /*
      * Batch calculation of the anhysteresis-free
      * magnetization vector for a given field array \(H\).
      * */
    void calculate_Man_vec(const std::vector<double>& Hvec,
                           std::vector<double>& Man,
                           std::vector<double>& dMan) const;

public:
    /* AHM3 Constructor - Initializes the model parameters.
     * By default, it sets the step size to 0.001
     * and the test field to 800 A/m. */
    JAHM3(const HysteresisParams& params,
          double tin_ = 0.001,
          double Hbias_ = 0.0,
          const std::vector<double>& Hmvec_ = {800.0})
    : p(params)
    , tin(tin_)
    , Hbias(Hbias_)
    , Hmvec(Hmvec_)
    , currentMbias(0.0) {}

    /* The main method. Starts a 20-period integration cycle (warm-up),
     * performs "stitching" of the period ends (Drift Compensation),
     * and fills the Hst, Mst, and Bst vectors. */
    void calculate(SolverType type);

    /* Returns the value of \(dM/dH\) at a specific point.
     * Critical for coupling to a coil (JACoil),
     * as it determines the instantaneous inductance. */
    double get_dMdH_instant(double H, double M, double dH) const;

    /*
     * This method calculates the magnetization increment for one step.
    */
    double rk4_step(double h, double t, double M, double Hm) const;

    /* Dynamically recalculates the parameters p.Ms and p.k
     * depending on the heating (linear interpolation
     * between 25°C and 100°C). */
    void updateTemperature(HysteresisParams p25, HysteresisParams p100, double T_celsius);

    /* Analyzes the finished loop
     * and finds the coercivity (\(H_{c}\)),
     * residual induction (\(B_{r}\))
     * and saturation (\(B_{s}\)). */
    MaterialMetrics calculateMetrics() const;

    /* Integrates the loop area (integral \(\oint H \, dB\)),
     * returning the energy lost per cycle in J/m^{3}. */
    double calculateLosses() const;

    /* Converts cyclic energy into power (W/m^{3})
     * and can take into account additional eddy losses. */
    double calculateTotalPowerLoss(double frequency) const;

    /* Return references to data vectors for plotting graphs */
    const std::vector<double>& getB() const { return Bst; }
    const std::vector<double>& getH() const { return Hst; }
    const std::vector<double>& getM() const { return Mst; }

    /* A function for loading parameters for 3C94 ferrite.
     * Allows you to quickly load reference values ​​for 25°C
     * or 100°C. */
    static HysteresisParams get3C94Params(bool hot = false)
    {
        if (hot){ //100°C
            return HysteresisParams{22.0, 18.0, 0.22, 255000, 1.5e-4, 0.0, 0.0};
        } else { // 25°C
            return HysteresisParams{35.5, 28.2, 0.18, 374000, 1.1e-4, 0.0, 0.0};
        }
    }

    /* Generates a time scale for the X-axis of the graph
     * (from 0 to \(T\) in increments of tin). */
    std::vector<double> getTime() const;

    bool saveBHToFile(QString name);
    double getAlpha() const;
};

class JACoil
{
private:
    std::unique_ptr<JAHM3> m_model; // Link to precision physical model of JA
    CoilGeometry m_geo; // Geometry parameters of real device (N, S, Le)

    // State of the core "memory" (State Variables)
    double m_lastH = 0.0;              // Previous field value (A/m)
    double m_lastM = 0.0;              // Previous magnetization value (A/m)
    double m_currentI = 0.0;           // Current in the winding (A)

    // Thermal state
    double m_temperature;        // Current core termal value (°C)

public:
    /**
     * @brief Constructor of an inductor with a nonlinear core.
     * Starting state: a completely demagnetized, cold core.
     * Transfer ownership via std::move
     */
    JACoil(std::unique_ptr<JAHM3> m, const CoilGeometry& g, double initial_tmp)
        :m_model(std::move(m))
        ,m_geo(g)
        ,m_temperature(initial_tmp)
    {}

    // Remove the copying since unique_ptr cannot be copied.
    JACoil(const JACoil&) = delete;
    JACoil& operator=(const JACoil&) = delete;

    JACoil(JACoil&&) = default;
    JACoil& operator=(JACoil&&) = default;

    // 1. Transfer Current (diagram) -> H-Field (material)
    double currentToH(double I) const;

    /**
     * @brief Main simulation step. Updates the ferrite's "memory" state.
     * Called by the simulator at each time step dt when a new current is applied.
     */
    void updateState(double newI, double deltaTime);

    /**
     * @brief Calculates instantaneous differential inductance L(I).
     * @param dI Current increment (necessary to determine delta direction)
     */
    double getDynamicInductance(double dI) const;

    /**
     * @brief Calculates the instantaneous voltage across the coil when the current changes.
     * Used in circuit simulators (LTspice / MATLAB).
     */
    double getVoltage(double dI, double deltaTime) const;

    /**
     * @brief Check saturation margin (in percent) taking temperature into account.
     */
    double getSaturationMargin() const;

    // Set temperature value from external thermal simulation
    void setTemperature(double T_celsius);

    // Current inductive value in magnetic core
    double getB() const;

    // Reset (demagnetise)
    void reset();

    double getRelativePermeability() const;
};

}
#endif // JAHM3_H
