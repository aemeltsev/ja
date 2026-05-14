#include "jahm3.h"

// Helpful methods
inline double ja::JAHM3::langevin(double x) const
{
    if (std::abs(x) < 1e-4) return x / 3.0;
    return (1.0 / std::tanh(x)) - (1.0 / x);
}

inline double ja::JAHM3::dlangevin(double x) const
{
    if (std::abs(x) < 1e-4) return 1.0 / 3.0;
    double s = std::sinh(x);
    return (1.0 / (s * s)) - (1.0 / (x * x));
}

/**
* @brief Batch calculation of the anhysteresis-free magnetization Man(H) and its derivative.
* Used to visualize the ideal magnetization curve.
*
* @param Hvec Input field strength vector (A/m)
* @param Man Output magnetization vector (A/m)
* @param dMan Output differential susceptibility vector dMan/dHe
*/
void ja::JAHM3::calculate_Man_vec(const std::vector<double> &Hvec,
                                  std::vector<double> &Man,
                                  std::vector<double> &dMan) const
{
    Man.clear();
    dMan.clear();
    Man.reserve(Hvec.size());
    dMan.reserve(Hvec.size());

    for (double H : Hvec)
    {
        // For an ideal curve (without hysteresis) He = H / (1 - alpha * dMan/dHe)
        // But for a simplified display, He = H + alpha * Man is often used
        // Here we implement a direct calculation of the Langevin function:

        double x = H / p.a;
        double Man_val = 0.0;
        double dMan_val = 0.0;

        // Linearization near zero (singularity protection)
        if (std::abs(x) < 1e-4) {
            Man_val = p.Ms * (x / 3.0);
            dMan_val = p.Ms / (3.0 * p.a);
        } else {
            // Eq. (1) from letter
            double coth_x = 1.0 / std::tanh(x);
            Man_val = p.Ms * (coth_x - 1.0 / x);

            // Eq. (12) from letter (Langevin derivative)
            double s = std::sinh(x);
            dMan_val = (p.Ms / p.a) * (1.0 / (s * s) - 1.0 / (x * x));
        }

        Man.push_back(Man_val);
        dMan.push_back(dMan_val);
    }
}

/**
* @brief Basic calculation cycle of the J-A hysteresis model.
*
* The method implements the complete process of magnetization reversal simulation:
* from the initial
* state to a stable (steady-state) limit cycle.
*
* CALCULATION STEPS:
* 1. Time grid: A fixed period (T = 2.0 s) is set and the total simulation time
* is 40 s (20 periods). This is necessary to "warm up" the model so that
* the solution reaches a stable loop (steady-state mode).
*
* 2. Integration (Solver Loop): Loops through all time steps. At each iteration:
* - The current field H(t) is calculated as the sum of the sine wave (Hm) and the bias (Hbias).
* - The magnetization increment algorithm is selected depending on the 'solver' parameter:
* - Euler: Fast calculation of dM using the slope at the current point.
* - RK4: High-precision calculation via 4 stages (rk4_step), minimizing errors on bends.
*
* 3. Period Extraction: Only the last full period (e.g., interval [38 s, 40 s])
* is extracted from the total data set (40,000+ points).
*
* 4. Drift Compensation (Stitching): Eliminates numerical drift. Linearly distributes
* the microscopic error between the start and end points of the period so that
* the loop on the graph is perfectly closed (H[start] == ​​H[end], M[start] == ​​M[end]).
*
* 5. Centering and Calculating B:
* - The average magnetization value (currentMbias) is calculated to eliminate
* constant drift along the M-axis.
* - The final induction vector B = mu0 * (H + M_actual) is calculated.
*
* @param solver Solver type: ja::SolverType::Euler (fast) or ja::SolverType::RK4 (accurate).
*/
void ja::JAHM3::calculate(ja::SolverType type)
{
    /*
     * T = 2.0: We precisely fix the period of the sine wave.
     * stepsPerPeriod: We calculate how many whole steps tin fit into one period (for example, for tin = 0.001, this is exactly 2000).
     * numPeriods = 20: Magnetic materials have "memory." To stop the loop from "creeping" and become stable (stationary), we need to run the model through 15-20 magnetization reversal cycles.
     * totalSteps: The total number of iterations (for example, 40,000).
     * */
    double T = 2.0; // Period of the sine wave (s)
    size_t stepsPerPeriod = static_cast<size_t>(std::round(T / tin));
    const size_t numPeriods = 20;
    size_t totalSteps = stepsPerPeriod * numPeriods; // 20 periods to stabilize the loop

    /*
     * Data buffers. Instead of storing all 40,000 data points,
     * we reserve space only for the last period. This saves RAM.
     * Hm: We take the amplitude of the alternating field from the input data.
     * */
    std::vector<double> Hfull(totalSteps + 1);
    std::vector<double> Mfull(totalSteps + 1, 0.0);
    double Hm = Hmvec.back();

    /*
     * 2. The basic integration loop.
     * t = i * tin: This is the most accurate way to calculate time.
     * It protects against the accumulation of micro-errors like 0.00000000001,
     * which occurs with a simple addition t += tin.
     * */
    for (size_t i = 0; i < totalSteps; ++i) {
        double t = i * tin;
        Hfull[i] = Hm * std::sin(M_PI * t) + Hbias;

        double deltaM = 0.0;

        /*
         * Solver selection. This is the part that makes your program flexible.
         *
         * RK4 (Accuracy) calls the Runge-Kutta method here.
         * It internally invokes J-A physics four times, "feeling" the curve's bend.
         * This ensures we don't miss a sharp turn in the loop near saturation.
         *
         * Euler (Speed) uses Euler's method, which simply takes the slope at the current point
         * and "steps" along a straight line. This is four times faster,
         * but less accurate on sharp bends.
         * */
        if (type == ja::SolverType::RK4) {
            // RK4 branch: high precision (4 derivative measurements within a step)
            deltaM = rk4_step(tin, t, Mfull[i], Hm);
        }
        else {
            // Euler branch: high speed (1 measurement at the beginning of the step)
            double nextH = Hm * std::sin(M_PI * (t + tin)) + Hbias;
            double dH = nextH - Hfull[i];
            double dMdH = get_dMdH_instant(Hfull[i], Mfull[i], dH);
            deltaM = dMdH * dH;
        }

        /*
         * Status update. We add the calculated change in magnetization to the current value.
         * This is the solution to the differential equation.
         * This is how we obtain the next point on the graph.
         * */
        Mfull[i+1] = Mfull[i] + deltaM;
    }

    // We write down the last point of the field H
    Hfull[totalSteps] = Hm * std::sin(M_PI * totalSteps * tin) + Hbias;

    // 3. Selection of the last (stable) period
    size_t startIdx = totalSteps - stepsPerPeriod;
    Hst.clear(); Mst.clear(); Bst.clear();

    /* 4. Drift Compensation (Sewing the ends together for a perfect closure)
     * Calculate the micro-gaps between the beginning and end of the period
     * */
    double driftM = Mfull[totalSteps] - Mfull[startIdx];
    double driftH = Hfull[totalSteps] - Hfull[startIdx];

    for (size_t i = 0; i <= stepsPerPeriod; ++i)
    {
        size_t curr = startIdx + i;
        double alpha = static_cast<double>(i) / stepsPerPeriod;

        // We distribute the error linearly over the entire period
        double H_corr = Hfull[curr] - (driftH * alpha);
        double M_corr = Mfull[curr] - (driftM * alpha);

        Hst.push_back(H_corr);
        Mst.push_back(M_corr);
    }

    /* 5. Centering and final calculation of induction B
     * Find the average value of M to eliminate the constant bias
     * */
    double sumM = std::accumulate(Mst.begin(), Mst.end(), 0.0);
    currentMbias = sumM / Mst.size();

    for (size_t i = 0; i < Hst.size(); ++i)
    {
        /*
         *  Subtract the calculated bias (Mbias) and calculate B = mu0 * (H + M)
         *  IMPORTANT: If you use the formula B = mu0 * (He + M) from your article,
         *  He needs to be recalculated here: double He = Hst[i] + p.alpha * (Mst[i] - currentMbias);
         * */
        Bst.push_back(mu0 * (Hst[i] + (Mst[i] - currentMbias)));
    }
}

/**
* @brief Calculation of dynamic or differential magnetic permeability dM/dH.
*
* The method implements the classical differential form of the Jiles-Atherton (J-A) model. This allows
* to determine the magnetization derivative as a function of the current state of the system (H, M).
*
* Mathematical basis (for more information you can to see
* Jiye Zhao - State Space Representation of Giles–Atherton Hysteresis Model and Application pp. 3-5 pp. 3-5):
* 1. He = H + alpha * M — Effective field (Eq. 2), taking into account the domain interaction.
* 2. Man = Ms * Langevin(He/a) — Anhysteresis-free curve (Eq. 1), the target state of the system.
* 3. dM/dH = (term_rev + term_irr) / (1 - alpha * c * dMan/dHe) — Final equation (Eq. 14).
*
* Physical components:
* - Reversible part (term_rev): Describes the elastic bending of domain walls (parameter c).
* - Irreversible part (term_irr): Describes energy dissipation at structural defects (pinning k).
* - Denominator: Takes into account positive feedback via the interaction parameter alpha.
*
* Implementation features:
* - Condition (Man - M) * delta > 0 guarantees energy stability (physicality).
* - Linearization near zero (He < 1e-4) prevents division by zero (singularity).
*
* @param H Current magnetic field strength (A/m).
* @param M Current magnetization of the material (A/m) — the "memory" of the system.
* @param dH Field increment (used to determine the sign of delta = sgn(dH/dt)).
* @return double The value of the derivative dM/dH at a given point.
*/
double ja::JAHM3::get_dMdH_instant(double H, double M, double dH) const
{
    if (std::abs(dH) < 1e-12) return 0.0;

    /* 1. Direction of field change (delta = 1 for increasing, -1 for decreasing)
     * Field direction (Equation 6)
     * \(\delta = \text{sgn}(\frac{dH}{dt})\)
     * */
    double delta = (dH > 0) ? 1.0 : -1.0;

    /*
     * 2. Calculating the anhysteresis magnetization Man and its derivative
     * Using the Langevin formula: Man = Ms * (coth(He/a) - a/He)
     * Effective field (Equation 2)
     * \(H_e = H + \alpha M\)
     * Input parameter for all other functions, taking into account interdomain interaction.
     * */
    double He = H + p.alpha * M;
    double Man_val = 0.0;
    double dMandH = 0.0;

    if (std::abs(He) < eps) {
        // Linear approximation near zero to avoid division by zero
        Man_val = p.Ms * He / (3.0 * p.a);
        dMandH = p.Ms / (3.0 * p.a);
    } else {
        /*
         * Hysteresis curve (Equation 1 and 5)
         * \(M_{an}=M_{s}\left(\coth \left(\frac{H_{e}}{a}\right)-\frac{H_{e}}{a}\right)\)
         * Here, the "target" magnetization state for the current field is calculated.
         * */
        double x = He / p.a;
        double coth_x = 1.0 / std::tanh(x);
        Man_val = p.Ms * (coth_x - 1.0 / x);
        /* Derivative of the Langevin function: L'(x) = 1/sinh^2(x) - 1/x^2
         * - Derivative of an ideal curve (Equation 12)
         * \(\frac{dM_{an}}{dH_{e}}=\frac{M_{s}}{a}\left(\frac{a^{2}}{H_{e}^{2}}-csch^{2}\left(\frac{H_{e}}{a}\right)\right)\)
         * The code uses \(1/\sinh^2(x)\), which is mathematically identical to the cosecant \(csch^2(x)\) from the article.
         * */

        double sinh_x = std::sinh(x);
        dMandH = (p.Ms / p.a) * (1.0 / (sinh_x * sinh_x) - 1.0 / (x * x));
    }

    /* 3. Calculating the Irreversible Component (Mirr)
     * Physicality condition: dM_irr must always be aligned with (Man - M)
     * Irreversible Component (Equations 3 and 18)
     * \(\frac{dM_{irr}}{dH_e} = \frac{M_{an} - M_{irr}}{\delta k}\) (based on energy balance
     * */
    double phi = (Man_val - M);
    double dM_irr = 0.0;

    /* Energy stability check (prevents unphysical jumps)
     * In the code denominator, delta * p.k - p.alpha * phi is the result of transforming the derivative
     * with respect to \(H\) from the derivative with respect to \(H_{e}\). This is the critical moment,
     * which links the irreversible change directly to the external field.
     * */
    if (phi * delta > 0) {
        dM_irr = phi / (delta * p.k - p.alpha * phi);
    } else {
        dM_irr = 0.0;
    }

    /* 4. Total derivative taking into account the reversible component (p.c.)
     * dM/dH = (1-c)*dMirr/dH + c*dMan/dH
     * Reversible component and total magnetization (Equation 4, 12, 21)
     * \(M = (1-c)M_{irr} + cM_{an}\)
     * This is an implementation of the differential form of equation (14/15) from the article.
     * The coefficient \(c\) determines what part of the change occurs instantaneously (elastically),
     * and what part is lossy (viscous).
     * */
    double dMdH = (1.0 - p.c) * dM_irr + p.c * dMandH;

    return dMdH;
}

/**
* @brief Numerical integration of the magnetization equation
* using the fourth-order Runge-Kutta method (RK4).
*
* The method calculates the magnetization increment (delta M) in one time step (h).
* Unlike the Euler method, RK4 performs four trial calculations of the derivative
* within an interval, ensuring O(h^4) accuracy and high stability in
* steep sections of the hysteresis loop (near the coercivity Hc).
*
* Algorithm (according to Flowchart):
* 1. k1 = f(t, M) * dH
* 2. k2 = f(t + h/2, M + k1/2) * dH
* 3. k3 = f(t + h/2, M + k2/2) * dH
* 4. k4 = f(t + h, M + k3) * dH
*
* The final increment is calculated as a weighted average:
* deltaM = (k1 + 2*k2 + 2*k3 + k4) / 6
*
* @param h Time integration step (tin).
* @param t Current simulation time (s).
* @param M Current magnetization value (A/m) — "memory" state.
* @param Hm External magnetic field amplitude (A/m).
* @return double Calculated magnetization increment dM per step h.
*/
double ja::JAHM3::rk4_step(double h, double t, double M, double Hm) const
{

    // h this is tin (time step)

    /*
     * Allows us to calculate the exact value of the field strength \(H\)
     * at any point in time within a step \(h\) (at the beginning, middle, or end).
     * To determine the slope at the middle of a step (\(t + h/2\)),
     * we need to know what the external field will be at that precise micro-moment.
     * */
    auto getH = [&](double time) {
        return Hm * std::sin(M_PI * time) + Hbias;
    };

    /*
     * Fixing the dH_full direction. This is a critical detail for the J-A model.
     * The \(\delta\) (delta) parameter in the model is defined as
     * the sign of the field change (\(dH\)).
     * If we calculate \(\delta\) separately for each stage (k1-k4),
     * then at the peaks of the sinusoid (where the field changes direction),
     * RK4 may begin to "flip" between the hysteresis branches,
     * leading to numerical instability. We fix the direction for the entire step.
     * */
    double dH_full = getH(t + h) - getH(t);

    /*
     * The RK4 method approximates the curve not with a single straight line (like Euler),
     * but with a parabola.
     *
     * k1 (Start): Calculate the magnetization increment based on the current state.
     * This is the "base" trajectory.
     *
     * k2 (Prediction 1): Step into the middle of the interval, using slope k1.
     * We check: "What will the magnetic susceptibility be there
     * if we slightly magnetize the material?"
     *
     * k3 (Prediction 2): Step into the middle again, but now use slope k2.
     * This is a refinement of the trajectory.
     *
     * k4 (Final): Step to the end of the step, using slope k3.
     *
     * At each stage, we call get_dMdH_instant, which returns the "slope" of the curve,
     * and multiply it by the change in field (\(dH_{full}\)),
     * obtaining the magnetization increment \(dM\).
     * */
    double k1 = get_dMdH_instant(getH(t), M, dH_full) * dH_full;

    double k2 = get_dMdH_instant(getH(t + h/2.0), M + k1/2.0, dH_full) * dH_full;

    double k3 = get_dMdH_instant(getH(t + h/2.0), M + k2/2.0, dH_full) * dH_full;

    double k4 = get_dMdH_instant(getH(t + h), M + k3, dH_full) * dH_full;

    /*
     * Weighted average. This is the mathematical core of RK4.
     * We sum all the weighted predictions. The extremes (k1 and k4) are given a weight of 1.
     * The midpoints (k2 and k3) are given a weight of 2,
     * as they more accurately describe the function's behavior within the interval.
     * Dividing by 6 returns the weighted average increment.
     * */
    return (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

/**
* @brief Dynamically updates the parameters of the model depending on temperature.
*
* The method uses linear interpolation between two sets of reference data
* (typically 25°C and 100°C) to correct the material properties. This allows simulating the
* effect of "sag" in magnetic properties as the core heats up.
*
* @param p25 Structure of parameters measured at room temperature (25°C).
* @param p100 Structure of parameters measured at operating temperature (100°C).
* @param T_celsius Current calculated material temperature.
*/
void ja::JAHM3::updateTemperature(ja::HysteresisParams p25, ja::HysteresisParams p100, double T_celsius)
{
    // 1. Calculate the interpolation factor (0.0 for 25°C, 1.0 for 100°C)
    double factor = (T_celsius - 25.0) / (100.0 - 25.0);

    // 2. Protect against extrapolation beyond reasonable limits.
    // A limit of 1.5 corresponds to approximately 137.5°C, which is close to the Curie point for many ferrites.
    factor = std::clamp(factor, 0.0, 1.5);

    // 3. Interpolate the saturation magnetization (Ms) and the shape parameter (a).
    // When heated, the Ms of ferrites drops significantly, which lowers the saturation threshold Bs.
    p.Ms = p25.Ms + factor * (p100.Ms - p25.Ms);
    p.a = p25.a + factor * (p100.a - p25.a);

    // 4. Interpolation of the pinning coefficient (k).
    // Typically, k decreases with heating, which makes the hysteresis loop "narrow"
    // and reduces static losses per cycle.
    p.k = p25.k + factor * (p100.k - p25.k);

    // 5. Correction of the reversibility (c) and interaction (alpha) parameters.
    // These quantities are less dependent on T, but taking them into account improves the accuracy
    // of the approximation of the "knee" of the magnetization curve.
    p.c = p25.c + factor * (p100.c - p25.c);
    p.alpha = p25.alpha + factor * (p100.alpha - p25.alpha);
}

/**
* @brief Analyze the calculated loop to extract
* the material's physical characteristics.
* The method finds the coercivity (Hc), remanent flux density (Br),
* and saturation flux density (Bs).
*
* @return MaterialMetrics A structure with the calculated values ​​of Hc, Br, and Bs.
*/
ja::MaterialMetrics ja::JAHM3::calculateMetrics() const
{
    MaterialMetrics m = {0.0, 0.0, 0.0};

    // Protection against empty vector (if the calculation has not been performed yet)
    if (Hst.empty() || Bst.empty()) return m;

    /* 1. Calculating B_s (Saturation Induction)
     * Bs is defined as the induction amplitude: (max - min) / 2.
     * This allows for correct calculation even in the presence of vertical displacement.
     * */
    auto [minB, maxB] = std::minmax_element(Bst.begin(), Bst.end());
    m.Bs = (*maxB - *minB) / 2.0;

    // Iterate over all points of the stationary period to find intersections of the axes
    for (size_t i = 1; i < Hst.size(); ++i)
    {
        /* 2. Finding Hc (Coercive Force)
         * Find the point where the curve intersects the H-axis
         * (induction B becomes equal to 0).
         * The condition B[i-1] * B[i] <= 0 means that at this step
         * the induction has changed sign.
         * */
        if (Bst[i-1] * Bst[i] <= 0)
        {
            /* Linear interpolation to find the exact zero between grid nodes.
             * t is the weighting factor (how close the zero is to point i-1).
             * 1e-18 is added to avoid division by zero.
             * */
            double t = std::abs(Bst[i-1]) / (std::abs(Bst[i-1]) + std::abs(Bst[i]) + 1e-18);
            double H_zero = Hst[i-1] + t * (Hst[i] - Hst[i-1]);

            // We take the maximum, since there are two intersections (positive and negative).
            m.Hc = std::max(m.Hc, std::abs(H_zero));
        }

        /* 3. Finding Br (Residual Induction)
         * Find the point where the curve intersects the B-axis
         * (the H-field strength becomes 0).
         * */
        if (Hst[i-1] * Hst[i] <= 0)
        {
            // Linear interpolation to find the exact value of B when H = 0.
            double t = std::abs(Hst[i-1]) / (std::abs(Hst[i-1]) + std::abs(Hst[i]) + 1e-18);
            double B_zero = Bst[i-1] + t * (Bst[i] - Bst[i-1]);

            m.Br = std::max(m.Br, std::abs(B_zero));
        }
    }
    return m;
}

/**
* @brief Calculation of static magnetic energy losses per cycle (J/m³).
*
* The method calculates the area of ​​the hysteresis loop,
* which physically corresponds to
* the energy converted to heat due to magnetic domain friction (pinning)
* during one complete magnetization reversal cycle.
* The integral $\int H \, dB$ sums all the horizontal stripes within the loop;
* the sum of these stripes gives the exact area of ​​the figure.
* That is, each time a domain wall "jumps" over a defect in
* the crystal (pinning), part of the energy of
* the magnetizing field $H$ is irreversibly converted into heat.
* The wider the loop (the larger the $k$ parameter of the model),
* the greater the area and the higher the losses.
*
* @return double Energy loss per unit volume (J/m³).
*/
double ja::JAHM3::calculateLosses() const
{
    // Protection: at least 2 points are required to form the area.
    if(Hst.size() < 2) return 0.0;

    double area = 0.0;

    // Numerical integration over a closed contour using the trapezoidal method.
    // Mathematically, this is an approximation of the integral: W = \oint H dB
    for(size_t i = 1; i < Hst.size(); ++i)
    {
        /*
         * (Hst[i] + Hst[i-1]) * 0.5 — average field strength per step.
         * (Bst[i] - Bst[i-1]) — magnetic induction increment (segment height).
         * The output gives the area of ​​a narrow horizontal trapezoid.
         * */
        area += (Hst[i] + Hst[i-1]) * (Bst[i] - Bst[i-1]) * 0.5;
    }

    /*
     * The integral may be biased depending on the direction of traversal.
     * (clockwise or counterclockwise). In energy physics, losses are
     * dimensions are always positive, hence the pregnant module.
     * */
    return std::abs(area);
}

/**
* @brief Calculates the total specific power loss in the material (W/m³).
*
* The method combines static hysteresis losses and dynamic
* eddy losses according to the Bertotti loss separation.
*
* @param frequency Operating frequency of magnetization reversal (Hz).
* @return double Total power loss density per unit volume (W/m³).
*/
double ja::JAHM3::calculateTotalPowerLoss(double frequency) const
{
    /* 1. Hysteresis component (Ph)
     * Calculated as the static loop energy (J/m³) multiplied by the frequency (1/s).
     * Physically, these are losses due to domain wall friction.
     * Increases linearly with frequency.
     * */
    double Ph = calculateLosses() * frequency;

    /* 2. Calculating the induction amplitude
     * For dynamic losses, the maximum induction (Bs) in a cycle is critical.
    */
    MaterialMetrics m = calculateMetrics();

    /* 3. Eddy current losses (Pe)
     * According to classical electrodynamics, eddy currents are induced by an alternating
     * magnetic field. These losses are proportional to the square of the frequency and the square of the induction.
     * k_eddy is a coefficient depending on the specific electrical resistance
     * of the material and its geometry (sheet thickness or ferrite grain size).
     * */
    double k_eddy = 0.5e-4; // Empirical constant for power ferrites (MnZn)

    // Eq: Pe = k_eddy * f² * B²
    double Pe = k_eddy * std::pow(frequency * m.Bs, 2);

    // 4. Summation
    // The final result is the thermal power that the core must dissipate.
    return Ph + Pe;
}

/**
* @brief Generates a time vector for a stationary period.
*
* Creates a time scale from 0 to T, where T is the duration of one period.
* This is necessary for the correct display of the B(t) and H(t) waveforms in the GUI,
* so that the graphs start at time zero.
*
* @return std::vector<double> Vector of timestamps with a step of tin.
*/
std::vector<double> ja::JAHM3::getTime() const
{
    // The size of the time vector must strictly match the size of the Hst and Bst vectors
    std::vector<double> t;
    if (Hst.empty()) return t;

    t.reserve(Hst.size());
    for (size_t i = 0; i < Hst.size(); ++i) {
        // Calculate the time of the i-th point relative to the beginning of the period
        // t = i * tin ensures the absence of accumulated summation error
        t.push_back(static_cast<double>(i) * tin);
    }
    return t;
}

bool ja::JAHM3::saveBHToFile(QString name)
{
    QFile file(name);

    // Trying to open file. If we cannot open file return false.
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        // qWarning() output system reason of error to debug console
        qWarning() << "Cannot open file for writing:" << file.errorString();
        return false;
    }

    QTextStream stream(&file);
    // If size of output vectors is equal
    if (Hst.size() == Bst.size()) {
        for(size_t i = 0; i < Hst.size(); ++i) {
            stream << Hst[i] << "," << Bst[i] << "\n";
        }
    }

    // Check if the disk is full and if the data was written successfully
    if (stream.status() != QTextStream::Ok) {
        qWarning() << "Write data to stream error.";
        return false;
    }

    file.close();
    return true; // Everything went well.
}


