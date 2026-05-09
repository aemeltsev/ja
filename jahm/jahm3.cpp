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

void ja::JAHM3::calculate_Man_vec(const std::vector<double> &Hvec, std::vector<double> &Man, std::vector<double> &dMan) const
{
    // TODO

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
