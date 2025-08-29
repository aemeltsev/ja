# ja :curly_loop:

 The Jiles-Atherton (J-A) model is a commonly used physics-based model in
 describing the hysteresis characteristics of ferromagnetic materials.

 * the secant method used for solving anhysteretic magnetization with some optimized for faster convergence
 * for estimate magnetization was employed the Fourth Order Runge-Kutta method

# **General scheme of the algorithm:**

Create a time vector `tH` to simulate a periodic signal. The `tH` vector is filled with time values ​​from 0 to 40 with a step of `tin`. If `tin = 0.004`, then the vector will contain 10000 points.
Create a magnetic field strength vector `Hpvec`, changing according to a sinusoidal law.
Each element of `Hpvec` is calculated as `Hbias * sin(π * t)`, where `Hbias` is the bias amplitude.
The `find_special_indices` function finds indices where the time signal is a multiple of 0.5 (`remainder(tH[i], 0.5) == 0`). These indices are used to determine special points where the magnetization is zeroed.
Initialize the anhysteresis-free magnetization vectors `Man` and its derivative `dMandHvec`.

Calculate the anhysteresis-free magnetization `Man` and its derivative `dMandHvec` using the secant method.

```cpp
secant_method(Hpvec, Man, dMandHvec, H0point);
```
Check for belonging to special points and if the current index belongs to `H0point`, then `Man[iH] = 0` and the derivative is calculated analytically.
Iteration process for the remaining points. `prop0` and `prop1` are initial approximations for the secant method.
The secant method iterations are performed until the specified accuracy or maximum number of iterations is reached.
```cpp
while (std::abs(fMan_sec1) >= eps * std::abs(Man_sec1) && isec < max_iter) {
double Man_sec2 = Man_sec1 - fMan_sec1 * (Man_sec1 - Man_sec0) / (fMan_sec1 - fMan_sec0);
// Update values ​​for the next iteration
}
```
The derivative is calculated using the analytical formula.
```cpp
double T = (Hp + alpha * Man_sec1) / a;
double Tsin2 = std::pow(std::sinh(T), 2);
double T2 = std::pow(T, 2);
double denominator = (a * Tsin2 * T2 / Ms) - (alpha * (Tsin2 - T2));
dMandHvec[iH] = (Tsin2 - T2) / denominator;
```

Prepare `M1vec` vector to store total magnetization. Calculate the differential equation of total magnetization using the 4th order Runge-Kutta method.

```cpp
rk4_method(Hpvec, Man, dMandHvec, M1vec, 0.0);
```

Determine the directions of change of magnetic field.
```cpp
double delta = (Hvec[n] > Hvec[n-1]) ? 1.0 : -1.0;
```
Calculate the coefficients `k1`, `k2`, `k3`, `k4`.
```cpp
double k1 = (delta * k * c * dMandH + phi1) / (delta * k - alpha * phi1);
double M1_k2 = M1vec[n-1] + h * k1 / 2.0;
double k2 = (delta * k * c * dMandH + phi1) / (delta * k - alpha * phi1);
// Similarly for k3 and k4
```

Update the magnetization values.
```cpp
M1vec[n] = M1vec[n-1] + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
```

Find the maximum magnetization value to adjust subsequent calculations.

```cpp
Mbias = *std::max_element(M1vec.begin(), M1vec.end());
```

We iteratively calculate for each value of the magnetic field amplitude `Hm` from the vector `Hmvec`.

```cpp
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
```

Recalculating the anhysteresis-free magnetization and the total magnetization.

```cpp
std::vector<double> Man_new(Hpvec_new.size(), 0.0);
std::vector<double> dMandHvec_new(Hpvec_new.size(), 0.0);
secant_method(Hpvec_new, Man_new, dMandHvec_new, H0point);

std::vector<double> M1vec_new(Hpvec_new.size(), Mbias);
rk4_method(Hpvec_new, Man_new, dMandHvec_new, M1vec_new, Mbias);
```

We obtain the magnetization and magnetic field values ​​for the stationary section (the last 10% of the signal).

```cpp
size_t numst_start = static_cast<size_t>(std::floor(Hpvec_new.size() * 0.9));
std::vector<double> tst(tH.begin() + numst_start, tH.end());
Hst.assign(Hpvec_new.begin() + numst_start, Hpvec_new.end());
Mst.assign(M1vec_new.begin() + numst_start, M1vec_new.end());
for (auto& val : Mst) {
val -= Mbias;
}
```

Calculate magnetic induction on a stationary section.

```cpp
Bst.resize(Hst.size());
for (size_t i = 0; i < Hst.size(); ++i) {
Bst[i] = mu0 * (Hst[i] + Mst[i]);
}
```

# Literature
* Gustav Mörée - Review of Hysteresis Models for Magnetic Materials
* Guangming Xue - Numerical Solving Method for Jiles-Atherton Model and Influence Analysis of the Initial Magnetic Field on Hysteresis
* Nowicki M.-Modeling the Hysteresis Loop of Ultra-High Permeability Amorphous Alloy for Space Applications
* Кружаев А.В.-Компьютерное моделирование и экспериментальное исследование переходных процессов в однофазном трансформаторе напряжения
* Chowdhury J.-Real-time Physical Modeling for Analog Tape Machines
* Upadhaya B.-A constraint-based optimization technique for estimating physical parameters of Jiles–Atherton hysteresis model
* Jastrzebski R.-Comparison of macroscopic descriptions of magnetization curves
* Knypinski L.-Application of a PSO algorithm for identification of the parameters of Jiles-Atherton hysteresis model
