#include "jahm3.h"

// Вспомогательные методы
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

