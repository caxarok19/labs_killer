#include "funcs.h"
#include<fstream>


Eigen::Vector<double, 1> lin(const Eigen::Vector<double, 1> state, const double t) {
    Eigen::Vector<double, 1> res;
    res = Eigen::Vector < double, 1>::Zero();
    res[0] = t * t * t;
    return(res);
}
double ans1(double t) { return (t * t * t * t / 4); }
Eigen::Vector<double, 2> osc(const Eigen::Vector<double, 2> state, const double t) {
    Eigen::Vector<double, 2> res = { state[1], -state[0] };
    return(res);
}
double ans2(double t) { return (std::sin(t)); }


Eigen::Vector<double, 4> func_orb(const Eigen::Vector<double, 4> state, const double t) {
    const double mu = 0.012277471;
    const double eta = 1 - mu;
    double A = std::pow((state[0] + mu) * (state[0] + mu) + state[2] * state[2], 3.0 / 2.0);
    double B = std::pow((state[0] - eta) * (state[0] - eta) + state[2] * state[2], 3.0 / 2.0);
    Eigen::Vector<double, 4> res = { state[1],
    state[0] + 2 * state[3] - eta * (state[0] + mu) / A - mu * (state[0] - eta) / B,
    state[3],
    state[2] - 2 * state[1] - eta * state[2] / A - mu * state[2] / B };
    return (res);
}

 

Eigen::Vector<double, 6> orbit(const Eigen::Vector<double, 6> state, const double t) {
    Eigen::Vector<double, 6> res;
    const double mu = 398600.441589;
    const double e = 3.0 * mu * 6371.0 * 6371.0 * 1.08263 / 2.0;
    const double r = std::sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
    const double A = -mu / std::pow(r, 3);
    const double B = 5.0 * e / std::pow(r, 7);
    const double C = e / std::pow(r, 5);

    res = { state[3], state[4], state[5],
            A * state[0] + B * state[0] * state[2] * state[2] + C * state[0],
            A * state[1] + B * state[1] * state[2] * state[2] + C * state[1],
            A * state[2] + B * state[2] * state[2] * state[2] - C * state[2] };

    return res;
}

СauchyProblem<decltype(lin), 1> linear;
СauchyProblem<decltype(osc), 2> oscillator;
СauchyProblem<decltype(orbit), 6> Earth;
СauchyProblem<decltype(func_orb), 4>Orbit;

 

int main() {

    std::ofstream data3("data_III.txt");
    data3.precision(16);

    auto res = integrate<decltype(func_orb), BDF4, СauchyProblem<decltype(func_orb), 4>, RK4Table>(
        func_orb, { Eigen::Vector<double, 4>(0.994, 0.0, 0.0, -2.00158510637908252240537862224), 0 }, 100, {0.001, 1e-24,10}, Orbit);

    for (unsigned int i = 0; i < res.size(); i++) {
        data3 << res[i].state[0] <<  " " << res[i].state[2] << std::endl;
    }

    data3.close();

    return 0;
};
