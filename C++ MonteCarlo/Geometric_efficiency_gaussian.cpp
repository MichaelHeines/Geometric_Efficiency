/*
IMPORTANT NOTE: This program determines the geometric efficiency for a (possibly annular) circular detector - gaussian source system using Montecarlo
integration. The compiler needs to be at least c++17 due to the use of the cylindrical bessel functions of the std library added in that version.

It is critical that the minimum distance z_low is never set to 0, as the integral will not converge causing an unusable point.
Furthermore, the distances are all given in units of the detector radius. A file can be written out with columns showing the distance
z, the geometric efficiency using a gaussian source distribution, and the geometric efficiency of the point source distribution (for
comparing purposes). In the annular detector option, the rescaled distance is with respect to the outer detector.
*/

#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <string>

// Creates a linspace starting from start_in, going to end_in with num_in total steps
std::vector<double> linspace(double start_in, double end_in, int num_in) {
    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i) {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end);
    return linspaced;
}


// Calculates the point source approximation for a vector of distances z with lenght nr_points
std::vector<double> point_source(std::vector<double> z, int nr_points) {
    std::vector<double> approx(nr_points);
    double temp;

    for (int i = 0; i < nr_points; i++) {
        temp = .5 - (.5 * z[i]) / (sqrt(1 + pow(z[i], 2)));
        approx[i] = temp;
    }

    return approx;
}


// Provides the integrand for the mc integration of the circular-circular system
double integrand(double k, double r) {
    double f = std::cyl_bessel_j(1, k) * std::cyl_bessel_j(0, k * r) * r;
    return f;
}


// Calculates the integral over r using mc integration
double mc_int_r(int n_r, double k, double sigma) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> gaus(0, sigma);
    double I_r = 0;
    double r;

    for (int i = 0; i < n_r; i++) {
        r = gaus(gen);

        if (r < 0) {
            r = -r;
        }

        I_r += (1 / double(n_r)) * integrand(k, r);
    }

    return I_r;
}


// Calculates the integral over k (following r) using mc integration with n_k samples (that is increased for small distances)
double mc_int_k(int n_k0, int n_r, double z, double sigma) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> exp(z);
    double I_k = 0;
    double pi = 3.14159265358979323846;
    double n_k;
    double k;
    double prefactor = sqrt(pi / 2) / (2 * z * sigma);

    if (z >= 0.8) {
        n_k = n_k0;
    }
    else if (z >= 0.5) {
        n_k = 10 * n_k0;
    }
    else if (z >= 0.3) {
        n_k = 30 * n_k0;
    }
    else {
        n_k = 100 * n_k0;
    }

    for (int i = 0; i < n_k; i++) {
        k = exp(gen);
        I_k += prefactor * (1 / n_k) * mc_int_r(n_r, k, sigma);
    }

    return I_k;
}


// Writes a file for the geometric efficiency as a function of distance of the circular-cirular system
void write_geo_file(double sigma, double z_low, double z_high, int n_r, int n_k, int nr_points, const std::string file_name) {
    double result;
    std::ofstream myFile(file_name);
    std::vector<double> e_geo(nr_points);
    std::vector<double> z = linspace(z_low, z_high, nr_points);
    std::vector<double> e_ps = point_source(z, nr_points);
    myFile << "z/rd ; gaussian ; point source \n";

    for (int i = 0; i < nr_points; i++) {
        result = mc_int_k(n_k, n_r, z[i], sigma);
        e_geo[i] = result;
        std::cout << (z[i] - z_low) / (z_high - z_low) << "\t" << result << std::endl;
        myFile << z[i] << "; " << e_geo[i] << "; " << e_ps[i] << "\n";
    }

    myFile.close();
}


// Writes a file for the geometric efficiency as a function of distance of the angular circular - gaussian system
void write_annular_file(double sigma, double det_ratio, double z_low, double z_high, int n_r, int n_k, int nr_points, const std::string file_name) {
    double result_inner;
    double result_outer;
    std::ofstream myFile(file_name);
    std::vector<double> e_geo(nr_points);
    std::vector<double> z_outer = linspace(z_low, z_high, nr_points);                                   // In units of the outer radius
    std::vector<double> z_inner = linspace(z_low * det_ratio, z_high * det_ratio, nr_points);               // In units of the outer radius
    std::vector<double> e_ps_outer = point_source(z_outer, nr_points);
    std::vector<double> e_ps_inner = point_source(z_inner, nr_points);
    std::vector<double> e_ps(nr_points);

    for (int i = 0; i < nr_points; i++) {
        e_ps[i] = e_ps_outer[i] - e_ps_inner[i];
    }

    myFile << "z/rd ; gaussian annular ; point source \n";

    for (int i = 0; i < nr_points; i++) {
        result_outer = mc_int_k(n_k, n_r, z_outer[i], sigma);
        result_inner = mc_int_k(n_k, n_r, z_inner[i], sigma);
        e_geo[i] = result_outer - result_inner;
        std::cout << (z_outer[i] - z_low) / (z_high - z_low) << "\t" << e_geo[i] << std::endl;
        myFile << z_outer[i] << "; " << e_geo[i] << "; " << e_ps[i] << "\n";
    }

    myFile.close();
}


// Caculates the Mean and Gaussian error at a specific distance
double error_estimate(double sigma, double z, int n_r, int n_k, int n_std) {
    double result;
    double mean = 0;
    double mean_sq = 0;
    double progress = 0;
    double std;
    std::vector<double> e_geo(n_std);

    for (int i = 0; i < n_std; i++) {
        result = mc_int_k(n_k, n_r, z, sigma);
        e_geo[i] = result;
        progress += 1 / double(n_std);
        std::cout << progress << std::endl;
    }

    for (int i = 0; i < n_std; i++) {
        mean += e_geo[i] / double(n_std);
        mean_sq += e_geo[i] * e_geo[i] / double(n_std);
    }

    std = sqrt(mean_sq - mean * mean);
    std::cout << mean << "\t" << 100 * std / mean << std::endl;
    return std;
}


int main() {
    double sigma = 1;
    double det_ratio = 1.5;
    double z = 1;
    double z_low = 1;
    double z_high = 3;
    int n_r = pow(10, 2);
    int n_k = pow(10, 3);
    int nr_points = 100;
    int n_std = 100;
    const std::string file_name_geo = "geo_test_cpp.csv";
    const std::string file_name_annular = "annular_test_cpp.csv";

    write_geo_file(sigma, z_low, z_high, n_r, n_k, nr_points, file_name_geo);
    write_annular_file(sigma, det_ratio, z_low, z_high, n_r, n_k, nr_points, file_name_annular);
    error_estimate(sigma, z, n_r, n_k, n_std);
    return 0;
}


