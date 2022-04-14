#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <fstream>
#define pi 3.14159265358979323846

class position {
    public:
        std::vector<double> x, y;

        position(std::vector<double> x_in, std::vector<double> y_in){
            x = x_in;
            y = y_in;
        }

        void add_vec(std::vector<double> x_diff, std::vector<double> y_diff){
            int size_1 = x.size();
            int size_2 = x_diff.size();

            if (size_1 != size_2){
                std::cerr << "ERROR: vector addtion with different size vectors" << std::endl;
                exit(0);
            }

            for (int i = 0; i < size_1; i++){
                x[i] += x_diff[i];
                y[i] += y_diff[i];
            }
        }

        std::vector<double> calculate_rsq(){
            int size = x.size();
            std::vector<double> r_sq(size);

            for (int i = 0; i < size; i++){
                r_sq[i] = x[i]*x[i] + y[i]*y[i];
            }
            return r_sq;
        }

        // On empty position
        void generate_isotropic(double z, int seed){
            std::default_random_engine generator{seed};
            std::uniform_real_distribution<double> phi_distr(0, 2*pi);
            std::uniform_real_distribution<double> theta_create_distr(0, 1);
            double phi, theta;

            for (int i = 0; i < x.size(); i++){
                phi = phi_distr(generator);
                theta = acos(1 - 2 * theta_create_distr(generator));
                x[i] = z * tan(theta) * cos(phi);
                y[i] = z * tan(theta) * sin(phi);
            }
        }

        // On empty position
        void generate_circular_distr(double r_s, int seed){
            std::default_random_engine generator{seed};
            std::uniform_real_distribution<double> phi_distr(0, 2*pi);
            std::uniform_real_distribution<double> r_create_distr(0, 1);
            double phi, r;

            for (int i = 0; i < x.size(); i++){
                phi = phi_distr(generator);
                r = r_s * sqrt(r_create_distr(generator));
                x[i] = r * cos(phi);
                y[i] = r * sin(phi);
            }
        }

        // On empty position
        void generate_gaussian_distr(double sigma, int seed){
            std::default_random_engine generator{seed};
            std::uniform_real_distribution<double> phi_distr(0, 2*pi);
            std::normal_distribution<double> r_distr(0, sigma);
            double phi, r;

            for (int i = 0; i < x.size(); i++){
                phi = phi_distr(generator);
                r = r_distr(generator);
                x[i] = r * cos(phi);
                y[i] = r * sin(phi);
            }
        }

};


std::vector<double> linspace(double min, double max, int nr_points){
    std::vector<double> out(nr_points);
    double delta = (max - min)/(1.0 * (nr_points - 1));

    for (int i = 0; i < nr_points; i++){
        out[i] = min + delta*i;
    }
    return out;
}


std::vector<double> point_source(std::vector<double> z) {
    int size = z.size();
    std::vector<double> ps(size);
    double temp;

    for (int i = 0; i < size; i++) {
        temp = 50 - (50 * z[i]) / (sqrt(1 + pow(z[i], 2)));
        ps[i] = temp;
    }
    return ps;
}


double geom_eff_point(double z, double source, int n, int seed, std::string type){
    int N_hit = 0;
    std::vector<double> x1(n), x2(n), y1(n), y2(n);
    position generate_source(x1, y1);
    position generate_emission(x2, y2);

    if (type == "circular"){
        generate_source.generate_circular_distr(source, seed);
    } else if (type == "gaussian"){
        generate_source.generate_gaussian_distr(source, seed);
    } else{
        std::cerr << "ERROR: Not a valid source type. Choose 'circular' or 'gaussian'" << std::endl;
        exit(0);
    }
    
    seed++;
    generate_emission.generate_isotropic(z, seed);
    generate_source.add_vec(generate_emission.x, generate_emission.y);
    std::vector<double> r_final = generate_source.calculate_rsq();

    for (int i = 0; i < n; i++){
        if (r_final[i] <= 1){
            N_hit++;
        }
    }
    return 50.0*N_hit/n;
}


void write_geo_file(std::vector<double> z, std::vector<double> efficiencies, std::vector<double> rel_ers, std::string filename) {
    std::vector<double> e_ps = point_source(z);
    std::ofstream myFile(filename);
    myFile << "z/rd \t point source \t Model \t Relative uncertainty \n";

    for (int i = 0; i < z.size(); i++) {
        myFile << z[i] << "\t" << e_ps[i] << "\t" << efficiencies[i] << "\t" << rel_ers[i]<< "\n";
    }
    
    std::cout << "Wrote output file" << std::endl;
    myFile.close();
}


int main(int argc, char **argv){
    int seed = 15763027;
    std::string type;
    double z_min, z_max, source, result;
    int n_points, power;
    std::string filename;

    if (argc != 2){
        std::cerr << "ERROR: input option 'circular' or 'gaussian' for the source distribution" << std::endl;
        exit(0);
    } else{
        type = argv[1];
    }

    std::cout << "z_min/rd:" << std::endl;
    std::cin >> z_min;
    std::cout << "z_max/rd:" << std::endl;
    std::cin >> z_max;
    std::cout << "number of points:" << std::endl;
    std::cin >> n_points;
    std::cout << "source/rd:" << std::endl;
    std::cin >> source;
    std::cout << "Power:" << std::endl;
    std::cin >> power;
    std::cout << "Filename:" << std::endl;
    std::cin >> filename;

    int n_perpoint = pow(10, power);
    std::vector<double> z = linspace(z_min, z_max, n_points);
    std::vector<double> efficiencies(n_points), rel_ers(n_points);

    for (int i = 0; i < n_points; i++){
        efficiencies[i] = geom_eff_point(z[i], source, n_perpoint, seed, type);
        rel_ers[i] = 100 / sqrt(efficiencies[i] * n_perpoint/2);
        std::cout << z[i]/z[n_points-1] << "\t" << efficiencies[i] << "\t" << rel_ers[i] << std::endl; 
    }

    write_geo_file(z, efficiencies, rel_ers, filename);
    return 1;
}
