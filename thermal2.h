// OFHC thermal parameters
const double cp_cu = 0.091;  // J/kg/K
const double k_cu = 156;  // W/m/K  (RRR = 100)
const double d_cu = 8960;  // kg/m^3

// helium-3 thermal parameters
double cp_3he(double T) {
    return 561.57 * pow(T, 2) - 284.09 * T + 1131.9;
}

double k_3he(double T) {
    return (3.1403 * T + 5.4571) / 1000.0;
}

double d_3he(double T) {
    return -2.2451 * pow(T, 2) + 3.0238 * T + 81.05;
}

const double kg_3he = 40;
double k_kap_3he(double T) {
	return kg_3he * 20 * pow(T, 3);
}

// helium-4 thermal parameters
double cp_4he(double T) {
    return 4908.3 * pow(T, 2) - 9973.3 * T + 5102.6;
}

double d_4he(double T) {
    return 0.5724 * pow(T, 2) - 1.2547 * T + 145.83;
}

const double kg_4he = 40;
double k_kap_4he(double T) {
	return kg_4he * 20 * pow(T, 3);
}

const double T_lambda = 2.1768;
double f_inv(double T) {
    return 9.469 * pow(10, 14) * pow(pow((1.2/2.1768), 5.7) * (1 - pow((1.2/2.1768), 5.7)), 3);
}

