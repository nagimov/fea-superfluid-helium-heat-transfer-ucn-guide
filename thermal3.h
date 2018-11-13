// OFHC thermal parameters
const double cp_cu = 0.091;  // J/kg/K
const double k_cu = 156;  // W/m/K  (RRR = 100)
const double d_cu = 8960;  // kg/m^3

// helium-3 thermal parameters
double cp_3he(double T) {
    return 283.73 * pow(T, 3) - 502.41 * pow(T, 2) + 948.86 * T + 699.05;
}

double k_3he(double T) {
    return (-0.1021 * pow(T, 3) + 0.9494 * pow(T, 2) + 1.2804 * T + 6.3904) / 1000.0;
}

double d_3he(double T) {
    return -0.6514 * pow(T, 3) + 0.1976 * pow(T, 2) + 0.1933 * T + 82.044;
}

const double kg_3he = 40;
double k_kap_3he(double T) {
	return kg_3he * 20 * pow(T, 3);
}

// helium-4 thermal parameters
double cp_4he(double T) {
    return 3380.9 * pow(T, 3) - 9287.4 * pow(T, 2) + 9151.7 * T - 3126.1;
}

double d_4he(double T) {
    return 0.3617 * pow(T, 3) - 0.9466 * pow(T, 2) + 0.7926 * T + 144.95;
}

const double kg_4he = 40;
double k_kap_4he(double T) {
	return kg_4he * 20 * pow(T, 3);
}

const double T_lambda = 2.1768;
double f_inv(double T) {
    return 9.469 * pow(10, 14) * pow(pow((1.2/2.1768), 5.7) * (1 - pow((1.2/2.1768), 5.7)), 3);
}

