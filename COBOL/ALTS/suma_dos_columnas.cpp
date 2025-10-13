#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

int main() {
    double total_col1 = 0.0;
    double total_col2 = 0.0;
    std::ifstream file("DATOS.DAT");

    if (!file.is_open()) {
        std::cerr << "Error al abrir DATOS.DAT" << std::endl;
        return 1;
    }

    double num_col1, num_col2;
    while (file >> num_col1 >> num_col2) {
        total_col1 += num_col1;
        total_col2 += num_col2;
    }

    file.close();

    std::cout << "TOTAL COLUMNA 1: " << std::fixed << std::setprecision(2) << std::setw(15) << total_col1 << std::endl;
    std::cout << "TOTAL COLUMNA 2: " << std::fixed << std::setprecision(2) << std::setw(15) << total_col2 << std::endl;

    return 0;
}
