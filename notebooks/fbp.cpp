#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include "FBP5.0.h" // Make sure to include your C code header file

// Function prototypes
void calculate_firespread(const std::vector<std::vector<std::string>>& data, std::vector<std::vector<float>>& results);
std::vector<std::vector<std::string>> read_csv(const std::string& filename);
void write_csv(const std::string& filename, const std::vector<std::vector<float>>& results);

int main() {
    // Read the CSV file
    std::string input_filename = "/Users/minho/Documents/GitHub/Cell2Fire/notebooks/DataGenerated_Small.csv";
    std::vector<std::vector<std::string>> data = read_csv(input_filename);
    
    // Prepare results container
    std::vector<std::vector<float>> results;

    // Calculate rates of spread
    calculate_firespread(data, results);
    
    // Write the results to CSV file
    std::string output_filename = "/Users/minho/Documents/GitHub/Cell2Fire/notebooks/DataGenerated_Small_ros.csv";
    write_csv(output_filename, results);

    return 0;
}

// Read CSV file and return data as a vector of vectors
std::vector<std::vector<std::string>> read_csv(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<std::string>> data;
    std::string line;

    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }
        data.push_back(row);
    }

    file.close();
    return data;
}

// Write results to CSV file
void write_csv(const std::string& filename, const std::vector<std::vector<float>>& results) {
    std::ofstream file(filename);

    file << "Head Rate of Spread,Flank Rate of Spread,Back Rate of Spread\n";
    for (const auto& row : results) {
        file << std::fixed << std::setprecision(2) << row[0] << "," << row[1] << "," << row[2] << "\n";
    }

    file.close();
}

// Calculate head, flank, and back rate of spread
void calculate_firespread(const std::vector<std::vector<std::string>>& data, std::vector<std::vector<float>>& results) {
    fuel_coefs fuels[18];
    setup_const(fuels);

    for (const auto& row : data) {
        if (row.size() < 16) continue; // Skip incomplete rows

        inputs inp;
        main_outs out;
        snd_outs sec;
        fire_struc head, flank, back;

        // Populate inputs structure from CSV data
        inp.fueltype = row[0];
        inp.ffmc = std::stof(row[8]);
        inp.ws = std::stof(row[9]);
        inp.waz = std::stof(row[10]);
        inp.bui = std::stof(row[11]);
        inp.ps = std::stof(row[12]);
        inp.saz = std::stof(row[13]);
        inp.pc = std::stof(row[14]);
        inp.pdf = std::stof(row[15]);
        inp.gfl = 0; // Default
        inp.cur = 0; // Default
        inp.time = 0; // Default
        inp.pattern = 0; // Default
        inp.jd = 0; // Default
        inp.jd_min = 0; // Default
        inp.lat = std::stof(row[5]);
        inp.lon = std::stof(row[6]);
        inp.elev = std::stof(row[7]);

        // Initialize pointers to fire structures
        fire_struc *hptr = &head;
        fire_struc *fptr = &flank;
        fire_struc *bptr = &back;

        // Call calculate function
        calculate(&inp, fuels, &out, &sec, hptr, fptr, bptr);

        // Store results
        results.push_back({
            hptr->rss, // Head rate of spread
            fptr->rss, // Flank rate of spread
            bptr->rss  // Back rate of spread
        });
    }
}
