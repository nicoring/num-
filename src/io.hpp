#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "matrix.hpp"

namespace num {

template <typename T> std::vector<T> parse_row(std::string line, const char delimiter) {
    std::istringstream line_stream(line);
    std::string cell;
    T value;
    std::vector<T> result;
    while (std::getline(line_stream, cell, delimiter)) {
        std::istringstream stream(cell);
        stream >> value;
        result.push_back(value);
    }
    return result;
}

template <typename T>
matrix<T> read_csv(const std::string filename, const char delimiter = ',',
                   const bool has_header = true) {
    std::vector<std::vector<T>> data;

    std::string line;
    std::ifstream input_file(filename);
    if (input_file.is_open()) {
        if (has_header) {
            input_file.ignore(1);
        }
        while (std::getline(input_file, line)) {
            std::vector<T> parsed_row = parse_row<T>(line, delimiter);
            data.push_back(parsed_row);
        }
        input_file.close();
    } else {
        std::cout << "Unable to open file" << std::endl;
    }
    matrix<T> mat(data);
    return mat;
}
} // namespace num
