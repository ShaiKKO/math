// Copyright (c) 2025 Boost.Math Contributors
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "benchmark_base.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <iomanip>
#include <cmath>

namespace bmb = boost::math::benchmark;

// CSV parser for benchmark results
class benchmark_reader {
public:
    static std::vector<bmb::benchmark_result> read_csv(const std::string& filename) {
        std::vector<bmb::benchmark_result> results;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        std::string line;
        // Skip header
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            bmb::benchmark_result result;
            if (parse_csv_line(line, result)) {
                results.push_back(result);
            }
        }
        
        file.close();
        return results;
    }
    
private:
    static bool parse_csv_line(const std::string& line, bmb::benchmark_result& result) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < 10) {
            return false;
        }
        
        try {
            result.method_name = tokens[0];
            result.function_name = tokens[1];
            result.dimension = std::stoull(tokens[2]);
            result.tolerance = std::stod(tokens[3]);
            result.execution_time = std::stod(tokens[4]);
            result.relative_error = std::stod(tokens[5]);
            result.function_calls = std::stoull(tokens[6]);
            result.memory_usage = std::stod(tokens[7]);
            result.efficiency_score = std::stod(tokens[8]);
            result.converged = (tokens[9] == "true");
            return true;
        } catch (...) {
            return false;
        }
    }
};

// Analysis results structure
struct analysis_results {
    struct method_summary {
        std::string name;
        std::size_t total_runs;
        std::size_t successful_runs;
        double success_rate;
        double avg_time;
        double avg_error;
        double avg_efficiency;
        double best_efficiency;
        double worst_efficiency;
        std::size_t total_evaluations;
        
        // Dimension-specific performance
        std::map<std::size_t, double> avg_time_by_dim;
        std::map<std::size_t, double> avg_error_by_dim;
        
        // Function-specific performance
        std::map<std::string, double> avg_time_by_func;
        std::map<std::string, double> avg_error_by_func;
    };
    
    std::map<std::string, method_summary> methods;
    std::set<std::size_t> dimensions;
    std::set<std::string> functions;
    std::set<double> tolerances;
    
    // Best performers
    std::map<std::string, std::string> best_by_function;
    std::map<std::size_t, std::string> best_by_dimension;
    std::map<double, std::string> best_by_tolerance;
};

// Comprehensive analyzer
class results_analyzer {
public:
    static analysis_results analyze(const std::vector<bmb::benchmark_result>& results) {
        analysis_results analysis;
        
        // Collect unique values
        for (const auto& r : results) {
            analysis.dimensions.insert(r.dimension);
            analysis.functions.insert(r.function_name);
            analysis.tolerances.insert(r.tolerance);
        }
        
        // Group by method
        std::map<std::string, std::vector<const bmb::benchmark_result*>> by_method;
        for (const auto& r : results) {
            by_method[r.method_name].push_back(&r);
        }
        
        // Analyze each method
        for (const auto& [method, method_results] : by_method) {
            analysis_results::method_summary summary;
            summary.name = method;
            summary.total_runs = method_results.size();
            
            // Count successful runs
            summary.successful_runs = std::count_if(
                method_results.begin(), method_results.end(),
                [](const auto* r) { return r->converged; }
            );
            summary.success_rate = double(summary.successful_runs) / summary.total_runs;
            
            // Compute averages
            double total_time = 0, total_error = 0, total_efficiency = 0;
            summary.best_efficiency = 0;
            summary.worst_efficiency = std::numeric_limits<double>::max();
            summary.total_evaluations = 0;
            
            for (const auto* r : method_results) {
                if (r->converged) {
                    total_time += r->execution_time;
                    total_error += r->relative_error;
                    total_efficiency += r->efficiency_score;
                    summary.total_evaluations += r->function_calls;
                    
                    summary.best_efficiency = std::max(summary.best_efficiency, r->efficiency_score);
                    if (r->efficiency_score > 0) {
                        summary.worst_efficiency = std::min(summary.worst_efficiency, r->efficiency_score);
                    }
                    
                    // By dimension
                    summary.avg_time_by_dim[r->dimension] += r->execution_time;
                    summary.avg_error_by_dim[r->dimension] += r->relative_error;
                    
                    // By function
                    summary.avg_time_by_func[r->function_name] += r->execution_time;
                    summary.avg_error_by_func[r->function_name] += r->relative_error;
                }
            }
            
            if (summary.successful_runs > 0) {
                summary.avg_time = total_time / summary.successful_runs;
                summary.avg_error = total_error / summary.successful_runs;
                summary.avg_efficiency = total_efficiency / summary.successful_runs;
                
                // Normalize by-dimension averages
                std::map<std::size_t, std::size_t> count_by_dim;
                std::map<std::string, std::size_t> count_by_func;
                
                for (const auto* r : method_results) {
                    if (r->converged) {
                        count_by_dim[r->dimension]++;
                        count_by_func[r->function_name]++;
                    }
                }
                
                for (auto& [dim, time] : summary.avg_time_by_dim) {
                    time /= count_by_dim[dim];
                    summary.avg_error_by_dim[dim] /= count_by_dim[dim];
                }
                
                for (auto& [func, time] : summary.avg_time_by_func) {
                    time /= count_by_func[func];
                    summary.avg_error_by_func[func] /= count_by_func[func];
                }
            }
            
            analysis.methods[method] = summary;
        }
        
        // Find best performers for each category
        for (const auto& func : analysis.functions) {
            double best_efficiency = 0;
            std::string best_method;
            
            for (const auto& r : results) {
                if (r.function_name == func && r.converged && 
                    r.efficiency_score > best_efficiency) {
                    best_efficiency = r.efficiency_score;
                    best_method = r.method_name;
                }
            }
            
            if (!best_method.empty()) {
                analysis.best_by_function[func] = best_method;
            }
        }
        
        // Best by dimension
        for (const auto& dim : analysis.dimensions) {
            double best_efficiency = 0;
            std::string best_method;
            
            for (const auto& r : results) {
                if (r.dimension == dim && r.converged && 
                    r.efficiency_score > best_efficiency) {
                    best_efficiency = r.efficiency_score;
                    best_method = r.method_name;
                }
            }
            
            if (!best_method.empty()) {
                analysis.best_by_dimension[dim] = best_method;
            }
        }
        
        // Best by tolerance
        for (const auto& tol : analysis.tolerances) {
            double best_efficiency = 0;
            std::string best_method;
            
            for (const auto& r : results) {
                if (std::abs(r.tolerance - tol) < 1e-10 && r.converged && 
                    r.efficiency_score > best_efficiency) {
                    best_efficiency = r.efficiency_score;
                    best_method = r.method_name;
                }
            }
            
            if (!best_method.empty()) {
                analysis.best_by_tolerance[tol] = best_method;
            }
        }
        
        return analysis;
    }
    
    static void print_summary(const analysis_results& analysis) {
        std::cout << "\n================================================\n";
        std::cout << "         BENCHMARK ANALYSIS SUMMARY\n";
        std::cout << "================================================\n\n";
        
        // Overall statistics
        std::cout << "Dataset Statistics:\n";
        std::cout << "  Dimensions tested: ";
        for (auto d : analysis.dimensions) std::cout << d << " ";
        std::cout << "\n  Functions tested: ";
        for (const auto& f : analysis.functions) std::cout << f << " ";
        std::cout << "\n  Tolerances tested: ";
        for (auto t : analysis.tolerances) std::cout << std::scientific << t << " ";
        std::cout << "\n\n";
        
        // Method performance summary
        std::cout << "Method Performance Summary:\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << std::setw(15) << "Method" 
                  << std::setw(10) << "Success"
                  << std::setw(12) << "Avg Time"
                  << std::setw(12) << "Avg Error"
                  << std::setw(15) << "Avg Efficiency"
                  << std::setw(12) << "Total Evals\n";
        std::cout << std::string(80, '-') << "\n";
        
        for (const auto& [name, summary] : analysis.methods) {
            std::cout << std::setw(15) << name
                      << std::setw(9) << std::fixed << std::setprecision(1) 
                      << (summary.success_rate * 100) << "%"
                      << std::setw(12) << std::fixed << std::setprecision(4) 
                      << summary.avg_time
                      << std::setw(12) << std::scientific << std::setprecision(2) 
                      << summary.avg_error
                      << std::setw(15) << std::scientific << std::setprecision(2) 
                      << summary.avg_efficiency
                      << std::setw(12) << summary.total_evaluations << "\n";
        }
        
        // Best performers
        std::cout << "\n\nBest Performers by Category:\n";
        std::cout << std::string(50, '-') << "\n";
        
        std::cout << "\nBy Function Type:\n";
        for (const auto& [func, method] : analysis.best_by_function) {
            std::cout << "  " << std::setw(20) << func << ": " << method << "\n";
        }
        
        std::cout << "\nBy Dimension:\n";
        for (const auto& [dim, method] : analysis.best_by_dimension) {
            std::cout << "  Dimension " << std::setw(3) << dim << ": " << method << "\n";
        }
        
        std::cout << "\nBy Tolerance:\n";
        for (const auto& [tol, method] : analysis.best_by_tolerance) {
            std::cout << "  Tolerance " << std::scientific << std::setprecision(0) 
                      << tol << ": " << method << "\n";
        }
    }
    
    static void generate_performance_report(const analysis_results& analysis,
                                           const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        file << "# Cubature Benchmark Performance Report\n\n";
        file << "Generated: " << get_timestamp() << "\n\n";
        
        // Executive Summary
        file << "## Executive Summary\n\n";
        
        // Find overall best method
        std::string best_overall;
        double best_avg_efficiency = 0;
        for (const auto& [name, summary] : analysis.methods) {
            if (summary.avg_efficiency > best_avg_efficiency) {
                best_avg_efficiency = summary.avg_efficiency;
                best_overall = name;
            }
        }
        
        file << "**Best Overall Method**: " << best_overall 
             << " (avg efficiency: " << std::scientific << std::setprecision(2) 
             << best_avg_efficiency << ")\n\n";
        
        // Detailed Performance Analysis
        file << "## Detailed Performance Analysis\n\n";
        
        for (const auto& [name, summary] : analysis.methods) {
            file << "### " << name << "\n\n";
            file << "- **Success Rate**: " << std::fixed << std::setprecision(1) 
                 << (summary.success_rate * 100) << "%\n";
            file << "- **Average Execution Time**: " << std::fixed << std::setprecision(4) 
                 << summary.avg_time << " seconds\n";
            file << "- **Average Relative Error**: " << std::scientific << std::setprecision(2) 
                 << summary.avg_error << "\n";
            file << "- **Total Function Evaluations**: " << summary.total_evaluations << "\n";
            file << "- **Efficiency Range**: [" << std::scientific << std::setprecision(2) 
                 << summary.worst_efficiency << ", " << summary.best_efficiency << "]\n\n";
            
            // Dimension scaling
            file << "#### Performance by Dimension\n\n";
            file << "| Dimension | Avg Time (s) | Avg Error |\n";
            file << "|-----------|--------------|----------|\n";
            for (const auto& [dim, time] : summary.avg_time_by_dim) {
                file << "| " << std::setw(9) << dim 
                     << " | " << std::fixed << std::setprecision(4) << time
                     << " | " << std::scientific << std::setprecision(2) 
                     << summary.avg_error_by_dim.at(dim) << " |\n";
            }
            file << "\n";
            
            // Function-specific performance
            file << "#### Performance by Function Type\n\n";
            file << "| Function | Avg Time (s) | Avg Error |\n";
            file << "|----------|--------------|----------|\n";
            for (const auto& [func, time] : summary.avg_time_by_func) {
                file << "| " << std::setw(15) << func 
                     << " | " << std::fixed << std::setprecision(4) << time
                     << " | " << std::scientific << std::setprecision(2) 
                     << summary.avg_error_by_func.at(func) << " |\n";
            }
            file << "\n";
        }
        
        // Recommendations
        file << "## Recommendations\n\n";
        file << "Based on the benchmark results, we recommend:\n\n";
        
        // Low-dimensional problems
        std::string best_low_dim;
        for (const auto& [dim, method] : analysis.best_by_dimension) {
            if (dim <= 5) {
                best_low_dim = method;
                break;
            }
        }
        if (!best_low_dim.empty()) {
            file << "1. **Low-dimensional problems (d ≤ 5)**: Use " << best_low_dim << "\n";
        }
        
        // High-dimensional problems
        std::string best_high_dim;
        for (const auto& [dim, method] : analysis.best_by_dimension) {
            if (dim >= 15) {
                best_high_dim = method;
            }
        }
        if (!best_high_dim.empty()) {
            file << "2. **High-dimensional problems (d ≥ 15)**: Use " << best_high_dim << "\n";
        }
        
        // High accuracy requirements
        std::string best_high_acc;
        for (const auto& [tol, method] : analysis.best_by_tolerance) {
            if (tol <= 1e-6) {
                best_high_acc = method;
                break;
            }
        }
        if (!best_high_acc.empty()) {
            file << "3. **High accuracy requirements (ε ≤ 1e-6)**: Use " << best_high_acc << "\n";
        }
        
        file << "\n## Conclusion\n\n";
        file << "The benchmark suite demonstrates that different integration methods excel in ";
        file << "different scenarios. The choice of method should be guided by:\n\n";
        file << "- Problem dimension\n";
        file << "- Required accuracy\n";
        file << "- Function smoothness\n";
        file << "- Available computational resources\n";
        
        file.close();
    }
    
private:
    static std::string get_timestamp() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }
};

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <benchmark_results.csv> [output_report.md]\n";
            return 1;
        }
        
        std::string input_file = argv[1];
        std::string output_file = "performance_report.md";
        if (argc > 2) {
            output_file = argv[2];
        }
        
        std::cout << "Reading benchmark results from: " << input_file << "\n";
        
        // Read benchmark results
        auto results = benchmark_reader::read_csv(input_file);
        std::cout << "Loaded " << results.size() << " benchmark results\n";
        
        if (results.empty()) {
            std::cerr << "No results found in file\n";
            return 1;
        }
        
        // Analyze results
        auto analysis = results_analyzer::analyze(results);
        
        // Print summary to console
        results_analyzer::print_summary(analysis);
        
        // Generate detailed report
        std::cout << "\nGenerating detailed report to: " << output_file << "\n";
        results_analyzer::generate_performance_report(analysis, output_file);
        
        // Generate ranking
        std::cout << "\n================================================\n";
        std::cout << "         OVERALL METHOD RANKING\n";
        std::cout << "================================================\n\n";
        
        // Sort methods by average efficiency
        std::vector<std::pair<std::string, double>> ranking;
        for (const auto& [name, summary] : analysis.methods) {
            ranking.push_back({name, summary.avg_efficiency});
        }
        std::sort(ranking.begin(), ranking.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        std::cout << std::setw(5) << "Rank" 
                  << std::setw(20) << "Method"
                  << std::setw(20) << "Avg Efficiency\n";
        std::cout << std::string(45, '-') << "\n";
        
        int rank = 1;
        for (const auto& [method, efficiency] : ranking) {
            std::cout << std::setw(5) << rank++
                      << std::setw(20) << method
                      << std::setw(20) << std::scientific << std::setprecision(2) 
                      << efficiency << "\n";
        }
        
        std::cout << "\n================================================\n";
        std::cout << "Analysis complete. Report saved to: " << output_file << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}