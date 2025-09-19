// include/phatpsy/core/parameters.hpp
#pragma once
#include <map>
#include <string>
#include <variant>

namespace phatpsy {

using ParameterValue = std::variant<bool, int, double, std::string>;

class PhatpsyParameters {
private:
    std::map<std::string, ParameterValue> params_;
    
public:
    PhatpsyParameters();
    
    // Parameter access
    template<typename T>
    T get(const std::string& key) const {
        auto it = params_.find(key);
        if (it == params_.end()) {
            throw std::runtime_error("Parameter " + key + " not found");
        }
        return std::get<T>(it->second);
    }
    
    template<typename T>
    void set(const std::string& key, const T& value) {
        params_[key] = value;
    }
    
    // PHATPSY-specific parameters
    bool debug_mode() const { return get<bool>("debug"); }
    bool unrestricted() const { return get<bool>("unrestricted"); }
    bool use_ewmo() const { return get<bool>("use_ewmo"); }
    int max_scf_iterations() const { return get<int>("max_scf_iter"); }
    double scf_convergence() const { return get<double>("scf_convergence"); }
    double fermi_temperature() const { return get<double>("fermi_temperature"); }
    
    void set_defaults();
    void read_from_file(const std::string& filename);
    void write_to_file(const std::string& filename) const;
};

} // namespace phatpsy
