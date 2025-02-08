// NACA Airfoils.h : Include file for standard system include files,
// or project specific include files.

#ifndef NACA_AIRFOILS_HPP
#define NACA_AIRFOILS_HPP

#include <iostream>
#include <concepts>
#include <vector>
#include <cmath>
#include <numbers>
#include <charconv>
#include <string>

namespace naca{

	enum class airfoil_spacing {
		linear,
		cosine
	};
	enum class airfoil_trailing_edge {
		closed,
		open
	};

	template<typename precision>
	requires std::floating_point<precision>
	class airfoil final {
	public:
		constexpr explicit airfoil(std::string_view input_model, std::size_t input_total_points, airfoil_spacing input_spacing, airfoil_trailing_edge input_trailing_edge) noexcept :
			model{ input_model },
			points(input_total_points + !(input_total_points % 2)),
			spacing{ input_spacing },
			trailing_edge{ input_trailing_edge } {
			calculate_coordinates();
		}
		constexpr void set_total_points(std::size_t new_size) noexcept {
			points.resize(new_size);
		}
		constexpr void set_model(std::string_view input_model) noexcept {
			model = input_model;
		}

		constexpr std::string_view get_model() const noexcept {
			return model;
		}
		constexpr std::size_t get_total_points() const noexcept {
			return points.size();
		}

		struct point;
		constexpr const std::vector<point>& get_points() const noexcept {
			return points;
		}
	private:
		std::string model;
		struct point { precision x, y; };
		std::vector<point> points;
		airfoil_spacing spacing;
		airfoil_trailing_edge trailing_edge;
		constexpr void calculate_coordinates() noexcept {

			precision M{}, P{}, T{};
			std::from_chars(model.data(), model.data() + 1, M);
			std::from_chars(model.data() + 1, model.data() + 2, P);
			std::from_chars(model.data() + 2, model.data() + model.size(), T);
			M *= static_cast<precision>(0.01), P *= static_cast<precision>(0.1), T *= static_cast<precision>(0.01);

			for (std::size_t counter{ 0 }; counter < static_cast<std::size_t>(points.size() * 0.5); ++counter) {

				const precision x{ spacing == airfoil_spacing::cosine ?
				static_cast<precision>((1 - std::cos(std::numbers::pi_v<precision> -counter * (std::numbers::pi_v<precision> / ((points.size() - 1) * 0.5)))) * 0.5) :
				static_cast<precision>(1 - counter / ((points.size() - 1) * 0.5)) };

				const precision yc{ M / (P * P) * (2 * P * x - x * x) * static_cast<precision>(x < P) + static_cast<precision>(x >= P) * M / ((1 - P) * (1 - P)) * (1 - 2 * P + 2 * P * x - x * x) };
				const precision gradient{ 2 * M / (P * P) * (P - x) * static_cast<precision>(x < P) + static_cast<precision>(x >= P) * 2 * M / ((1 - P) * (1 - P)) * (P - x) };
				const precision yt{ static_cast<precision>(5 * T * (0.2969 * std::sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - (0.1036 * static_cast<precision>(trailing_edge == airfoil_trailing_edge::closed) + 0.1015 * static_cast<precision>(trailing_edge == airfoil_trailing_edge::open)) * x * x * x * x)) };

				const precision theta{ std::atan(gradient) };
				points[counter] = { x + yt * std::sin(theta), yc - yt * std::cos(theta) };
				points[points.size() - counter - 1] = { x - yt * std::sin(theta),  yc + yt * std::cos(theta) };
			}

			points[(points.size() - 1) * 0.5] = { 0, 0 };
		}
	};
}

#endif

// TODO: Reference additional headers your program requires here.
