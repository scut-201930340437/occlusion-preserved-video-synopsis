#pragma once
#include <vector>
#include <algorithm>
#include "segment.h"
#include "defines.h"

//double compute_energy_outside(Context & context, std::vector<std::vector<Segment*>> & all_segs, int synoplen, std::vector<std::vector<double>>& activity);

double compute_energy_outside(Context& context);
double compute_energy_chrono_reverse(std::vector<Segment*> & fronts);
double compute_energy_speed(std::vector<std::vector<Segment*>> & all_segs);
double compute_energy(Context& context, std::vector<double>& tube_energy, double& enout, double& encolli, double& enscale, double & enself);
double compute_energy_chrono_reverse3(Context& context);
double compute_energy_chrono_reverse4(Context& context);
double compute_energy_collision3(Context& context, std::vector<double>& tube_energy,
	std::vector<std::vector<std::tuple<cv::Vec6d, int, int>>>& framebb,
	std::vector<int>& frame_count);
double compute_energy_collision4(Context& context, std::vector<double>& tube_energy, double& cost_circuit_self);
double compute_energy_collision4_vis(Context& context);
double compute_energy_collision5(Context& context, std::vector<double>& tube_energy, double& cost_circuit_self);
double compute_energy_real(Context& context, std::vector<double>& tube_energy, double& enout, double& encolli, double& enscale, double& enself);