#pragma once

#include "defines.h"
#include "segment.h"
#include <windows.h>
#include<string>

int random_bbb(cv::RNG& rng, Segment* seg, bool follow);
int random_aaa(cv::RNG& rng, Segment* seg, bool follow);
void mcmc(Context& context, std::vector<std::vector<Segment*>>& bestCopy);
double mcmc2(Context& context, int& iter_best, double& min_s, double& max_s, bool can_preserve_order);
double mcmc_chro_colli(Context& context, int& iter_best, double& min_s, double& max_s,double limit_time);
//void mcmc_chrono(Context& context, std::vector<std::vector<Segment*>>& bestCopy, double& min_s, double& max_s);
bool allin_range(Context& context);
bool chronological_order_preserved(Context& context, int slack);
double mcmc3(Context& context, std::vector<std::vector<Segment*>>& bestCopy, int& iter_best, double& min_s, double& max_s);
int mcmc_chrono(Context& context, std::vector<std::vector<Segment*>>& bestCopy, double& min_s, double& max_s);
//int mcmc_beg_chrono(Context& context, std::vector<std::vector<Segment*>>& bestCopy, double& min_s, double& max_s);
LPCWSTR stringToLPCWSTR(std::string orig);