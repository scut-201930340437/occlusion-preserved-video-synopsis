#pragma once
#include "defines.h"
#include <iostream>

void read_group(std::string filepath, std::map<int, std::vector<int>>& group, int tube_num);
void fusion(Context& context, std::map<int, int>& tubes_group_id,double speed);
double cal_colli_area(cv::Vec6d& a, cv::Vec6d& b);