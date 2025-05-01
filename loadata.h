#pragma once
#include "segment.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <opencv2/opencv.hpp>

using namespace std;

void read_activity(string filepath, std::vector<std::vector<double>>& activity, int tube_num);
void read_boundingbox(string filepath, std::vector<std::vector<cv::Vec6d>>& bbox, int tube_num, double scale, int video_height,
	int video_width);
void read_bestCopy(string filepath, std::vector<std::vector<Segment*>>& all_segs);
void read_occlu(string filepath, std::unordered_map<int, std::vector<cv::Point3i>>& orig_colli, std::unordered_map<int, std::set<int>>& occlu_graph);