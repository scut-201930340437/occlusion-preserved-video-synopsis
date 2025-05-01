#pragma once
#include <vector>
#include <string>
#include <opencv2/opencv.hpp>
#include "segment.h"
#include "defines.h"
using namespace std;

void visualize(Context & context, std::map<int,int>&tubes_group_id, double speed, int chrono_slack);
void visualize_for_ours(Context& context, int chrono_slack);
void visualize_for_ours_step(Context& context);
void visualize_root(string filepath, std::vector<Segment*>& segs, std::vector<std::vector<cv::Vec6d>>& bbox, double scale,
	int root_id, int video_len);
void visualize_A_frame(cv::Mat& frame, int frame_id, Context & context);
cv::Rect get_bbox_of_obj_by_frame(std::vector<std::vector<cv::Vec6d>>& bbox, int tubeid, int frameid);

void visualize2(Context& context, std::vector<std::vector<Segment*>> best);
void visualize1(Context& context, std::vector<std::vector<Segment*>> best);