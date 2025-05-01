#pragma once
#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <opencv2/opencv.hpp>

class Segment;
struct EdgeNode;


struct Context
{
	cv::Mat segCollisionMat;
	std::vector<int> collisionVec;
	std::string filepath;
	int seg_len;
	int seg_avg_len = 0;
	std::vector<Segment*> heads;
	std::vector < std::vector<Segment*>> all_segs;
	std::vector<std::pair<EdgeNode*, EdgeNode*>> all_edges;

	int tube_num;

	std::vector<int> circulated_tube;

	std::vector<std::vector<double>> activity;
	std::vector<std::vector<cv::Vec6d>> bbox;

	std::unordered_map<int, std::vector<cv::Point3i>>orig_colli;

	std::unordered_map<int, std::set<int>> occlu_graph;

	std::unordered_map<int,std::unordered_map<int,double>>tube_speed;

	std::vector<Segment*> roots;

	std::vector<bool>visualize_state;

	std::map<int, std::vector<std::pair<int, cv::Vec6d>>>bbox_per_frame; // for compared methods

	int synoplen;
	int srclen;
	double scale;
	int chrono_slack;
	int seg_len_lbound;
	double fps;
	int algo_type;
	int mcmc2_iteration_num;
	bool resizing_component;
	bool remainspeed_component;
	bool non_colliding_divide_component;
	double sigma;
	double sigma_decay;
	int video_width;
	int video_height;
};

void release_context(Context& context);




inline int mymin(int a, int b)
{
	return a < b ? a : b;
}

inline int mymax(int a, int b)
{
	return a > b ? a : b;
}

inline int mymin(int a, int b, int c)
{
	return mymin(mymin(a, b), c);
}

inline int mymax(int a, int b, int c)
{
	return mymax(mymax(a, b), c);
}



inline double mymin(double a, double b)
{
	return a < b ? a : b;
}

inline double mymax(double a, double b)
{
	return a > b ? a : b;
}

inline double mymin(double a, double b, double c)
{
	return mymin(mymin(a, b), c);
}

inline double mymax(double a, double b, double c)
{
	return mymax(mymax(a, b), c);
}