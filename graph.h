#pragma once
#include <vector>
#include "segment.h"
#include <string>
#include <fstream>
#include <stack>
#include "mcmc.h"
#include "defines.h"

void check_consistency(Context& context);

void reset_traverse(std::vector<std::vector<Segment*>>& all_segs);
void reset_traverse(std::vector<Segment*> & heads);
void reset_last_next_break(std::vector<Segment*>& heads);
void reset_colli_break(std::vector<Segment*>& heads);

//void dfs_change_aabb(Segment* seg, int oaa, int obb);
void bfs_change_aabb(Context& context, Segment* seg, int oaa, int obb, bool show);
void dfs_change_aabb(Segment* seg, int oaa, int obb, bool print_self = false);

cv::Mat create_segment_graph(std::string filename, std::vector<std::vector<Segment*>>& all_segs,
	std::vector<std::pair<EdgeNode*, EdgeNode*>>& all_edges,
	int seg_len,
	int & srclen);

cv::Mat create_segment_graph2(std::string filename, 
	std::vector<Segment*>& heads,
	std::vector<std::vector<Segment*>>& all_segs,
	std::vector<std::pair<EdgeNode*, EdgeNode*>>& all_edges,
	int seg_len,
	int& srclen, 
	int algo_type,
	bool non_colliding_dividing_component,
	std::unordered_map<int, std::vector<cv::Point3i>>&orig_colli);

//cv::Mat create_segment_simple_graph(std::string filename,
//	std::vector<Segment*>& heads,
//	std::vector<std::vector<Segment*>>& all_segs,
//	std::vector<std::pair<EdgeNode*, EdgeNode*>>& all_edges,
//	int seg_len,
//	int& srclen);

void dfs_connected_graph_ofA(Segment* A, std::vector<Segment*>& all);
void dfs_jump_chain(Segment* A, std::vector<std::pair<int, int>>& chain);
bool check_circulation(std::vector<std::pair<int, int>>& chain);
void bfs_traverseMST_and_show(Context& context, Segment* seg, bool show,
	std::vector<std::pair<Segment*, Segment*>>& influ_pairs);

void dfs_traverse(Segment* seg);


void extract_segment_from(int beg, int end, int tube_id, int seg_len, std::vector<Segment*>& segs);
void get_front_segment(Segment* root, Segment*& front);
void get_back_segment(Segment* root, Segment*& back);
void get_front_segment2(Segment* root, Segment*& front);
void get_back_segment2(Segment* root, Segment*& back);

void ShowTubeSegmentsAndCollisions(Segment* head);

void find_circle(Segment* begin, Segment* current,
	std::stack<Segment*>& stack,
	std::vector<std::vector<Segment*>>& paths);

void find_circle(std::vector<Segment*> & heads, Segment* begin, std::vector<std::vector<Segment*>>& paths);
void gather_all_segments_in_circulation(std::vector<Segment*>& heads, Segment* seg, std::set<Segment*>& asic);

void dfs_collect_graph_segments(Segment* seg, std::vector<Segment*>& graph_segs);
void find_all_colli_segs(Segment* seg, std::vector<Segment*>& tmp_segs, std::vector<int>& a, std::vector<int>& b);