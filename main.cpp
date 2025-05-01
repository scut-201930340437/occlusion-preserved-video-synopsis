#include <iostream>
#include "graph.h"
#include "defines.h"
#include "loadata.h"
#include "visualize.h"
#include "mcmc.h"
#include "fast_forward.h"
#include "video_stitch.h"
#include "energy.h"
#include "frame_seq_fusion.h"
#include "octree_tube_rearrange.h"
#include "pccva.h"
#include<windows.h>
#include<string>
#include <dshow.h>
#include <shlobj.h>
#include<opencv2/opencv.hpp>
#include <mmsystem.h>
#include <stdio.h>
#include<tchar.h>
#pragma comment(lib, "winmm.lib")
#pragma comment(lib,"Strmiids.lib")
#pragma comment(lib, "strmiids")
#define _CRT_SECURE_NO_WARNINGS
using namespace std;

void init_graph(Context& context, double& min_s, double& max_s)
{
	/* The first step: build a graph */
	/*context.segCollisionMat = create_segment_graph2(context.filepath + "/tube_beg_end.txt",
		context.heads,
		context.all_segs,
		context.all_edges,
		context.seg_len,
		context.srclen,
		context.algo_type,
		context.non_colliding_divide_component,
		context.orig_colli);*/

	context.tube_num = int(context.all_segs.size());

	read_activity(context.filepath + "/activity", context.activity, context.tube_num);
	read_boundingbox(context.filepath + "/boundingbox", context.bbox, context.tube_num, context.scale, context.video_height, context.video_width);

	for (size_t i = 0; i < context.heads.size(); i++)
	{
		Segment* head = context.heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			std::cout << "(" << next->tube_id_ << ":" << next->a_ << "," << next->b_;

			for (size_t j = 0; j < next->colli_segs_.size(); j++)
			{
				std::cout << "," << next->colli_segs_[j]->tube_id_;
			}
			std::cout << ")" << std::endl;
			next = next->next_;
		}
		std::cout << std::endl;
	}


	//std::cout << context.collisionMat << std::endl;

	context.tube_num = int(context.all_segs.size());

	//计算建图时间
	/*end = clock();
	std::ofstream ofilestime0(context.filepath + "/step0time.txt");
	ofilestime0 << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	start = clock();*/


	/*context.collisionVec.resize(context.tube_num, 0);
	for (int i = 0; i < context.collisionMat.rows; i++)
	{
		for (int j = 0; j < context.collisionMat.cols; j++)
		{
			context.collisionVec[i] += context.collisionMat.at<uchar>(i, j);
		}

		std::cout << context.collisionVec[i] << std::endl;
	}

	context.circulated_tube.resize(context.tube_num, 0);*/



	/* Find circulation */
	/*for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			std::cout << i << "," << j;
			reset_traverse(context.all_segs);
			std::vector<std::pair<int, int>> chain;
			dfs_jump_chain(context.all_segs[i][j], chain);
			bool ret = check_circulation(chain);
			if (ret)
			{
				std::cout << "    yes" << std::endl;
				context.circulated_tube[i] = 1;
			}
			else
			{
				std::cout << "    no" << std::endl;
			}
		}
	}*/

	/* Find roots */
	reset_traverse(context.heads);
	context.roots.clear();

	//std::vector<std::pair<Segment*, Segment*>> influ_pairs;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::cout << i << std::endl;
		if (!context.all_segs[i][0]->is_traversed)
		{
			context.roots.push_back(context.all_segs[i][0]);
			//bfs_traverseMST_and_show(context, context.all_segs[i][0], true, influ_pairs);
			dfs_traverse(context.all_segs[i][0]);
		}
	}

	/* Change root to front */
	for (size_t i = 0; i < context.roots.size(); i++)
	{
		Segment* front = NULL;
		reset_traverse(context.heads);
		get_front_segment2(context.roots[i], front);
		Segment* back = NULL;
		reset_traverse(context.heads);
		get_back_segment2(context.roots[i], back);

		double graph_length = back->b_ - front->a_ + 1;

		std::cout << "Length:" << front->a_ << "," << back->b_ << "," << back->b_ - front->a_ + 1 << std::endl;

		context.roots[i] = front;

		if (back->b_ - front->a_ + 1 > context.synoplen)  //时间上的缩放，全部放入摘要长度中
		{
			double scale = double(context.synoplen - 100) / graph_length;
			std::vector<Segment*> graph_segs;
			reset_traverse(context.heads);
			dfs_collect_graph_segments(front, graph_segs);
			for (size_t kk = 0; kk < graph_segs.size(); kk++)
			{
				int oaa = graph_segs[kk]->aa_;
				int obb = graph_segs[kk]->bb_;

				graph_segs[kk]->bb_ = graph_segs[kk]->aa_ + int(graph_segs[kk]->len_ * scale);
				reset_traverse(context.heads);
				dfs_change_aabb(graph_segs[kk], oaa, obb);
				check_consistency(context);
			}
		}
	}



	/* Initialize segs' positions */
	reset_traverse(context.heads);
	for (size_t i = 0; i < context.roots.size(); i++)
	{
		/*double tmp = double(context.synoplen) / double(context.roots.size());
		tmp *= i;*/
		double tmp = 0;

		int oaa = context.roots[i]->aa_;
		int obb = context.roots[i]->bb_;
		int le = obb - oaa;

		context.roots[i]->aa_ = int(tmp);
		context.roots[i]->bb_ = context.roots[i]->aa_ + le;

		dfs_change_aabb(context.roots[i], oaa, obb);
		check_consistency(context);
	}

	if (!allin_range(context))
	{
		std::cout << "Init not complete: not all objects in range! ---------------------------------" << std::endl;
	}

	//可视化第一阶段
	//std::vector<std::vector<Segment*>> BCopy;
	//for (size_t i = 0; i < context.all_segs.size(); i++)
	//{
	//	std::vector<Segment*> tmp;
	//	for (size_t j = 0; j < context.all_segs[i].size(); j++)
	//	{
	//		tmp.push_back(new Segment(*context.all_segs[i][j]));
	//	}
	//	BCopy.push_back(tmp);
	//}
	//visualize1(context, BCopy);

	//int tube_num = int(context.all_segs.size());
	//std::vector<std::pair<int, int>> revert_pairs1;
	////std::vector<int> revert_degree;
	//for (int i = 0; i < tube_num - 1; i++)
	//{
	//	for (int j = i + 1; j < tube_num; j++)
	//	{
	//		/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
	//		{
	//			continue;
	//		}*/

	//		Segment* s1 = context.all_segs[i][0];
	//		Segment* s2 = context.all_segs[j][0];
	//		int ds1 = s1->a_ - s2->a_;
	//		int ds2 = s1->aa_ - s2->aa_;
	//		if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
	//		{
	//			revert_pairs1.push_back(std::pair<int, int>(i, j));
	//			//revert_degree.push_back(abs(ds2));
	//		}
	//	}
	//}

	//double enout1; // energy outside
	//double encolli1;
	//double enscale1;
	//std::vector<double> tube_energy1;
	//double enself1;
	//double nenergy1 = compute_energy(context, tube_energy1, enout1, encolli1, enscale1, enself1);

	//std::cout<<"CDN:" <<revert_pairs1.size() << " Collision:" << encolli1 << std::endl;
	//std::ofstream ofiles1(context.filepath + "/step1.txt");
	//ofiles1 << revert_pairs1.size()<<" "<<encolli1 << "\n";


	/*end = clock();
	std::ofstream ofilestime1(context.filepath + "/step1time.txt");
	ofilestime1 << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	start = clock();*/




	/* Run mcmc_chrono, ensuring chornological order */
	std::vector<std::vector<Segment*>> bestCopy;
	//mcmc_chrono(context, bestCopy, min_s, max_s); -----------------------

	for (size_t i = 0; i < bestCopy.size(); i++)
	{
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			context.all_segs[i][j]->aa_ = bestCopy[i][j]->aa_;
			context.all_segs[i][j]->bb_ = bestCopy[i][j]->bb_;
			delete bestCopy[i][j];
		}
	}

	bestCopy.clear();

	if (!chronological_order_preserved(context, context.chrono_slack))
	{
		std::cout << "Init wrong: not all objects chronological order preserved ! ---------------------------------" << std::endl;
	}
	else
	{
		std::cout << "Init correct" << std::endl;
	}

	//可视化第二阶段
	//std::vector<std::vector<Segment*>> CCopy;
	//for (size_t i = 0; i < context.all_segs.size(); i++)
	//{
	//	std::vector<Segment*> tmp;
	//	for (size_t j = 0; j < context.all_segs[i].size(); j++)
	//	{
	//		tmp.push_back(new Segment(*context.all_segs[i][j]));
	//	}
	//	CCopy.push_back(tmp);
	//}
	//visualize2(context, CCopy);

	////int tube_num2 = int(context.all_segs.size());
	//std::vector<std::pair<int, int>> revert_pairs2;
	////std::vector<int> revert_degree;
	//for (int i = 0; i < tube_num - 1; i++)
	//{
	//	for (int j = i + 1; j < tube_num; j++)
	//	{
	//		/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
	//		{
	//			continue;
	//		}*/

	//		Segment* s1 = context.all_segs[i][0];
	//		Segment* s2 = context.all_segs[j][0];
	//		int ds1 = s1->a_ - s2->a_;
	//		int ds2 = s1->aa_ - s2->aa_;
	//		if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
	//		{
	//			revert_pairs2.push_back(std::pair<int, int>(i, j));
	//			//revert_degree.push_back(abs(ds2));
	//		}
	//	}
	//}

	//double enout2; // energy outside
	//double encolli2;
	//double enscale2;
	//std::vector<double> tube_energy2;
	//double enself2;
	//double nenergy2 = compute_energy(context, tube_energy2, enout2, encolli2, enscale2, enself2);

	//std::ofstream ofiles2(context.filepath + "/step2.txt");
	//ofiles2 << revert_pairs2.size() << " " << encolli2 << "\n";
	/*end = clock();
	std::ofstream ofilestime2(context.filepath + "/step2time.txt");
	ofilestime2 << (double)(end - start) / CLOCKS_PER_SEC << "\n";*/


}

bool init_graph(Context & context,double& min_s, double& max_s, clock_t& start, clock_t& end)
{
	/* The first step: build a graph */
	context.segCollisionMat = create_segment_graph2(context.filepath+"/tube_beg_end.txt",
		context.heads,
		context.all_segs,
		context.all_edges,
		context.seg_len,
		context.srclen,
		context.algo_type,
		context.non_colliding_divide_component,
		context.orig_colli);

	context.tube_num = int(context.all_segs.size());

	read_activity(context.filepath + "/activity", context.activity, context.tube_num);
	read_boundingbox(context.filepath + "/boundingbox", context.bbox, context.tube_num, context.scale, context.video_height, context.video_width);

	std::cout << context.heads.size() << std::endl;
	for (size_t i = 0; i < context.heads.size(); i++)
	{
		Segment* head = context.heads[i];
		Segment* next = head->next_;

		while (next != NULL)
		{
			//std::cout << "(" << next->tube_id_ << ":" << next->a_ << "," << next->b_;

			for (size_t j = 0; j < next->colli_segs_.size(); j++)
			{
				//std::cout << ","<< next->colli_segs_[j]->tube_id_;
			}
			//std::cout << ")" << std::endl;
			next = next->next_;
		}
		//std::cout << std::endl;
	}
	//std::cout << context.collisionMat << std::endl;
	
	context.tube_num = int(context.all_segs.size());
	

	//计算建图时间
	end = clock();
	std::ofstream ofilestime0(context.filepath + "/step0time.txt");
	ofilestime0 << (double)(end - start) / CLOCKS_PER_SEC  << "\n";
	std::cout << "step 0 time: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
	start = clock();
	

	/*context.collisionVec.resize(context.tube_num, 0);
	for (int i = 0; i < context.collisionMat.rows; i++)
	{
		for (int j = 0; j < context.collisionMat.cols; j++)
		{
			context.collisionVec[i] += context.collisionMat.at<uchar>(i, j);
		}

		std::cout << context.collisionVec[i] << std::endl;
	}

	context.circulated_tube.resize(context.tube_num, 0);*/
	


	/* Find circulation */
	/*for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			std::cout << i << "," << j;
			reset_traverse(context.all_segs);
			std::vector<std::pair<int, int>> chain;
			dfs_jump_chain(context.all_segs[i][j], chain);
			bool ret = check_circulation(chain);
			if (ret)
			{
				std::cout << "    yes" << std::endl;
				context.circulated_tube[i] = 1;
			}
			else
			{
				std::cout << "    no" << std::endl;
			}
		}
	}*/

	/* Find roots */
	reset_traverse(context.heads);
	context.roots.clear();

	//std::vector<std::pair<Segment*, Segment*>> influ_pairs;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::cout << i << std::endl;
		if (!context.all_segs[i][0]->is_traversed)
		{
			context.roots.push_back(context.all_segs[i][0]);
			//bfs_traverseMST_and_show(context, context.all_segs[i][0], true, influ_pairs);
			dfs_traverse(context.all_segs[i][0]);
		}
	}

	/* Change root to front */
	for (size_t i = 0; i < context.roots.size(); i++)
	{
		Segment* front = NULL;
		reset_traverse(context.heads);
		get_front_segment2(context.roots[i], front);
		Segment* back = NULL;
		reset_traverse(context.heads);
		get_back_segment2(context.roots[i], back);

		int graph_length = back->b_ - front->a_ + 1;

		/*std::vector<Segment*> segs;
		reset_traverse(context.heads);
		dfs_collect_graph_segments(front, segs);
		int min_a = 10000000;
		int max_b = -1;
		for (int kk = 0; kk<int(segs.size()); ++kk)
		{
			min_a = min(min_a, segs[kk]->a_);
			max_b = max(max_b, segs[kk]->b_);
		}

		int graph_length = max_b - min_a + 1;*/

		//std::cout << "Length:" << front->a_ << "," << back->b_ << "," << graph_length << std::endl;

		context.roots[i] = front;

		if (graph_length > context.synoplen)  //时间上的缩放，全部放入摘要长度中
		{
			//double scale = double(context.synoplen - 100) / graph_length;
			double scale = double(context.synoplen - 150) / graph_length; // 调整synop
			//double scale = 0.01;
			std::vector<Segment*> graph_segs;
			reset_traverse(context.heads);
			dfs_collect_graph_segments(front, graph_segs);
			for (size_t kk = 0; kk < graph_segs.size(); kk++)
			{
				int oaa = graph_segs[kk]->aa_;
				int obb = graph_segs[kk]->bb_;
				
				graph_segs[kk]->bb_ = graph_segs[kk]->aa_ + int(graph_segs[kk]->len_ * scale);
				reset_traverse(context.heads);
				dfs_change_aabb(graph_segs[kk], oaa, obb);
				check_consistency(context);
			}
		}
	}

	

	/* Initialize segs' positions */
	reset_traverse(context.heads);
	for (size_t i = 0; i < context.roots.size(); i++)
	{
		/*double tmp = double(context.synoplen) / double(context.roots.size());
		tmp *= i;*/
		//double tmp = 0;

		int oaa = context.roots[i]->aa_;
		int obb = context.roots[i]->bb_;
		int le = obb - oaa;
		
		context.roots[i]->aa_ = 0;
		context.roots[i]->bb_ = context.roots[i]->aa_ + le;

		dfs_change_aabb(context.roots[i], oaa, obb);
		check_consistency(context);
	}

	if (!allin_range(context))
	{
		std::cout << "Init not complete: not all objects in range! ---------------------------------" << std::endl;
	}

	//可视化第一阶段
	//std::vector<std::vector<Segment*>> BCopy;
	//for (size_t i = 0; i < context.all_segs.size(); i++)
	//{
	//	std::vector<Segment*> tmp;
	//	for (size_t j = 0; j < context.all_segs[i].size(); j++)
	//	{
	//		tmp.push_back(new Segment(*context.all_segs[i][j]));
	//	}
	//	BCopy.push_back(tmp);
	//}
	//visualize1(context, BCopy);

	//int tube_num = int(context.all_segs.size());
	//std::vector<std::pair<int, int>> revert_pairs1;
	////std::vector<int> revert_degree;
	//for (int i = 0; i < tube_num - 1; i++)
	//{
	//	for (int j = i + 1; j < tube_num; j++)
	//	{
	//		/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
	//		{
	//			continue;
	//		}*/

	//		Segment* s1 = context.all_segs[i][0];
	//		Segment* s2 = context.all_segs[j][0];
	//		int ds1 = s1->a_ - s2->a_;
	//		int ds2 = s1->aa_ - s2->aa_;
	//		if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
	//		{
	//			revert_pairs1.push_back(std::pair<int, int>(i, j));
	//			//revert_degree.push_back(abs(ds2));
	//		}
	//	}
	//}

	//double enout1; // energy outside
	//double encolli1;
	//double enscale1;
	//std::vector<double> tube_energy1;
	//double enself1;
	//double nenergy1 = compute_energy(context, tube_energy1, enout1, encolli1, enscale1, enself1);

	//std::cout<<"CDN:" <<revert_pairs1.size() << " Collision:" << encolli1 << std::endl;
	//std::ofstream ofiles1(context.filepath + "/step1.txt");
	//ofiles1 << revert_pairs1.size()<<" "<<encolli1 << "\n";
	

	end = clock();
	std::ofstream ofilestime1(context.filepath + "/step1time.txt");
	ofilestime1 << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	

	/*ofstream before_53(context.filepath + "/before_53.txt");
	std::cout << "************ 53 tube **************" << std::endl;
	before_53 << "************ 53 tube **************\n";

	Segment* head = context.heads[52];
	head = head->next_;
	while (head != NULL)
	{
		std::cout << head->tube_id_ << " " << head->a_ << " " << head->b_ << " " << head->aa_ << " " << head->bb_ << std::endl;
		before_53 << head->tube_id_ << " " << head->a_ << " " << head->b_ << " " << head->aa_ << " " << head->bb_ << "\n";

		head = head->next_;
	}

	std::cout << "************ 53 tube **************" << std::endl;
	before_53 << "************ 53 tube **************\n";*/

	// 计算 seg 平均长度
	context.seg_avg_len = 0;
	int seg_num = 0;
	for (size_t i = 0; i < context.heads.size(); i++)
	{
		Segment* head = context.heads[i];
		head = head->next_;
		while (head != NULL)
		{
			seg_num++;
			seg_num += head->b_ - head->a_;
			head = head->next_;
		}
	}
	context.seg_avg_len /= seg_num;

	// 计算各tube的速度
	//context.tube_speed.resize(context.bbox.size());
	for (size_t p1 = 0; p1 < context.all_segs.size(); ++p1)
	{
		int tube_beg = context.bbox[p1][0][0];
		for (size_t p2 = 0; p2 < context.all_segs[p1].size(); ++p2) {
			Segment* seg = context.all_segs[p1][p2];
			if (seg->b_ - seg->a_ + 1 >= context.seg_len_lbound)
			{
				double speed = 0.0;
				for (int i = tube_beg; i <= seg->b_ - 10; ++i)
				{
					int x1 = context.bbox[seg->tube_id_][i - tube_beg][1], y1 = context.bbox[seg->tube_id_][i - tube_beg][2];
					int x2 = context.bbox[seg->tube_id_][i + 10 - tube_beg][1], y2 = context.bbox[seg->tube_id_][i + 10 - tube_beg][2];
					double dis = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) + 1e-15;
					speed += dis / 10.0;
				}
				speed /= double(seg->b_ - 10 - tube_beg + 1);
				context.tube_speed[seg->tube_id_][seg->a_] = speed;
			}
		}
	}

	/* Run mcmc_chrono, ensuring chornological order */
	start = clock();

	std::vector<std::vector<Segment*>> bestCopy;
	int iter_chrono = mcmc_chrono(context, bestCopy, min_s, max_s);

	for (size_t i = 0; i < bestCopy.size(); i++)
	{
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			context.all_segs[i][j]->aa_ = bestCopy[i][j]->aa_;
			context.all_segs[i][j]->bb_ = bestCopy[i][j]->bb_;
			delete bestCopy[i][j];
		}
	}

	bestCopy.clear();

	bool can_preserve_order = true;
	
	if (!chronological_order_preserved(context, context.chrono_slack))
	{
		std::cout << "Init wrong: not all objects chronological order preserved ! ---------------------------------" << std::endl;
		can_preserve_order = false;
	}
	else
	{
		std::cout << "Init correct" << std::endl;
	}

	/*可视化第二阶段*/
	//std::vector<std::vector<Segment*>> CCopy;
	//for (size_t i = 0; i < context.all_segs.size(); i++)
	//{
	//	std::vector<Segment*> tmp;
	//	for (size_t j = 0; j < context.all_segs[i].size(); j++)
	//	{
	//		tmp.push_back(new Segment(*context.all_segs[i][j]));
	//	}
	//	CCopy.push_back(tmp);
	//}
	//visualize2(context, CCopy);

	////int tube_num2 = int(context.all_segs.size());
	//std::vector<std::pair<int, int>> revert_pairs2;
	////std::vector<int> revert_degree;
	//for (int i = 0; i < tube_num - 1; i++)
	//{
	//	for (int j = i + 1; j < tube_num; j++)
	//	{
	//		/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
	//		{
	//			continue;
	//		}*/

	//		Segment* s1 = context.all_segs[i][0];
	//		Segment* s2 = context.all_segs[j][0];
	//		int ds1 = s1->a_ - s2->a_;
	//		int ds2 = s1->aa_ - s2->aa_;
	//		if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
	//		{
	//			revert_pairs2.push_back(std::pair<int, int>(i, j));
	//			//revert_degree.push_back(abs(ds2));
	//		}
	//	}
	//}

	//double enout2; // energy outside
	//double encolli2;
	//double enscale2;
	//std::vector<double> tube_energy2;
	//double enself2;
	//double nenergy2 = compute_energy(context, tube_energy2, enout2, encolli2, enscale2, enself2);

	//std::ofstream ofiles2(context.filepath + "/step2.txt");
	//ofiles2 << revert_pairs2.size() << " " << encolli2 << "\n";
	end = clock();
	std::ofstream ofilestime2(context.filepath + "/step2time.txt");
	ofilestime2 << (double)(end - start) / CLOCKS_PER_SEC <<" "<<iter_chrono<< "\n"; //------------------------------
	


	//start = clock();
	//// Run beg chrono
	//std::vector<std::vector<Segment*>> bestCopy2;
	//iter_chrono = mcmc_beg_chrono(context, bestCopy2, min_s, max_s);

	//for (size_t i = 0; i < bestCopy2.size(); i++)
	//{
	//	for (size_t j = 0; j < bestCopy2[i].size(); j++)
	//	{
	//		context.all_segs[i][j]->aa_ = bestCopy2[i][j]->aa_;
	//		context.all_segs[i][j]->bb_ = bestCopy2[i][j]->bb_;
	//		delete bestCopy2[i][j];
	//	}
	//}

	/*bestCopy2.clear();

	if (!chronological_order_preserved(context, context.chrono_slack))
	{
		std::cout << "Init wrong: not all objects chronological order preserved ! ---------------------------------" << std::endl;
	}
	else
	{
		std::cout << "Init correct" << std::endl;
	}

	end = clock();
	std::ofstream ofilestime3(context.filepath + "/step3time.txt");
	ofilestime3 << (double)(end - start) / CLOCKS_PER_SEC << " " << iter_chrono << "\n";*/

	return can_preserve_order;
}


//int main_gogo3()
//{
//	std::string filepath = "C:/Users/nieyongwei/Desktop/speedsynop/gogo-3";
//
//	std::ofstream ofile(filepath + "/comare_with_algo1.csv", std::ios_base::trunc |std::ios_base::out);
//	
//	int synlen[11] = { 3000, 2800, 2600,2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000 };
//	for (int i = 0; i < 1; i++)
//	{
//		for (int j = 0; j < 1; j++)
//		{
//			Context context;
//			context.filepath = filepath;
//			context.scale = 0.5;    //gogo
//			context.seg_len = 400;
//			context.fps = 25;
//
//			context.chrono_slack = 200;
//			/* The following two hyper-params may need to be adjusted simultaneously */
//			context.algo_type = 0; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
//			context.resizing_component = false;	// use the resizing strategy or not
//			context.non_colliding_divide_component = false;
//
//			context.mcmc2_iteration_num = 500000;
//			context.sigma_decay = 1;
//
//			context.video_height = int(1024 * context.scale);
//			context.video_width = int(576 * context.scale);
//
//			context.synoplen = synlen[i];
//			context.sigma = synlen[i];
//
//			srand(time(NULL));
//
//			init_graph(context);
//			std::vector<std::vector<Segment*>> bestCopy;
//			int iter_best = context.mcmc2_iteration_num;
//			double best_energy = mcmc2(context, bestCopy,iter_best);
//
//			for (size_t i = 0; i < bestCopy.size(); i++)
//			{
//				for (size_t j = 0; j < bestCopy[i].size(); j++)
//				{
//					delete bestCopy[i][j];
//				}
//			}
//			bestCopy.clear();
//
//			release_context(context);
//			ofile << best_energy << ",";
//		}
//		ofile << std::endl;
//	}
//
//	return 0;
//}

//int main_gogo_full2()
//{
//	std::string filepath = "C:/Users/nieyongwei/Desktop/speedsynop/gogo-full2";
//
//	std::ofstream ofile(filepath + "/results.csv", std::ios_base::trunc | std::ios_base::out);
//
//	int synlen[4] = { 5000, 4500, 4000,3500};
//	for (int i = 0; i < 4; i++)
//	{
//		time_t seed = time(NULL);
//
//		srand(seed);
//		for (int j = 0; j < 5; j++)
//		{
//			Context context;
//			context.filepath = filepath;
//			context.scale = 0.25;    //gogo
//			context.seg_len = 400;
//			context.fps = 25;
//
//			context.chrono_slack = 200;
//			/* The following two hyper-params may need to be adjusted simultaneously */
//			context.algo_type = 0; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
//			context.resizing_component = false;	// use the resizing strategy or not
//			context.non_colliding_divide_component = false;
//
//			context.mcmc2_iteration_num = 500000;
//			context.sigma_decay = 1;
//
//			context.video_height = int(1280 * context.scale);
//			context.video_width = int(720 * context.scale);
//
//			context.synoplen = synlen[i];
//			context.sigma = synlen[i];
//
//			srand(time(NULL));
//
//			init_graph(context);
//			std::vector<std::vector<Segment*>> bestCopy;
//			int iter_best = context.mcmc2_iteration_num;
//			double best_energy = mcmc2(context, bestCopy,iter_best);
//
//			for (size_t i = 0; i < bestCopy.size(); i++)
//			{
//				for (size_t j = 0; j < bestCopy[i].size(); j++)
//				{
//					delete bestCopy[i][j];
//				}
//			}
//			bestCopy.clear();
//
//			release_context(context);
//			ofile << best_energy << ",";
//		}
//		ofile << std::endl;
//
//		srand(seed);
//		for (int j = 0; j < 5; j++)
//		{
//			Context context;
//			context.filepath = filepath;
//			context.scale = 0.25;    //gogo
//			context.seg_len = 400;
//			context.fps = 25;
//
//			context.chrono_slack = 200;
//			/* The following two hyper-params may need to be adjusted simultaneously */
//			context.algo_type = 1; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
//			context.resizing_component = false;	// use the resizing strategy or not
//			context.non_colliding_divide_component = false;
//
//			context.mcmc2_iteration_num = 500000;
//			context.sigma_decay = 1;
//
//			context.video_height = int(1280 * context.scale);
//			context.video_width = int(720 * context.scale);
//
//			context.synoplen = synlen[i];
//			context.sigma = synlen[i];
//
//			srand(time(NULL));
//
//			init_graph(context);
//			std::vector<std::vector<Segment*>> bestCopy;
//			int iter_best = context.mcmc2_iteration_num;
//			double best_energy = mcmc2(context, bestCopy, iter_best);
//
//			for (size_t i = 0; i < bestCopy.size(); i++)
//			{
//				for (size_t j = 0; j < bestCopy[i].size(); j++)
//				{
//					delete bestCopy[i][j];
//				}
//			}
//			bestCopy.clear();
//
//			release_context(context);
//			ofile << best_energy << ",";
//		}
//		ofile << std::endl;
//	}
//
//	return 0;
//}




int main_st2()
{
	std::string filepath = "E:/VideoSyn_scopic/video6";

	std::ofstream ofile(filepath + "/LiG_algo2_biaozhuncha.csv", std::ios_base::trunc | std::ios_base::out);

	const int len = 1;
	int synlen[len] = {380};
	//int synlen[3] = { 500,400,300 };
	//int slack[9] = { 400,350,300,250,200,150,100,50,0 };
	/*double mins[9] = {0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8};
	double maxs[9] = {50,25,16.7,12.5,10,5,2.5,1.7,1.25};*/

	for (int i = 0; i < len; i++)
	{
		ofile << synlen[i] << ",";
		for (int j = 0; j <3; j++)
		{
			Context context;
			context.filepath = filepath;
			context.scale = 1;
			context.seg_len = 400;
			context.fps = 30;

			context.chrono_slack = 120;
			//context.chrono_slack = slack[j];
			/* The following two hyper-params may need to be adjusted simultaneously */
			context.algo_type = 2; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
			context.resizing_component = false;	// use the resizing strategy or not
			context.remainspeed_component = false;
			context.non_colliding_divide_component = true;

			context.mcmc2_iteration_num = 500000;
			context.sigma_decay = 1;

			context.video_height = int(960 * context.scale);
			context.video_width = int(540 * context.scale);

			context.synoplen = synlen[i];
			context.sigma = synlen[i];
			/*double min_s = mins[j];
			double max_s = maxs[j];*/

			double min_s = 0.2;
			double max_s = 5;

			srand(time(NULL));

			init_graph(context,min_s,max_s);
			std::vector<std::vector<Segment*>> bestCopy;

			int iter_best = context.mcmc2_iteration_num;
			//double best_energy = mcmc2(context, iter_best,min_s,max_s); ------------------
			double best_energy = 0;
			/*if (iter_best < context.mcmc2_iteration_num) {
				std::vector<double> tube_energy2;
				double enself2;
				double best_enco = compute_energy_collision5(context, tube_energy2, enself2);
				int colli = best_enco;
				if (colli > 0) {
					iter_best = context.mcmc2_iteration_num - iter_best - 1;
					best_energy = mcmc3(context, bestCopy, iter_best, min_s, max_s);
				}
			}*/

			for (size_t i = 0; i < bestCopy.size(); i++)
			{
				for (size_t j = 0; j < bestCopy[i].size(); j++)
				{
					delete bestCopy[i][j];
				}
			}
			bestCopy.clear();

			release_context(context);
			int t_best_energy = best_energy;
			ofile << t_best_energy-1 << ",";
		}
		ofile << std::endl;
	}

	return 0;
}

int main_junxun()
{
	srand(time(NULL));
	Context context;
	//context.filepath = "C:/Users/nieyongwei/Desktop/speedsynop/st2";
	context.filepath = "E:/VideoSyn_scopic/video3";
	//context.filepath = "C:/Users/nieyongwei/Desktop/speedsynop/greedy";
	

	context.scale = 1;    //gogo
	context.seg_len = 400;
	context.fps = 30;

	context.synoplen = 250;
	context.chrono_slack = 40;
	/* The following two hyper-params may need to be adjusted simultaneously */
	context.algo_type = 0; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
	context.resizing_component = false;	// use the resizing strategy or not
	context.remainspeed_component = false;
	context.non_colliding_divide_component = true; //与方法2没关系

	context.mcmc2_iteration_num = 500000;
	//context.mcmc2_iteration_num = 1000;
	context.sigma = context.synoplen;
	context.sigma_decay = 1;


	context.video_height = int(1920 * context.scale);
	context.video_width = int(1080 * context.scale);
	/*std::string filepath = "C:/Users/nieyongwei/Desktop/speedsynop/st2";
	stich_two_videos(filepath + "/synopsis-fast_forward-1000.mp4",
		filepath + "/synopsis.mp4",
		filepath + "/stitch-1000-2.mp4");
	return 0;*/

	double min_s = 0.1;
	double max_s = 10000;

	init_graph(context, min_s, max_s);

	/*fast_forward(context.filepath + "/src.mp4", context.filepath + "/synopsis-fast_forward.mp4", 1000, context.scale);
	release_context(context);
	return 0;*/

	std::vector<std::vector<Segment*>> bestCopy;

	
	int iter_best = context.mcmc2_iteration_num;
	//mcmc2(context, iter_best, min_s, max_s); //-------------------------------


	std::cout << "Visualization" << std::endl;
	//visualize(context);
	std::cout << "Visualization end!" << std::endl;

	std::ofstream ofile(context.filepath + "/bestCopy.txt");
	for (size_t i = 0; i < bestCopy.size(); i++)
	{
		ofile << i << ":\n";
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			ofile << j << ":" << bestCopy[i][j]->a_ << "," << bestCopy[i][j]->b_ << ","
				<< bestCopy[i][j]->aa_ << "," << bestCopy[i][j]->bb_ << '\n';
		}
	}

	for (size_t i = 0; i < bestCopy.size(); i++)
	{
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			delete bestCopy[i][j];
		}
	}
	bestCopy.clear();

	release_context(context);

	return 0;
}

int main_video1gogo()
{
	srand(time(NULL));
	Context context;
	//context.filepath = "C:/Users/nieyongwei/Desktop/speedsynop/st2";
	context.filepath = "E:/VideoSyn_scopic/video1";
	

	context.scale = 1;    //gogo
	context.seg_len = 200;
	context.fps = 25;

	context.synoplen = 3000;
	context.chrono_slack = 500;
	/* The following two hyper-params may need to be adjusted simultaneously */
	context.algo_type = 0; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
	context.resizing_component = false;	// use the resizing strategy or not
	context.remainspeed_component = false;
	context.non_colliding_divide_component = true; //与方法2没关系

	context.mcmc2_iteration_num = 500000;
	//context.mcmc2_iteration_num = 1000;
	context.sigma = context.synoplen;
	context.sigma_decay = 1;


	context.video_height = int(1280 * context.scale);
	context.video_width = int(720 * context.scale);

	double min_s = 0.08;
	double max_s = 10000;

	init_graph(context,min_s,max_s);

	/*fast_forward(context.filepath + "/src.mp4", context.filepath + "/synopsis-fast_forward.mp4", 1000, context.scale);
	release_context(context);
	return 0;*/

	std::vector<std::vector<Segment*>> bestCopy;

	
	int iter_best = context.mcmc2_iteration_num;
	//mcmc2(context, iter_best, min_s, max_s); // ---------------


	std::cout << "Visualization" << std::endl;
	//visualize(context);
	std::cout << "Visualization end!" << std::endl;

	std::ofstream ofile(context.filepath + "/bestCopy.txt");
	for (size_t i = 0; i < bestCopy.size(); i++)
	{
		ofile << i << ":\n";
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			ofile << j << ":" << bestCopy[i][j]->a_ << "," << bestCopy[i][j]->b_ << ","
				<< bestCopy[i][j]->aa_ << "," << bestCopy[i][j]->bb_ << '\n';
		}
	}

	for (size_t i = 0; i < bestCopy.size(); i++)
	{
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			delete bestCopy[i][j];
		}
	}
	bestCopy.clear();

	release_context(context);

	return 0;
}

int main()
{
	clock_t start, end;
	start = clock();
	srand(time(NULL));
	Context context;
	//context.filepath = "C:/Users/nieyongwei/Desktop/speedsynop/st2";

	//
	//context.filepath = "E:/VideoSyn_scopic/video2";
	//context.filepath = "D:/VideoSyn_Pro/synopsis_dataset/video-st";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/src_scut_secSeg1";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/src1/640x368";
	
	//context.filepath = path;

	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/ours/src1/1-34";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/ours/src_scut_sec/res_qual";
	context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/ours/src_scut_fir_1/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/ours/src_scut_fir_2/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/ours/src_scut_fir_2/644-693/675-693";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/ours/MOT17/02-SDP";

	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/size_speed/src_scut_sec/324-456/444-456";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/size_speed/MOT20/02";

	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/object_interaction-based/src_scut_sec/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/object_interaction-based/src_scut_fir_1/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/object_interaction-based/src_scut_fir_2/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/object_interaction-based/MOT20/05";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/object_interaction-based/dancetrack/0080";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/object_interaction-based/src_scut_fir_2/716-774";

	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/octree/src_scut_sec/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/octree/src_scut_fir_1/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/octree/src_scut_fir_2/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/octree/MOT20/05";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/octree/dancetrack/0080";

	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/pccva/src_scut_sec/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/pccva/src_scut_fir_1/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/pccva/src_scut_fir_2/res_qual";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/pccva/MOT20/03";
	//context.filepath = "D:/pytorchPro/Towards-Realtime-MOT/input/pccva/dancetrack/0080";


	context.scale = 1;    //gogo
	context.seg_len = 400;

	context.fps = 30; // -------

	//context.synoplen = int(double(83 - 72 + 1) * 12.17951807228916 + 0.5); // 摘要长度 ------
	//context.synoplen = int(double(250 - 224 + 1) * 5.304839764551995 + 0.5); // 摘要长度 ------
	//context.synoplen = 459 * 3 * 0.75;
	context.synoplen = 785;

	context.chrono_slack = 400; // 时序损失，允许两个物体先后顺序逆转的帧数
	context.seg_len_lbound = 15;
	/* The following two hyper-params may need to be adjusted simultaneously */
	context.algo_type = 0; // 0: our algorithm   1: full collision in collision region   2: full collisions over all segments
	context.resizing_component = false;	// use the resizing strategy or not
	context.remainspeed_component = false;
	context.non_colliding_divide_component = false; //与方法2没关系

	context.mcmc2_iteration_num = 3000;
	context.sigma = context.synoplen;
	context.sigma_decay = 1;

	double limit_time = 191; //----------

	// 视频宽高
	context.video_height = int(720 * context.scale); // pccva 需要width和height参数 --------------------------------------------
	context.video_width = int(1280 * context.scale);

	// 速度
	double speed = 0.3899; //----------
	//double speed = 0.1588;

	double min_s = 0.1;
	double max_s = 8.0;
	
	// 建立FOG->消除循环冲突->确保所有tube保留->保持时间顺序
	bool can_preserve_order = true; // 只用于方法2 // -----------------------------------
	//can_preserve_order = init_graph(context,min_s,max_s,start,end);   //----------------------------------

	/*fast_forward(context.filepath + "/src.mp4", context.filepath + "/synopsis-fast_forward.mp4", 1000, context.scale);
	release_context(context);
	return 0;*/

	//std::vector<std::vector<Segment*>> bestCopy;
	
	start = clock();
	int iter_best = context.mcmc2_iteration_num;
	//mcmc2(context, iter_best, min_s, max_s, can_preserve_order);  //-------------------------------------------

	//iter_best=(int)mcmc_chro_colli(context, iter_best, min_s, max_s,limit_time); // 同时优化
	

	// 其他方法
	std::map<int, int> tubes_group_id;
	//fusion(context, tubes_group_id, speed); // fusion 方法

	//octree_tube_rerrange(context, tubes_group_id, speed); // octree 方法

	//pccva(context, tubes_group_id, speed); // pccva 方法

	end = clock();
	std::ofstream ofilestime3(context.filepath + "/step3time.txt", ofstream::app);
	std::cout<<"step 3 time: "<< (double)(end - start) / CLOCKS_PER_SEC << std::endl;
	ofilestime3 << (double)(end - start) / CLOCKS_PER_SEC << "\n";


	// 可视化------------------------------- 直接从结果文件生成摘要视频
	context.tube_num = 44 - 1 + 1;
	read_boundingbox(context.filepath + "/boundingbox", context.bbox, context.tube_num, context.scale, context.video_height, context.video_width);
	std::map<int, std::vector<int>> groups;
	read_bestCopy(context.filepath, context.all_segs);

	read_occlu(context.filepath, context.orig_colli, context.occlu_graph);
	

	// object方法和pccva
	/*read_group(context.filepath, groups, context.tube_num);
	for (auto it : groups)
	{
		for (size_t i = 0; i < it.second.size(); ++i)
		{
			tubes_group_id[it.second[i]] = it.first;
		}
	}*/

	// octree方法
	/*for (int i = 1; i <= context.tube_num; ++i)
	{
		tubes_group_id[i] = i;
	}*/
	
	/* 加速 （如果需要的话）*/
	// octree 只需要调整tube的长度
	//for (int i = 0; i < int(context.bbox.size()); ++i)
	//{
	//	int beg = context.bbox[i][0][0];
	//	int tmp = int(context.bbox[i].size()) - 1;
	//	int end = context.bbox[i][tmp][0];
	//	// 计算加速后的新结束帧号
	//	int newEnd = beg + round(speed * (end - beg + 1)) - 1;
	//	std::vector<cv::Vec6d> bboxes_cur_tube;
	//	// 获取加速后，当前tube的各帧的bbox
	//	for (int kk = beg; kk <= newEnd; ++kk)
	//	{
	//		int k;
	//		if (beg == newEnd)
	//		{
	//			k = beg;
	//		}
	//		else
	//		{
	//			k = (int)(((double)(kk - beg) / (newEnd - beg)) * (end - beg) + beg);
	//		}
	//		// 更新帧号
	//		cv::Vec6d bbox = context.bbox[i][k - beg];
	//		bbox[0] = kk;
	//		// 记录原始帧号
	//		bbox[5] = k;
	//		bboxes_cur_tube.push_back(bbox);
	//	}
	//	context.bbox[i].assign(bboxes_cur_tube.begin(), bboxes_cur_tube.end());
	//}
	
	// object and pccva 还需要调整同一个group中tube的间隔
	//for (int group_id = 0; group_id < int(groups.size()); ++group_id)
	//{
	//	int tube_id = groups[group_id][0] - 1;
	//	int a_0 = context.bbox[tube_id][0][0];
	//	// 后续tube
	//	for (size_t j = 1; j < groups[group_id].size(); ++j)
	//	{
	//		int tube_id = groups[group_id][j] - 1;
	//		int a = context.bbox[tube_id][0][0];
	//		// 原本的tube的间隔也要加速
	//		int gap = int(double(a - a_0) * speed + 0.5);
	//		int aa = a_0 + gap;
	//		int tmp = int(context.bbox[tube_id].size()) - 1;
	//		int b = context.bbox[tube_id][tmp][0];
	//		int bb = b - (a - aa);
	//		std::vector<cv::Vec6d> bboxes_cur_tube;
	//		for (int kk = aa; kk <= bb; ++kk)
	//		{
	//			// 更新帧号
	//			cv::Vec6d bbox = context.bbox[tube_id][kk - aa];
	//			bbox[0] = kk;
	//			bbox[5] = kk + a - aa;
	//			bboxes_cur_tube.push_back(bbox);
	//		}
	//		context.bbox[tube_id].assign(bboxes_cur_tube.begin(), bboxes_cur_tube.end());
	//	}
	//}


	std::cout << "Visualization_step" << std::endl; // ---------------------------------------------
	visualize_for_ours(context, context.chrono_slack);
	//visualize_for_ours_step(context);
	//visualize(context, tubes_group_id, speed, context.chrono_slack);
	std::cout << "Visualization end!" << std::endl;

	std::cout << "compute CA." << endl; // only for compared methods //--------------------------------------------------------
	/*std::map<int, std::set<int>> colli_relation;
	double colli_area = 0.0;
	for (auto it : context.bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size()) - 1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				if (tubes_group_id[it.second[i].first] != tubes_group_id[it.second[j].first])
				{
					double frame_colli = cal_colli_area(it.second[i].second, it.second[j].second);
					colli_area += frame_colli;
					if (frame_colli > 0)
					{
						colli_relation[it.second[i].first].insert(it.second[j].first);
					}
				}
			}
		}
	}*/

	/*int colli_num = 0;
	for (auto it : colli_relation)
	{
		colli_num += int(it.second.size());
	}*/

	std::cout << "compute CDN." << endl;

	// 统计时序损失
	int chrono_cost = 0;
	for (int i = 0; i < int(context.all_segs.size()) - 1; i++)
	{
		if (int(context.all_segs[i].size()) > 0)
		{
			for (size_t j = i + 1; j < context.all_segs.size(); j++)
			{
				if (int(context.all_segs[j].size()) > 0)
				{
					Segment* s1 = context.all_segs[i][0];
					Segment* s2 = context.all_segs[j][0];
					int d1 = s1->a_ - s2->a_;
					int d2 = s1->aa_ - s2->aa_;
					if (d1 * d2 < 0 && abs(d2) > context.chrono_slack)
					{
						++chrono_cost;
					}
				}
			}
		}
	}


	std::ofstream ofile_metrics(context.filepath + "/metrics.txt", ofstream::app);
	std::ofstream ofile(context.filepath + "/bestCopy.txt");
	std::ofstream ofile_bestCopyVis(context.filepath + "/bestCopyVis.txt");
	ofile_metrics << "synoplen: " << context.synoplen << "\n";
	ofile_metrics << "chrono_slack: " << context.chrono_slack << "\n";
	ofile_metrics << "mcmc2_iteration_num: " << iter_best << "\n";
	//ofile_metrics << "Collision: " << colli_area << "\n";  //--------------------------------------------------
	ofile_metrics << "chrono_cost: " << chrono_cost << "\n";
	//ofile_metrics << "colli_num: " << colli_num << "\n";
	// ---------
	
	// 统计平均速度
	double avg_speed = 0.0;
	int tubes_num = 0;
	int segs_num = 0;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		ofile << i << ":\n";
		/*if (context.visualize_state[i])
		{
			ofile_bestCopyVis << i << ";\n";
		}*/
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			ofile << j << ":" << context.all_segs[i][j]->a_ << "," << context.all_segs[i][j]->b_ << ","
				<< context.all_segs[i][j]->aa_ << "," << context.all_segs[i][j]->bb_ << '\n';
			/*if (context.visualize_state[i])
			{
				ofile_bestCopyVis << j << ":" << context.all_segs[i][j]->a_ << "," << context.all_segs[i][j]->b_ << ","
					<< context.all_segs[i][j]->aa_ << "," << context.all_segs[i][j]->bb_ << '\n';
			}*/
			segs_num++;
		}
		// 统计速度
		int tmp = int(context.all_segs[i].size()) - 1;
		avg_speed += 1.0 * (context.all_segs[i][tmp]->bb_ - context.all_segs[i][0]->aa_ + 1) / (context.all_segs[i][tmp]->b_ - context.all_segs[i][0]->a_ + 1);
		tubes_num++;
	}
	ofile_metrics <<"speed sum: "<<avg_speed << " avg_speed: " << avg_speed / tubes_num << " segs num: " << tubes_num << "\n";
	ofile_metrics << "segs_num: " << segs_num << "\n";
	ofile_metrics << "\n";


	/*for (size_t i = 0; i < bestCopy.size(); i++)
	{
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			ofile << bestCopy[i][j]->aa_ << " " << bestCopy[i][j]->bb_ << '\n';
		}
	}*/

	/*for (size_t i = 0; i < bestCopy.size(); i++)
	{
		for (size_t j = 0; j < bestCopy[i].size(); j++)
		{
			delete bestCopy[i][j];
		}
	}
	bestCopy.clear();*/

	release_context(context);

	return 0;
}