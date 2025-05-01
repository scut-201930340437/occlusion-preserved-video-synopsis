#include <iostream>
#include <fstream>
#include <cmath>
#include "loadata.h"
#include "segment.h"
#include "frame_seq_fusion.h"

void read_group(string filepath, std::map<int, std::vector<int>>& group, int tube_num)
{
	string filename = filepath + "/result_groups.txt";
	std::ifstream infile(filename.c_str());

	std::string line;
	while (getline(infile, line) && line != "")
	{
		double tube_id, group_id;
		std::istringstream(line) >> tube_id >> group_id;

		group[group_id].emplace_back(tube_id);
	}
}

double cal_colli_area(cv::Vec6d& a, cv::Vec6d& b)
{
	double x1 = a[1], y1 = a[2], w1 = a[3], h1 = a[4];
	double x2 = b[1], y2 = b[2], w2 = b[3], h2 = b[4];

	double endx = max(x1 + w1, x2 + w2);
	double startx = min(x1, x2);
	double width = w1 + w2 - (endx - startx);  // 重叠部分宽
	double endy = max(y1 + h1, y2 + h2);
	double starty = min(y1, y2);
	double height = h1 + h2 - (endy - starty);  // 重叠部分高


	if (width > 0.0 && height > 0.0) {
		//int area = width * height;  // 重叠部分面积
		return width * height;

	}
	else {
		// 不重叠：算出来的width或height小于等于0
		return 0.0;
	}
}


void fusion(Context& context, std::map<int,int>& tubes_group_id, double speed)
{
	double sigma = 0.000093; //0.0013
	
	context.tube_num = 715 - 456 + 1;
	context.scale = 1.0;

	std::map<int, std::vector<int>> groups;

	std::map<int, std::vector<Segment*>> best;

	read_boundingbox(context.filepath + "/boundingbox", context.bbox, context.tube_num, context.scale, context.video_height, context.video_width);
	read_group(context.filepath, groups, context.tube_num);

	std::cout << "read comple" << std::endl;

	// 加速
	for (int i = 0; i < int(context.bbox.size()); ++i)
	{
		int beg = context.bbox[i][0][0];
		int tmp = int(context.bbox[i].size()) - 1;
		int end = context.bbox[i][tmp][0];
		// 计算加速后的新结束帧号
		int newEnd = max(beg, beg + int(round(speed * (end - beg + 1))) - 1);
		
		std::vector<cv::Vec6d> bboxes_cur_tube;
		// 获取加速后，当前tube的各帧的bbox
		for (int kk = beg; kk <= newEnd; ++kk)
		{
			int k;
			if (beg == newEnd)
			{
				k = beg;
			}
			else
			{
				k = (int)(((double)(kk - beg) / (newEnd - beg)) * (end - beg) + beg);
			}
			// 更新帧号
			cv::Vec6d bbox = context.bbox[i][k - beg];
			bbox[0] = kk;
			// 记录原始帧号
			bbox[5] = k;
			bboxes_cur_tube.push_back(bbox);
		}
		context.bbox[i].assign(bboxes_cur_tube.begin(), bboxes_cur_tube.end());
	}
	for (int group_id = 0; group_id < int(groups.size()); ++group_id)
	{
		int tube_id = groups[group_id][0] - 1;
		int a_0 = context.bbox[tube_id][0][0];
		// 后续tube
		for (size_t j = 1; j < groups[group_id].size(); ++j)
		{
			int tube_id = groups[group_id][j] - 1;
			int a = context.bbox[tube_id][0][0];
			// 原本的tube的间隔也要加速
			int gap = int(double(a - a_0) * speed + 0.5);
			int aa = a_0 + gap;
			int tmp = int(context.bbox[tube_id].size()) - 1;
			int b = context.bbox[tube_id][tmp][0];
			int bb = b - (a - aa);
			std::vector<cv::Vec6d> bboxes_cur_tube;
			for (int kk = aa; kk <= bb; ++kk)
			{
				// 更新帧号
				cv::Vec6d bbox = context.bbox[tube_id][kk - aa];
				bbox[0] = kk;
				bbox[5] = kk + a - aa;
				bboxes_cur_tube.push_back(bbox);
			}
			context.bbox[tube_id].assign(bboxes_cur_tube.begin(), bboxes_cur_tube.end());
		}
	}
	

	std::map<int, std::pair<int, int>>FTS; // 各组起始帧
	int max_group_len = 0;
	for (auto it : groups)
	{
		// 获取每个tube的group id
		tubes_group_id[it.second[0]] = it.first;

		FTS[it.first].first = int(context.bbox[it.second[0] - 1][0][0]);

		int tmp = int(context.bbox[it.second[0] - 1].size()) - 1;

		FTS[it.first].second = int(context.bbox[it.second[0] - 1][tmp][0]);

		for (size_t i = 1; i < it.second.size(); ++i)
		{
			// 获取每个tube的group id
			tubes_group_id[it.second[i]] = it.first;

			FTS[it.first].first = std::min(FTS[it.first].first, int(context.bbox[it.second[i] - 1][0][0]));

			int tmp = int(context.bbox[it.second[i] - 1].size()) - 1;

			FTS[it.first].second = std::max(FTS[it.first].second, int(context.bbox[it.second[i] - 1][tmp][0]));

		}

		std::cout << "G beg end:" << FTS[it.first].first << " " << FTS[it.first].second << std::endl;

		max_group_len = max(max_group_len, FTS[it.first].second - FTS[it.first].first + 1);
	}
	std::cout << "max_group_len: " << max_group_len << std::endl;

	int synop_len = FTS[0].second - FTS[0].first + 1;
	std::cout << "synolen: " << synop_len << std::endl;

	std::map<int,std::vector<cv::Vec6d>> f_vel;

	int offset = context.bbox[groups[0][0] - 1][0][0];
	int pre_group_beg = 0; // 当前帧号

	// 处理第0组
	for (size_t i = 0; i < groups[0].size(); ++i)
	{
		std::vector<Segment*> tmp_seg;

		int tube_id = groups[0][i] - 1;

		int a = context.bbox[tube_id][0][0];

		int tmp = int(context.bbox[tube_id].size()) - 1;
		
		int b = context.bbox[tube_id][tmp][0];
		
		tmp_seg.push_back(new Segment(a, b, a - offset, b - offset, tube_id, b - a + 1));

		best[tube_id] = tmp_seg;

		for (size_t j = 0; j < context.bbox[tube_id].size(); ++j)
		{
			f_vel[context.bbox[tube_id][j][0] - offset].push_back(context.bbox[tube_id][j]);
		}
	}

	//std::cout << "group num:" << group.size() << std::endl;

	// 处理后续组
	for (int group_id = 1; group_id < int(groups.size()); ++group_id)
	{
		int group_len = FTS[group_id].second - FTS[group_id].first + 1;
		/*if (group_len > synop_len - 1 - pre_group_beg)
		{
			pl = group_len - pre_group_beg;
		}*/
		int pl = std::min(group_len, synop_len - pre_group_beg);// 可能位置的长度


		double vel_bboxes_area = 0.0;
		for (int frame_idx = pre_group_beg; frame_idx < pre_group_beg + pl; ++frame_idx)
		{
			// 计算帧容器中 pl 段 该帧 全部bbox的面积之和
			for (int bbox_idx = 0; bbox_idx < int(f_vel[frame_idx].size()); ++bbox_idx)
			{
				vel_bboxes_area += f_vel[frame_idx][bbox_idx][3] * f_vel[frame_idx][bbox_idx][4];
			}
			// 减去帧容器 pl 段 原有的碰撞
			for (int bbox_idx = 0; bbox_idx < int(f_vel[frame_idx].size()) - 1; ++bbox_idx)
			{
				for (int bbox_idx2 = bbox_idx + 1; bbox_idx2 < int(f_vel[frame_idx].size()); ++bbox_idx2)
				{
					vel_bboxes_area -= cal_colli_area(f_vel[frame_idx][bbox_idx], f_vel[frame_idx][bbox_idx2]);
				}
			}
		}
		

		// 遍历所有可能的posi
		bool work = false;
		int res_beg = 0;
		for (int tmp_beg = pre_group_beg; tmp_beg < std::min(pre_group_beg + pl, synop_len); ++tmp_beg)
		{
			double ciou = 0.0;

			int tube_id = groups[group_id][0] - 1;
			offset = context.bbox[tube_id][0][0] - tmp_beg;
			// 计算碰撞
			for (int frame_idx = tmp_beg; frame_idx < std::min(tmp_beg + group_len, synop_len); ++frame_idx)
			{
				double colli_area = 0.0, groub_bboxes_area = 0.0;
				
				// 计算 pl 段 该帧 中该group全部bbox的面积和
				for (int j = 0; j < int(groups[group_id].size()); ++j)
				{
					tube_id = groups[group_id][j] - 1;
					int aa_ = context.bbox[tube_id][0][0] - offset;

					int tmp = int(context.bbox[tube_id].size()) - 1;

					int bb_ = context.bbox[tube_id][tmp][0] - offset;

					// 当前tube可能在该帧产生碰撞
					if (frame_idx >= aa_ && frame_idx <= bb_)
					{
						groub_bboxes_area += context.bbox[tube_id][frame_idx - aa_][3] * context.bbox[tube_id][frame_idx - aa_][4];

						cv::Vec6d bbox1 = context.bbox[tube_id][frame_idx - aa_];
						for (int bbox_idx = 0; bbox_idx < int(f_vel[frame_idx].size()); ++bbox_idx)
						{
							colli_area += cal_colli_area(bbox1, f_vel[frame_idx][bbox_idx]);
						}
					}
				}

				// 减去 group 在 pl 段 内部的碰撞
				for (int j = 0; j < int(groups[group_id].size())-1; ++j)
				{
					int tube_id1 = groups[group_id][j] - 1;
					int aa_1 = context.bbox[tube_id1][0][0] - offset;
					int tmp1 = int(context.bbox[tube_id1].size()) - 1;
					int bb_1 = context.bbox[tube_id1][tmp1][0] - offset;

					// 当前 tube 可能在该帧与其他同 group 的 tube 产生碰撞
					if (frame_idx >= aa_1 && frame_idx <= bb_1)
					{
						for (int t = j + 1; t < int(groups[group_id].size()); ++t)
						{
							int tube_id2 = groups[group_id][t] - 1;
							int aa_2 = context.bbox[tube_id2][0][0] - offset;
							int tmp2 = int(context.bbox[tube_id2].size()) - 1;
							int bb_2 = context.bbox[tube_id2][tmp2][0] - offset;

							if (frame_idx >= aa_2 && frame_idx <= bb_2)
							{
								groub_bboxes_area -= cal_colli_area(context.bbox[tube_id1][frame_idx - aa_1], context.bbox[tube_id2][frame_idx - aa_2]);
							}
						}
					}
				}

				/*std::cout << "cgorup area " << groub_bboxes_area << std::endl;*/
				ciou += colli_area / (vel_bboxes_area + groub_bboxes_area - colli_area + 1e-15);
			}

			//std::cout << "ciou: " << ciou << std::endl;
			// 找到了放置posi
			if (ciou / pl <= sigma-1e-15)
			{
				//std::cout << "find" << std::endl;
				res_beg = tmp_beg;
				work = true;
				break;
			}
		}
		// 没有最优位置，添加到摘要后
		if (!work)
		{
			res_beg = synop_len;
		}
		// 更新上一组的位置
		pre_group_beg = res_beg;

		//std::cout << "res_beg: " << res_beg << std::endl;
		synop_len = std::max(synop_len, res_beg + group_len);

		// 记录结果
		offset = context.bbox[groups[group_id][0]-1][0][0] - res_beg;
		for (size_t j = 0; j < groups[group_id].size(); ++j)
		{
			std::vector<Segment*> tmp_seg;

			int tube_id = groups[group_id][j] - 1;

			int a = context.bbox[tube_id][0][0];
			
			int tmp = int(context.bbox[tube_id].size()) - 1;

			int b = context.bbox[tube_id][tmp][0];

			
			tmp_seg.push_back(new Segment(a, b, a - offset, b - offset, tube_id, b - a + 1));

			best[tube_id] = tmp_seg;

			for (size_t t = 0; t < context.bbox[tube_id].size(); ++t)
			{
				f_vel[context.bbox[tube_id][t][0] - offset].push_back(context.bbox[tube_id][t]);
			}

		}

		std::cout << "group " << group_id << " complete!" << std::endl;
	}

	// 更新synoplen
	context.synoplen = synop_len;
	// 输出
	for (auto it : best)
	{
		context.all_segs.push_back(it.second);
	}

	/*for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::vector<Segment*> tmp;
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			tmp.push_back(new Segment(*context.all_segs[i][j]));
		}
		bestCopy.push_back(tmp);
	}*/

	// 计算碰撞
	/*std::map<int, std::vector<pair<int, cv::Vec6d>>>bbox_per_frame;
	for (size_t i = 0; i < context.all_segs.size(); ++i)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); ++j)
		{
			Segment* seg = context.all_segs[i][j];
			for (int frame_idx = seg->aa_; frame_idx <= seg->bb_; ++frame_idx)
			{
				cv::Vec6d patch = context.bbox[i][frame_idx - seg->aa_];
				bbox_per_frame[frame_idx].push_back(pair<int, cv::Vec6d>(i,patch));
			}
		}
	}

	double colli_area = 0.0;
	for (auto it : bbox_per_frame)
	{
		for (int i = 0; i < int(it.second.size()) - 1; ++i)
		{
			for (int j = i + 1; j < int(it.second.size()); ++j)
			{
				// 如果两个tube属于不同的group，计算碰撞
				if (tubes_group_id[it.second[i].first] != tubes_group_id[it.second[j].first])
				{
					colli_area += cal_colli_area(it.second[i].second, it.second[j].second);
				}
			}
		}
	}*/

	// 输出colli_area
	std::ofstream ofile(context.filepath + "/metrics.txt",ofstream::app);
	//ofile << "colli_area: "<< colli_area << "\n";

	std::cout <<"synopsis_len: "<< context.synoplen << std::endl;
	std::cout << "group num: " << groups.size() << std::endl;

	ofile << "sigma: " << sigma << "\n";
}