#include <iostream>
#include <fstream>
#include "loadata.h"
#include "segment.h"
#include "pccva.h"
#include "frame_seq_fusion.h"

// 覆盖率排序
bool my_cmp1(pair<int, double>& a, pair<int, double>& b)
{
	if (a.second < b.second + 1e-25 && a.second > b.second - 1e-25)
	{
		return a.first > b.first;
	}
	return a.second < b.second - 1e-25;
}

// 票数排序
bool my_cmp2(pair<int, double>& a, pair<int, double>& b)
{
	if (a.second < b.second + 1e-25 && a.second > b.second - 1e-25)
	{
		return a.first < b.first;
	}
	return a.second > b.second - 1e-25;
}

double cal_colli_area_cube_tube(cv::Vec6d& a, cv::Vec6d& b, int D)
{
	double x1 = a[0], y1 = a[2], w1 = D, h1 = D;
	double x2 = b[1], y2 = b[2], w2 = b[3], h2 = b[4];

	double endx = max(x1 + w1, x2 + w2);
	double startx = min(x1, x2);
	double width = w1 + w2 - (endx - startx);  // 重叠部分宽
	double endy = max(y1 + h1, y2 + h2);
	double starty = min(y1, y2);
	double height = h1 + h2 - (endy - starty);  // 重叠部分高


	if (width > 0 && height > 0) {
		//int area = width * height;  // 重叠部分面积
		return width * height;

	}
	else {
		// 不重叠：算出来的width或height小于等于0
		return 0.0;
	}
}

bool bb_overlap_for_pccva(cv::Vec6d& a, cv::Vec6d& b)
{
	if (a[0] >= b[1] + b[3])
	{
		return false;
	}

	if (a[2] >= b[2] + b[4])
	{
		return false;
	}

	if (a[1] <= b[1])
	{
		return false;
	}

	if (a[3] <= b[2])
	{
		return false;
	}
	return true;
}

void pccva(Context& context, std::map<int, int>& tubes_group_id, double speed)
{
	std::map<int, std::vector<int>> groups;

	std::map<int, std::vector<Segment*>> best;

	int synop_len = context.synoplen;
	// cube 的边长
	int D = 15;
	double vol = double(D * D * D);
	// 覆盖率阈值
	double gamma = 0.75;

	// tube num
	context.tube_num = 702 - 1 + 1;
	context.scale = 1.0;
	
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
	int group_num_longer_synoplen = 0;
	int tube_num_longer_synoplen = 0;
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
		//std::cout << "G beg end:" << FTS[it.first].first << " " << FTS[it.first].second << std::endl;
		if (FTS[it.first].second - FTS[it.first].first + 1 > context.synoplen)
		{
			++group_num_longer_synoplen;
			tube_num_longer_synoplen += int(it.second.size());
		}
		max_group_len = max(max_group_len, FTS[it.first].second - FTS[it.first].first + 1);
	}
	std::cout << "group_num_longer_synoplen: " << group_num_longer_synoplen << std::endl;
	std::cout << "tube_num_longer_synoplen: " << tube_num_longer_synoplen << std::endl;
	std::cout << "max_group_len: " << max_group_len << std::endl;
	
	// 划分cube
	std::vector<cv::Vec6d>all_cubes;
	for (int i1 = 0, i2 = i1 + D - 1; i2 < context.video_height; i1 += D, i2 += D)
	{
		for (int j1 = 0, j2 = j1 + D - 1; j2 < context.video_width; j1 += D, j2 += D)
		{
			for (int t1 = 0, t2 = t1 + D - 1; t2 < synop_len; t1 += D, t2 += D)
			{
				cv::Vec6d cube;
				cube[0] = j1;
				cube[1] = j2;
				cube[2] = i1;
				cube[3] = i2;
				cube[4] = t1;
				cube[5] = t2;
				all_cubes.push_back(cube);
			}
		}
	}
	std::cout << "cube create comple! cube num: "<<all_cubes.size() << std::endl;


	// 每个cube的group和对应的覆盖率
	std::vector<vector<pair<int,double>>>groups_cover_for_each_cube(all_cubes.size());
	// 每个cube的 “候选” group的数量
	std::vector<int>candi_groups_num_for_each_cube(all_cubes.size(), 0);
	// 每个group的会给它投票的cube
	std::map<int, vector<int>>candi_cubes_for_each_group;
	// 找可能覆盖的cube
	for (int group_id = 0; group_id < int(groups.size()); ++group_id)
	{
		// 可行时间位置下界
		int lbg = 0;
		int group_len = FTS[group_id].second - FTS[group_id].first + 1;
		// 可行时间位置上界
		int ubg = synop_len - group_len;

		//std::cout << "lbg glen ubg" << lbg << " " << group_len << " " << ubg << std::endl;

		std::vector<double>max_coverage(all_cubes.size(), 0.0);

		for (int frame_beg = lbg; frame_beg <= ubg; frame_beg += 10) // ----------------
		{
			std::vector<double>cur_coverage(all_cubes.size(), 0.0);
			// 时间偏移
			//std::cout << "offset " << context.bbox[groups[group_id][0] - 1][0][0] << std::endl;
			int offset = int(context.bbox[groups[group_id][0] - 1][0][0]) - frame_beg;
			
			// 遍历当前组的每个tube
			for (int tube_idx = 0; tube_idx < int(groups[group_id].size()); ++tube_idx)
			{
				int tube_id = groups[group_id][tube_idx];
				//std::cout << "tube " << tube_id << " coverage begin" << std::endl;

				int tube_beg = int(context.bbox[tube_id - 1][0][0]) - offset;
				int tmp = int(context.bbox[tube_id - 1].size()) - 1;
				int tube_end = int(context.bbox[tube_id - 1][tmp][0]) - offset;
				// 遍历cube
				for (int cube_idx = 0; cube_idx < all_cubes.size(); ++cube_idx)
				{
					cv::Vec6d cube = all_cubes[cube_idx];
					int share_frame_beg = max(int(cube[4]), tube_beg);
					int share_frame_end = min(int(cube[5]), tube_end);
					
					// 当前tube与当前cube的重叠面积
					double colli_area = 0.0;
					// 当前tube与当前cube有共享帧
					for (int share_frame_idx = share_frame_beg; share_frame_idx <= share_frame_end; ++share_frame_idx)
					{
						cv::Vec6d tube_patch = context.bbox[tube_id - 1][share_frame_idx - tube_beg];
						colli_area += cal_colli_area_cube_tube(cube, tube_patch, D);
					}
					// 更新在当前时间，当前group对cube的覆盖率
					cur_coverage[cube_idx] += colli_area / vol;
					//std::cout << "coverage: " << cur_coverage[cube_idx] << std::endl;
				}
				//std::cout << "tube " << tube_id << " coverage comple" << std::endl;
			}

			// 更新当前group对各个cube的最大覆盖率
			for (int cube_idx = 0; cube_idx < int(all_cubes.size()); ++cube_idx)
			{
				max_coverage[cube_idx] = max(max_coverage[cube_idx], cur_coverage[cube_idx]);
			}
			//std::cout << "group " << group_id << " max_coverage get" << std::endl;
		}

		// 判断每个cube是否把当前group作为候选group
		for (int cube_idx = 0; cube_idx < int(all_cubes.size()); ++cube_idx)
		{
			groups_cover_for_each_cube[cube_idx].push_back(pair<int, double>(group_id, max_coverage[cube_idx]));
			// 覆盖率>=gamma，当前group是该cube的候选group
			if (max_coverage[cube_idx] >= gamma)
			{
				/*pair<int, double>tmp;
				tmp.first = group_id;
				tmp.second = max_coverage[cube_idx];*/
				
				++candi_groups_num_for_each_cube[cube_idx];
				candi_cubes_for_each_group[group_id].push_back(cube_idx);
			}
		}

		//std::cout << "compute coverage for group" << group_id << " comple" << std::endl;
	}

	// 对每个cube，排序它的候选group
	for (int cube_idx = 0; cube_idx<int(all_cubes.size()); ++cube_idx)
	{
		sort(groups_cover_for_each_cube[cube_idx].begin(), groups_cover_for_each_cube[cube_idx].end(), my_cmp1);
	}

	// 计算票数
	std::vector<pair<int,double>>votes;
	for (int group_id = 0; group_id < int(groups.size()); ++group_id)
	{
		votes.push_back(pair<int, double>(group_id, 0.0));

		for (size_t j = 0; j < candi_cubes_for_each_group[group_id].size(); ++j)
		{
			int cube_idx = candi_cubes_for_each_group[group_id][j];
			int candi_groups_num = candi_groups_num_for_each_cube[cube_idx];

			votes[group_id].second += 1.0 / double(candi_groups_num);
		}
	}

	// 根据票数排序
	sort(votes.begin(), votes.end(), my_cmp2);

	//std::cout << "votes comple " <<votes.size()<< std::endl;

	// 记录每个cube是否被它的主group覆盖
	std::vector<bool>cubes_state(all_cubes.size(), false);
	// 记录每个group是否已经被放置
	std::vector<bool>groups_state(groups.size(), false);
	// 根据票数开始重排
	int max_bb = 0;
	while (true)
	{
		std::vector<pair<int, double>>dropped_groups;

		for (size_t i = 0; i < votes.size(); ++i)
		{
			int group_id = votes[i].first;

			vector<int>Cp;
			// 找Cp
			for (int cube_idx = 0; cube_idx < int(all_cubes.size()); ++cube_idx)
			{
				// 该cube已经被它的主group占据
				if (cubes_state[cube_idx])
				{
					continue;
				}

				for (int tmp = int(groups_cover_for_each_cube[cube_idx].size()) - 1;tmp >= 0;--tmp)
				{
					int tmp_group_id = groups_cover_for_each_cube[cube_idx][tmp].first;
					// 覆盖率比当前group大的group都放置好了
					if (groups_state[tmp_group_id])
					{
						groups_cover_for_each_cube[cube_idx].erase(groups_cover_for_each_cube[cube_idx].begin() + tmp);
						continue;
					}
					else
					{
						// 当前group是当前cube的主group
						if (tmp_group_id == group_id)
						{
							Cp.push_back(cube_idx);
						}
						// 当前group不是当前cube的主group
						break;
					}
				}
			}

			// Lp中每个时间位置-覆盖了哪些Cp中的cube
			std::map<int, std::set<int>>Lp;
			// 找Lp
			// 可行时间位置下界
			int lbg = 0;
			int group_len = FTS[group_id].second - FTS[group_id].first + 1;
			// 可行时间位置上界
			int ubg = synop_len - group_len;

			for (int frame_beg = lbg; frame_beg <= ubg; ++frame_beg)
			{
				// 遍历Cp中的cube
				for (int cp_idx = 0; cp_idx < int(Cp.size()); ++cp_idx)
				{
					int cube_idx = Cp[cp_idx];
					if (Lp.count(frame_beg) > 0 && Lp[frame_beg].count(cube_idx) > 0)
					{
						continue;
					}

					// 时间偏移
					int offset = context.bbox[groups[group_id][0] - 1][0][0] - frame_beg;

					// 遍历当前组的每个tube
					bool work = false;
					for (int tube_idx = 0; tube_idx < int(groups[group_id].size()); ++tube_idx)
					{
						int tube_id = groups[group_id][tube_idx];

						int tube_beg = context.bbox[tube_id - 1][0][0] - offset;
						int tmp = int(context.bbox[tube_id - 1].size()) - 1;
						int tube_end = context.bbox[tube_id - 1][tmp][0] - offset;

						cv::Vec6d cube = all_cubes[cube_idx];
						int share_frame_beg = max(int(cube[4]), tube_beg);
						int share_frame_end = min(int(cube[5]), tube_end);

						// 当前tube与当前cube有共享帧
						for (int share_frame_idx = share_frame_beg; share_frame_idx <= share_frame_end; ++share_frame_idx)
						{
							cv::Vec6d tube_patch = context.bbox[tube_id - 1][share_frame_idx - tube_beg];
							if (bb_overlap_for_pccva(cube, tube_patch))
							{
								Lp[frame_beg].insert(cube_idx);
								work = true;
								break;
							}
						}

						if (work)
						{
							break;
						}

					}
				}
			}

			if (Lp.empty())
			{
				std::cout << "group " << group_id << " dont have Lp!" << std::endl;
				dropped_groups.push_back(votes[i]);
				continue;
			}

			// 最大覆盖Cp数
			int max_cover_Cp_num = -1;
			// 该组最终结果放置位置
			int res_frame_beg = -1;
			// 最大覆盖率
			double max_coverage = 0.0;
			for (auto it : Lp)
			{
				int frame_beg = it.first;
				int cover_Cp_num = int(it.second.size());
				// 更新最大覆盖Cp数 和 该组最终结果放置位置
				if (cover_Cp_num > max_cover_Cp_num)
				{
					max_cover_Cp_num = cover_Cp_num;
					res_frame_beg = frame_beg;

					// 更新最大覆盖率
					double coverage = 0.0;
					for (set<int>::iterator cp_it = it.second.begin(); cp_it != it.second.end(); cp_it++)
					{
						for (int tmp = int(groups_cover_for_each_cube[*cp_it].size()) - 1; tmp >= 0; --tmp)
						{
							if (groups_cover_for_each_cube[*cp_it][tmp].first == group_id)
							{
								coverage += groups_cover_for_each_cube[*cp_it][tmp].second;
								break;
							}
						}
					}

					max_coverage = coverage;
				}

				else
				{
					// 对于有相同覆盖Cp数量的时间位置，取覆盖率大的位置
					if (cover_Cp_num == max_cover_Cp_num)
					{
						double coverage = 0.0;
						for (set<int>::iterator cp_it = it.second.begin(); cp_it != it.second.end(); cp_it++)
						{
							for (int tmp = int(groups_cover_for_each_cube[*cp_it].size()) - 1; tmp >= 0; --tmp)
							{
								if (groups_cover_for_each_cube[*cp_it][tmp].first == group_id)
								{
									coverage += groups_cover_for_each_cube[*cp_it][tmp].second;
									break;
								}
							}
						}

						if (coverage > max_coverage - 1e-25)
						{
							max_coverage = coverage;
							res_frame_beg = frame_beg;
						}
					}
				}
			}

			for (std::set<int>::iterator it = Lp[res_frame_beg].begin(); it != Lp[res_frame_beg].end(); ++it)
			{
				cubes_state[*it] = true;
			}

			// 更新Cp中cube的主组
			/*for (size_t cp_idx = 0; cp_idx < Cp.size(); ++cp_idx)
			{
				int cube_idx = Cp[cp_idx];
				// 删掉最后一个group，也就是当前的group
				groups_cover_for_each_cube[cube_idx].pop_back();
			}*/

			// 记录结果
			int offset = context.bbox[groups[group_id][0] - 1][0][0] - res_frame_beg;
			
			for (size_t j = 0; j < groups[group_id].size(); ++j)
			{
				std::vector<Segment*> tmp_seg;

				int tube_id = groups[group_id][j] - 1;

				int a = context.bbox[tube_id][0][0];

				int tmp = int(context.bbox[tube_id].size()) - 1;

				int b = context.bbox[tube_id][tmp][0];

				tmp_seg.push_back(new Segment(a, b, a - offset, b - offset, tube_id, b - a + 1));

				best[tube_id+1] = tmp_seg;

				max_bb = max(max_bb, b - offset);
			}

			// group_id 被放置
			groups_state[group_id] = true;
			std::cout << "group " << group_id << " complete!" << std::endl;
		}

		if (dropped_groups.empty() || votes.size() == dropped_groups.size())
		{
			break;
		}

		votes.assign(dropped_groups.begin(), dropped_groups.end());
		std::cout << "one step rearrange complete" << std::endl;
	}
	std::cout << "max_bb: " << max_bb << std::endl;


	// 写入all_segs
	/*for (auto it : best)
	{
		context.all_segs.push_back(it.second);
	}*/
	int syn_tubes_num = 0;
	for (int tube_id = 1; tube_id <= context.tube_num; ++tube_id)
	{
		if (best.count(tube_id) > 0)
		{
			++syn_tubes_num;
			context.all_segs.push_back(best[tube_id]);
		}
		else
		{
			std::vector<Segment*> tmp_seg;
			tmp_seg.clear();
			context.all_segs.push_back(tmp_seg);
		}
	}

	std::cout << "write complete." <<" dropped tubes nums: "<< context.tube_num - syn_tubes_num << std::endl;

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
				bbox_per_frame[frame_idx].push_back(pair<int, cv::Vec6d>(i, patch));
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
	std::ofstream ofile(context.filepath + "/metrics.txt", ofstream::app);
	//ofile << "colli_area: " << colli_area << "\n";
	ofile << "D: " << D << "\n";
	ofile << "number of dropped tube: " << context.tube_num - syn_tubes_num << "\n";
	
	std::cout << "group num: " << groups.size() << std::endl;
}