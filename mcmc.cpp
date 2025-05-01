#include "mcmc.h"
#include "graph.h"
#include "energy.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

#define __DEBUG__

double random_scale(cv::RNG& rng, Segment* seg)
{
	double min_scale = 0.8;
	//double min_scale = 1;
	double max_scale = 1;

	double r = rng.gaussian(0.1);
	double ans = seg->scale_ + r;
	
	return mymin(mymax(ans, min_scale), max_scale);
}

int random_aa(cv::RNG& rng, Segment* seg, double sigma, int synlen, double min_s, double max_s, int fps,std::vector<std::vector<cv::Vec6d>>& bbox, Context &context)
{
	const double EPS = 1e-15;

	/*int tube_beg = bbox[seg->tube_id_][0][0];
	double speed = 0.0;
	for (int i = seg->a_; i <= seg->b_ - 10; ++i)
	{
		int x1 = bbox[seg->tube_id_][i - tube_beg][1], y1 = bbox[seg->tube_id_][i - tube_beg][2];
		int x2 = bbox[seg->tube_id_][i + 10 - tube_beg][1], y2 = bbox[seg->tube_id_][i + 10 - tube_beg][2];
		double dis = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) + 1e-15;
		speed += dis / 10.0;
	}
	speed /= double(seg->b_ - 10 - seg->a_ + 1);*/
	double speed = context.tube_speed[seg->tube_id_][seg->a_];
	
	int seg_len = seg->b_ - seg->a_ + 1;
	//double min_scale = 0.1;
	//double max_scale = 10000;
	
	//double min_scale = 0.1; // 4.5是定量结果用的参数
	double min_scale = speed / 4.5; // 3.0是相同collision的结果
	//double min_scale = 0.15;

	double max_scale = speed / (double(seg_len) / 100.0 + 1.0);
	//double max_scale = speed / (double(seg_len) / 130.0 + 1.0);
	//double max_scale = 1;
	//double max_scale = 8;

	/*if (seg_len > fps * 2)
	{
		max_scale = speed / (double(seg_len) / 100.0 + 2.1);
	}*/

	if (min_scale > 1.0 - EPS)
	{
		min_scale = 1.0;
	}

	if (max_scale < min_scale-EPS)
	{
		std::swap(min_scale, max_scale);
	}

	/*if (seg->b_ - seg->a_ + 1 < seg_avg_len)
	{
		max_scale = double(seg->b_ - seg->a_) / 10.0 - 1.0;
		if (max_scale > max_s - EPS)
		{
			max_scale = max_s;
		}
	}*/

	double r = rng.gaussian(sigma);
	int aa = int(seg->aa_ + r);
	int mmax = seg->bb_ - int(min_scale * seg->len_);
	int mmin = seg->bb_ - int(max_scale * seg->len_);
	return mymax(mymin(mmax, aa), mmin);
	/*mmax = mymin(mmax, synlen-1);
	mmin = mymax(0, mmin);
	if (mmax < mmin)
	{
		return seg->aa_;
	}
	else
	{
		return mymax(mmin, mymin(aa, mmax));
	}*/
	//return mymin(synlen, mymax(0, mm));
	//return rng.uniform(seg->bb_ - int(5 * seg->len_), seg->bb_ - int(0.1 * seg->len_));
}

int random_bb(cv::RNG& rng, Segment* seg, double sigma, int synlen, double min_s, double max_s, int fps, std::vector<std::vector<cv::Vec6d>>& bbox,Context &context)
{
	const double EPS = 1e-15;

	/*int tube_beg = bbox[seg->tube_id_][0][0];
	double speed = 0.0;
	for (int i = seg->a_; i <= seg->b_ - 10; ++i)
	{
		int x1 = bbox[seg->tube_id_][i - tube_beg][1], y1 = bbox[seg->tube_id_][i - tube_beg][2];
		int x2 = bbox[seg->tube_id_][i + 10 - tube_beg][1], y2 = bbox[seg->tube_id_][i + 10 - tube_beg][2];
		double dis = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) + 1e-15;
		speed += dis / 10.0;
	}
	speed /= double(seg->b_ - 10 - seg->a_ + 1);*/
	double speed = context.tube_speed[seg->tube_id_][seg->a_];

	int seg_len = seg->b_ - seg->a_ + 1;
	//double min_scale = 0.1;
	//double max_scale = 10000;
	
	double min_scale = speed / 4.5; // 3.0是相同collision的结果
	//double min_scale = 0.1;

	double max_scale = speed / (double(seg_len) / 100.0 + 1.0);
	//double max_scale = speed / (double(seg_len) / 130.0 + 1.0);
	//double max_scale = 8;

	/*if (seg_len > fps * 2)
	{
		max_scale = speed / (double(seg_len) / 90.0 + 2.1);
	}*/

	if (min_scale > 1.0 - EPS)
	{
		min_scale = 1.0;
	}

	if (max_scale < min_scale - EPS)
	{
		std::swap(min_scale, max_scale);
	}

	/*if (seg->b_ - seg->a_ + 1 < seg_avg_len)
	{
		max_scale = double(seg->b_ - seg->a_) / 10.0 - 1.0;
		if (max_scale > max_s - EPS)
		{
			max_scale = max_s;
		}
	}*/

	/*if (max_scale > max_s - EPS)
	{
		max_scale = max_s;
	}*/

	double r = rng.gaussian(sigma);
	int bb = int(seg->bb_ + r);
	int mmin = seg->aa_ + int(min_scale * seg->len_);
	int mmax = seg->aa_ + int(max_scale * seg->len_);
	return mymax(mymin(mmax, bb), mmin);
	/*mmax = mymin(mmax, synlen - 1);
	mmin = mymax(0, mmin);
	if (mmax < mmin)
	{
		return seg->bb_;
	}
	else
	{
		return mymax(mmin, mymin(bb, mmax));
	}*/

	//return mymin(synlen, mymax(0, mm));
	//return rng.uniform(seg->aa_ + int(0.1 * seg->len_), seg->aa_ + int(5 * seg->len_));
}

int random_aa_colli(cv::RNG& rng, Segment* seg, double sigma, int synlen)
{
	double min_scale = 0.55; //0.15
	double max_scale = 0.9; //0.2

	/*double min_scale = 0.2;
	double max_scale = 5;*/

	double r = rng.gaussian(sigma);
	int aa = int(seg->aa_ + r);
	int mmax = seg->bb_ - int(min_scale * seg->len_);
	int mmin = seg->bb_ - int(max_scale * seg->len_);
	return mymax(mymin(mmax, aa), mmin);
	/*mmax = mymin(mmax, synlen - 1);
	mmin = mymax(0, mmin);
	if (mmax < mmin)
	{
		return seg->aa_;
	}
	else
	{
		return mymax(mmin, mymin(aa, mmax));
	}*/
	//return seg->aa_;

	//return rng.uniform(seg->bb_ - int(5 * seg->len_), seg->bb_ - int(0.1 * seg->len_));
}

int random_bb_colli(cv::RNG& rng, Segment* seg, double sigma, int synlen)
{
	double min_scale = 0.55; //0.15
	double max_scale = 0.9; //0.2

	/*double min_scale = 0.2;
	double max_scale = 5;*/

	double r = rng.gaussian(sigma);
	int bb = int(seg->bb_ + r);
	int mmin = seg->aa_ + int(min_scale * seg->len_);
	int mmax = seg->aa_ + int(max_scale * seg->len_);
	return mymax(mymin(mmax, bb), mmin);
	/*mmax = mymin(mmax, synlen - 1);
	mmin = mymax(0, mmin);
	if (mmax < mmin)
	{
		return seg->bb_;
	}
	else
	{
		return mymax(mmin, mymin(bb, mmax));
	}*/
	//return mymin(synlen, mymax(0, mm));

	//return seg->bb_;

	//return rng.uniform(seg->aa_ + int(0.1 * seg->len_), seg->aa_ + int(5 * seg->len_));
}

int random_aaa(cv::RNG& rng, Segment* seg, bool follow)
{
	Segment* before = seg->last_;
	if (before->last_ == NULL) //head
		before = NULL;
	Segment* after = seg->next_;

	/*EdgeNode* edge = seg->firstedge;

	while (edge != NULL)
	{
		if (edge->seg->tube_id_ == seg->tube_id_)
		{
			if (edge->seg->a_ > seg->a_)
			{
				after = edge->seg;
			}
			else
			{
				before = edge->seg;
			}
		}
		edge = edge->next;
	}*/



	if (follow == true)
	{
		before = NULL;
	}

	int left1, right1;

	if (before == NULL && after == NULL)
	{
		left1 = -100000;
		right1 = seg->bb_;
	}
	else if (before == NULL && after != NULL)
	{
		left1 = -100000;
		right1 = mymin(seg->bb_, after->aa_);
	}
	else if (before != NULL && after == NULL)
	{
		left1 = before->aa_;
		right1 = mymin(before->bb_ + 1, seg->bb_);
	}
	else
	{
		left1 = before->aa_;
		right1 = mymin(before->bb_ + 1, after->aa_, seg->bb_);
	}

	double minspeed = 0.3;
	double maxspeed = 6.0;

	double maxlen = seg->len_ / minspeed;
	double minlen = seg->len_ / maxspeed;

	int left2 = int(seg->bb_ - maxlen);
	int right2 = int(seg->bb_ - minlen);

	int left = mymax(left1, left2);
	int right = mymin(right1, right2);

	if (left > right)
	{
		return seg->aa_;
	}
	else
	{
		return rng.uniform(left, right + 1);
	}
}


int random_bbb(cv::RNG& rng, Segment* seg, bool follow)
{
	Segment* before = seg->last_;

	if (before->last_ == NULL) // head
		before = NULL;

	Segment* after = seg->next_;

	/*EdgeNode* edge = seg->firstedge;

	while (edge != NULL)
	{
		if (edge->seg->tube_id_ == seg->tube_id_)
		{
			if (edge->seg->a_ > seg->a_)
			{
				after = edge->seg;
			}
			else
			{
				before = edge->seg;
			}
		}
		edge = edge->next;
	}*/

	if (follow == true)
	{
		after = NULL;
	}

	int left1, right1;

	if (before == NULL && after == NULL)
	{
		left1 = seg->aa_;
		right1 = 100000;
	}
	else if (before == NULL && after != NULL)
	{
		left1 = mymax(seg->aa_, after->aa_-1);
		right1 = after->bb_;
	}
	else if (before != NULL && after == NULL)
	{
		left1 = mymax(seg->aa_, before->bb_);
		right1 = 100000;
	}
	else
	{
		left1 = mymax(seg->aa_, before->bb_, after->aa_ - 1);
		right1 = after->bb_;
	}

	double minspeed = 0.3;
	double maxspeed = 6.0;

	double maxlen = seg->len_ / minspeed;
	double minlen = seg->len_ / maxspeed;

	int left2 = int(seg->aa_ + minlen);
	int right2 = int(seg->aa_ + maxlen);

	// left1 and right1 has to be satisfy, otherwise constraint is ruined

	int left = mymax(left1, left2);
	int right = mymin(right1, right2);

	//std::cout << left << "," << right << "               fdafda" << std::endl;

	if (left > right)
	{
		return seg->bb_;
	}
	else
	{
		int ret = rng.uniform(left, right + 1);
		//std::cout << ret << std::endl;
		//std::cout << double(ret - left) / (right - left) << std::endl;
		return ret;
		//return rng.uniform(left, right + 1);
	}
}
//
//void mcmc(Context& context, std::vector<std::vector<Segment*>>& bestCopy)
//{
//	double enout;
//	double encolli;
//	double enscale;
//	std::vector<double> tube_energy;
//	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };
//	
//	int root_num = int(context.roots.size());
//	int tube_num = int(context.all_segs.size());
//
//	for (size_t i = 0; i < context.all_segs.size(); i++)
//	{
//		std::vector<Segment*> tmp;
//		for (size_t j = 0; j < context.all_segs[i].size(); j++)
//		{
//			tmp.push_back(new Segment(*context.all_segs[i][j]));
//		}
//		bestCopy.push_back(tmp);
//	}
//
//	double hot = 1e23;
//	cv::RNG rng;
//
//	double enself;
//	double oenergy = compute_energy(context, tube_energy, enout, encolli, enscale, enself);
//	double best_energy = oenergy;
//
//
//	std::vector<std::vector<Segment*>> Acopy;
//	for (size_t j = 0; j < context.all_segs.size(); j++)
//	{
//		std::vector<Segment*> tmp;
//		for (size_t k = 0; k < context.all_segs[j].size(); k++)
//		{
//			tmp.push_back(new Segment(*context.all_segs[j][k]));
//		}
//		Acopy.push_back(tmp);
//	}
//
//
//	for (int i = 0; i < 100000; i++)
//	{
//		for (size_t j = 0; j < context.all_segs.size(); j++)
//		{
//			for (size_t k = 0; k < context.all_segs[j].size(); k++)
//			{
//				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
//				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
//			}
//		}
//
//		int tube_id = rand() % tube_num;
//
//		int seg_id = rand() % context.all_segs[tube_id].size();
//
//		Segment* seg = context.all_segs[tube_id][seg_id];
//		int oaa = seg->aa_;
//		int obb = seg->bb_;
//
//		int naa, nbb;
//		int C = rand() % 8;
//
//		if (!seg->colli_)
//		{
//			switch (C)
//			{
//			case 0:
//				naa = random_aa(rng, seg, false);
//				seg->aa_ = naa;
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				break;
//			case 1:
//				naa = random_aa(rng, seg, true);
//				seg->aa_ = naa;
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				reset_traverse(context.all_segs);
//				bfs_change_aabb(context, seg, oaa, obb, false);
//				break;
//			case 2:
//				nbb = random_bb(rng, seg, false);
//				seg->bb_ = nbb;
//				break;
//			case 3:
//				nbb = random_bb(rng, seg, true);
//				seg->bb_ = nbb;
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				reset_traverse(context.all_segs);
//				bfs_change_aabb(context, seg, oaa, obb, false);
//				break;
//			case 4:
//				naa = random_aa(rng, seg, false);
//				nbb = random_bb(rng, seg, false);
//				seg->aa_ = naa;
//				seg->bb_ = nbb;
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				break;
//			case 5:
//				naa = random_aa(rng, seg, true);
//				seg->aa_ = naa;
//				reset_traverse(context.all_segs);
//				bfs_change_aabb(context, seg, oaa, obb, false);
//				nbb = random_bb(rng, seg, false);
//				seg->bb_ = nbb;
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				break;
//			case 6:
//				
//				nbb = random_bb(rng, seg, true);
//				seg->bb_ = nbb;
//				reset_traverse(context.all_segs);
//				bfs_change_aabb(context, seg, oaa, obb, false);
//				naa = random_aa(rng, seg, false);
//				seg->aa_ = naa;
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				break;
//			case 7:
//				naa = random_aa(rng, seg, true);
//				seg->aa_ = naa;
//				reset_traverse(context.all_segs);
//				bfs_change_aabb(context, seg, oaa, obb, false);
//				nbb = random_bb(rng, seg, true);
//				seg->bb_ = nbb;
//				reset_traverse(context.all_segs);
//				bfs_change_aabb(context, seg, naa, obb, false);
//				//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//				break;
//			}
//		}
//		else // if seg collides, others must follow to prevent constraint broken
//		{
//			naa = random_aa(rng, seg, true);
//			seg->aa_ = naa;
//			reset_traverse(context.all_segs);
//			bfs_change_aabb(context, seg, oaa, obb, false);
//			nbb = random_bb(rng, seg, true);
//			seg->bb_ = nbb;
//			reset_traverse(context.all_segs);
//			bfs_change_aabb(context, seg, naa, obb, false);
//			////std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//		}
//		double enself;
//		double nenergy = compute_energy(context, tube_energy, enout, encolli, enscale,enself);
//		double eee1 = nenergy;
//		double eee2 = oenergy;
//
//		int divid = -1;
//		for (int tt = 0; tt < 20; tt++)
//		{
//			if (eee1 > division[tt])
//			{
//				divid = tt;
//			}
//		}
//
//		if (divid < 0 || divid >= 18)
//		{
//			std::cout << "Something wrong!" << std::endl;
//		}
//
//		divid++;
//		eee1 = eee1 / division[divid] + 10;
//		eee2 = eee2 / division[divid] + 10;
//
//		eee1 = pow(eee1, 25);
//		eee2 = pow(eee2, 25);
//
//		double a1 = exp(-1 / hot * eee1);
//		double a2 = exp(-1 / hot * eee2);
//
//		double alpha = a1 / a2;
//
//		alpha = alpha < 1 ? alpha : 1;
//
//		double rd = (rand() % 1001) / 1000.0;
//
//		if (i % 1000 == 0)
//		{
//			std::cout << i << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
//			if (rd < alpha)
//				std::cout << "  Accept" << std::endl;
//			else
//				std::cout << "  Refuse" << std::endl;
//		}
//
//		if (rd < alpha) // accept
//		{
//			oenergy = nenergy;
//
//			if (nenergy < best_energy)
//			{
//				best_energy = nenergy;
//
//				for (size_t u = 0; u < context.all_segs.size(); u++)
//				{
//					for (size_t v = 0; v < context.all_segs[u].size(); v++)
//					{
//						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
//						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
//						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
//						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
//					}
//				}
//			}
//		}
//		else
//		{
//			for (size_t j = 0; j < Acopy.size(); j++)
//			{
//				for (size_t k = 0; k < Acopy[j].size(); k++)
//				{
//					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
//					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
//				}
//			}
//
//			//naa = seg->aa_;
//			//nbb = seg->bb_;
//
//			//if (!seg->colli_)
//			//{
//			//	switch (C)
//			//	{
//			//	case 0:
//			//		
//			//		seg->aa_ = oaa;
//			//		break;
//			//	case 1:
//			//		seg->aa_ = oaa;
//			//		reset_traverse(context.all_segs);
//			//		bfs_change_aabb(context, seg, naa, nbb, false);
//			//		break;
//			//	case 2:
//			//		seg->bb_ = obb;
//			//		break;
//			//	case 3:
//			//		seg->bb_ = obb;
//			//		reset_traverse(context.all_segs);
//			//		bfs_change_aabb(context, seg, naa, nbb, false);
//			//		break;
//			//	case 4:
//			//		seg->aa_ = oaa;
//			//		seg->bb_ = obb;
//			//		break;
//			//	case 5:
//			//		seg->bb_ = obb;
//			//		seg->aa_ = oaa;
//			//		reset_traverse(context.all_segs);
//			//		bfs_change_aabb(context, seg, naa, obb, false);
//			//		
//			//		break;
//			//	case 6:
//			//		seg->aa_ = oaa;
//			//		seg->bb_ = obb;
//			//		reset_traverse(context.all_segs);
//			//		bfs_change_aabb(context, seg, oaa, nbb, false);
//			//		break;
//			//	case 7:
//			//		seg->bb_ = obb;
//			//		reset_traverse(context.all_segs);
//			//		bfs_change_aabb(context, seg, naa, nbb, false);
//			//		reset_traverse(context.all_segs);
//			//		seg->aa_ = oaa;
//			//		bfs_change_aabb(context, seg, naa, obb, false);
//			//		break;
//			//	}
//			//}
//			//else // if seg collides, others must follow to prevent constraint broken
//			//{
//			//	/*seg->bb_ = obb;
//			//	reset_traverse(context.all_segs);
//			//	bfs_change_aabb(context, seg, naa, nbb, false);
//			//	seg->aa_ = oaa;
//			//	reset_traverse(context.all_segs);
//			//	bfs_change_aabb(context, seg, naa, obb, false);*/
//			//}
//		}
//	}
//
//	// Delete copy
//	for (size_t j = 0; j < Acopy.size(); j++)
//	{
//		//std::cout << j << std::endl;
//		for (size_t k = 0; k < Acopy[j].size(); k++)
//		{
//			delete Acopy[j][k];
//		}
//	}
//
//	std::cout << "Best energy:" << best_energy << std::endl;
//}

bool allin_range(Context& context, int tube_id)
{
	if (context.circulated_tube[tube_id] == 1)
	{
		return true;
	}

	for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
	{
		if (context.all_segs[tube_id][i]->aa_ < 0 || context.all_segs[tube_id][i]->aa_ >= context.synoplen ||
			context.all_segs[tube_id][i]->bb_ < 0 || context.all_segs[tube_id][i]->bb_ >= context.synoplen)
		{
			return false;
		}
	}
	return true;
}

bool allin_range(Context& context)
{	
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		/*if (context.circulated_tube[i] == 1)
			continue;*/
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			if (context.all_segs[i][j]->aa_ < 0 || context.all_segs[i][j]->aa_ >= context.synoplen ||
				context.all_segs[i][j]->bb_ < 0 || context.all_segs[i][j]->bb_ >= context.synoplen)
			{
				return false;
			}
		}
	}
	
	return true;
}

bool chronological_order_preserved(Context& context, int slack)
{
	int tube_num = int(context.all_segs.size());
	for (int i = 0; i < tube_num - 1; i++)
	{
		for (int j = i + 1; j < tube_num; j++)
		{
			/*if (context.circulated_tube[i] == 1 ||
				context.circulated_tube[j] == 1)
				continue;*/

			Segment* s1 = context.all_segs[i][0];
			Segment* s2 = context.all_segs[j][0];

			int ds1 = s1->a_ - s2->a_;
			int ds2 = s1->aa_ - s2->aa_;
			if (ds1 * ds2 < 0 && abs(ds2) > slack)
			{
				return false;
			}
		}
	}
	return true;
}

bool beg_chronological_order_preserved(Context& context, int slack)
{
	int tube_num = int(context.all_segs.size());
	for (int i = 0; i < tube_num - 1; i++)
	{
		for (int j = i + 1; j < tube_num; j++)
		{
			/*if (context.circulated_tube[i] == 1 ||
				context.circulated_tube[j] == 1)
				continue;*/

			Segment* s1 = context.all_segs[i][0];
			Segment* s2 = context.all_segs[j][0];

			int ds1 = s1->a_ - s2->a_;
			int ds2 = s1->aa_ - s2->aa_;
			if (abs(abs(ds1) - abs(ds2)) > slack)
			{
				return false;
			}
		}
	}
	return true;
}

double show_max_scale(Context& context)
{
	double max_scale = -1;
	double min_len = 10000;
	int num_len_1 = 0;
	int tube_num = int(context.all_segs.size());
	for (int i = 0; i < tube_num - 1; i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			Segment* s = context.all_segs[i][j];
			/*if (s->len_ != s->b_ - s->a_+1)
			{
				std::cout << "wrong*******************" << std::endl;
			}*/
			if (!s->colli_segs_.empty())
			{
				if (s->len_ != s->bb_ - s->aa_ + 1)
				{
					std::cout << "wrong*******************" << std::endl;
				}
			}
			double scale = double(s->bb_ - s->aa_+1) / double(s->b_ - s->a_+1);
			if (scale > max_scale)
			{
				max_scale = scale;
			}
			if (s->len_ < min_len)
			{
				min_len = s->len_;
			}
			if(s->len_ == 1)
			{
				num_len_1 += 1;
			}
		}
	}
	std::cout << "Min len:" << min_len << std::endl;
	std::cout << "Num len 1:" << num_len_1 << std::endl;
	return max_scale;
}

int mcmc_chrono(Context& context, std::vector<std::vector<Segment*>>& bestCopy, double& min_s, double& max_s)
{
	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };

	int tube_num = int(context.all_segs.size());

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::vector<Segment*> tmp;
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			tmp.push_back(new Segment(*context.all_segs[i][j]));
		}
		bestCopy.push_back(tmp);
	}

	double hot = 1.5e23;
	cv::RNG rng;

	double oenergy = compute_energy_chrono_reverse3(context);
	double best_energy = oenergy;


	std::vector<std::vector<Segment*>> Acopy;
	for (size_t j = 0; j < context.all_segs.size(); j++)
	{
		std::vector<Segment*> tmp;
		for (size_t k = 0; k < context.all_segs[j].size(); k++)
		{
			tmp.push_back(new Segment(*context.all_segs[j][k]));
		}
		Acopy.push_back(tmp);
	}

	int iter_chrono = 1000000;
	for (int itern = 0; itern < 1; itern++)
	{
		for (size_t j = 0; j < context.all_segs.size(); j++)
		{
			for (size_t k = 0; k < context.all_segs[j].size(); k++)
			{
				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
			}
		}

		/* Two options:
		 * 1. Find a pair of tubes, and change their first segment's position
		 * 2. Find a tube, and change its aabb
		 */


		int opt = rand() % 1000;
		if (opt < 200) // the first option
		{
			std::vector<std::pair<int, int>> revert_pairs;
			std::vector<int> revert_degree;
			for (int i = 0; i < tube_num - 1; i++)
			{
				for (int j = i + 1; j < tube_num; j++)
				{
					/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
					{
						continue;
					}*/

					Segment* s1 = context.all_segs[i][0];
					Segment* s2 = context.all_segs[j][0];
					int ds1 = s1->a_ - s2->a_;
					int ds2 = s1->aa_ - s2->aa_;
					// 两个tube的first segment's position 逆转了，且逆转帧数超过 chrono_slack
					if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
					{
						revert_pairs.push_back(std::pair<int, int>(i, j));
						revert_degree.push_back(abs(ds2));
					}
				}
			}

			if (revert_pairs.empty())
			{
				std::cout << "Find the solution. There is no reverted tubes." << std::endl;
				iter_chrono = itern + 1;
				break;
			}

			/*if (fabs(oenergy) < 1e-15)
			{
				std::cout << "Find the solution. There is no reverted tubes." << std::endl;
				iter_chrono = itern + 1;
				break;
			}*/



			/*int se_pair_id = -1;
			int max_pair_revert = -100;
			for (size_t i = 0; i < colli_pairs.size(); i++)
			{
				if (revert_degree[i] > max_pair_revert)
				{
					max_pair_revert = revert_degree[i];
					se_pair_id = int(i);
				}
			}*/

			int se_pair_id = rand() % revert_pairs.size();

			if (itern % 1000 == 0)
			{
				std::cout << "Number of reverted pairs:" << revert_pairs.size() << "............................" << std::endl;
				//std::cout << max_pair_revert << std::endl;
			}

			int tubeid1 = revert_pairs[se_pair_id].first;
			int tubeid2 = revert_pairs[se_pair_id].second;


			int t1oaa = context.all_segs[tubeid1][0]->aa_;
			int t1obb = context.all_segs[tubeid1][0]->bb_;
			int t1len = t1obb - t1oaa;

			int t2oaa = context.all_segs[tubeid2][0]->aa_;
			int t2obb = context.all_segs[tubeid2][0]->bb_;
			int t2len = t2obb - t2oaa;

			int opt2 = rand() % 2;

			if (opt2 == 0)
			{
				reset_traverse(context.heads);
				context.all_segs[tubeid1][0]->aa_ = t2oaa; // tube1 first segment's position 改为tube2的
				context.all_segs[tubeid1][0]->bb_ = t2oaa + t1len;
				dfs_change_aabb(context.all_segs[tubeid1][0], t1oaa, t1obb);
				check_consistency(context);
			}
			else
			{
				reset_traverse(context.heads);
				context.all_segs[tubeid2][0]->aa_ = t1oaa; // tube2 first segment's position 改为tube1的
				context.all_segs[tubeid2][0]->bb_ = t1oaa + t2len;
				dfs_change_aabb(context.all_segs[tubeid2][0], t2oaa, t2obb);
				check_consistency(context);
			}
		}
		else
		{
			std::vector<int> revert_tubeids;

			for (int i = 0; i < tube_num; i++)
			{
				for (int j = 0; j < tube_num; j++)
				{
					if (i == j/* || context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1*/)
						continue;
					Segment* s1 = context.all_segs[i][0];
					Segment* s2 = context.all_segs[j][0];
					int ds1 = s1->a_ - s2->a_;
					int ds2 = s1->aa_ - s2->aa_;
					if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
					{
						revert_tubeids.push_back(i);
						break;
					}
				}
			}

			if (revert_tubeids.empty())
			{
				std::cout << "Find the solution2" << std::endl;
				iter_chrono = itern + 1;
				break;
			}

			int tube_id = revert_tubeids[rand() % revert_tubeids.size()];

			std::vector<int> numbers;
			for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
			{
				numbers.push_back(i);
			}
			std::random_shuffle(numbers.begin(), numbers.end());

			int change_num;
			if (context.all_segs[tube_id].size() == 1)
			{
				change_num = 1;
			}
			else
			{
				change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
			}


			for (int i = 0; i < change_num; i++)
			{

				bool print_self = false;

				/*if (itern == 43)
				{
					print_self = true;
				}*/

				int seg_id = numbers[i];

				Segment* seg = context.all_segs[tube_id][seg_id];

				if (seg->len_ < context.seg_len_lbound)	// skip segments of length lower than 15
				{
					continue;
				}

				int oaa = seg->aa_;
				int obb = seg->bb_;



				int naa, nbb;
				int C;

				if (context.remainspeed_component) {

					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000)); // rng.gaussian(1000) 从均值为0，标准差为1000的高斯分布随机采样一个数
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					/*C= rand() % 2;

					switch (C)
					{

					case 0:
						int segle;
						segle = seg->bb_ - seg->aa_;
						naa = seg->aa_ + int(rng.gaussian(1000));
						seg->aa_ = naa;
						seg->bb_ = naa + segle;
						reset_traverse(context.all_segs);
						dfs_change_aabb(seg, oaa, obb, false);
						check_consistency(context);
						break;
					case 1:
						double new_scale = random_scale(rng, seg);
						seg->scale_ = new_scale;
						reset_traverse(context.heads);
						dfs_change_aabb(seg, oaa, obb, false);
						check_consistency(context);
						break;
					}*/
				}
				else {
					C = rand() % 3;

					if (seg->colli_segs_.empty())
					{
						switch (C)
						{
						case 0:
							naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
							seg->aa_ = naa;
							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
							reset_traverse(context.heads);
							dfs_change_aabb(seg, oaa, obb, print_self);
							//std::cout << "xxxxx" << std::endl;
							check_consistency(context);
							break;
						case 1:
							nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
							seg->bb_ = nbb;
							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
							reset_traverse(context.heads);
							dfs_change_aabb(seg, oaa, obb, print_self);
							//std::cout << "xxxxx" << std::endl;
							check_consistency(context);
							break;
							//case 2:
							//	/*if (context.collisionVec[tube_id] == 0)
							//	{*/
							//	naa = random_aa(rng, seg, true);
							//	seg->aa_ = naa;
							//	reset_traverse(context.heads);
							//	dfs_change_aabb(seg, oaa, obb, print_self);
							//	//std::cout << "xxxxx" << std::endl;
							//	check_consistency(context);
							//	nbb = random_bb(rng, seg, false);
							//	seg->bb_ = nbb;
							//	//}
							//	//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
							//	break;
							//case 3:
							//	/*if (context.collisionVec[tube_id] == 0)
							//	{*/
							//	nbb = random_bb(rng, seg, true);
							//	seg->bb_ = nbb;
							//	reset_traverse(context.heads);
							//	dfs_change_aabb(seg, oaa, obb, print_self);
							//	//std::cout << "xxxxx" << std::endl;
							//	check_consistency(context);
							//	naa = random_aa(rng, seg, false);
							//	seg->aa_ = naa;
							//	//}
							//	//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
							//	break;
						case 2:
							naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
							seg->aa_ = naa;
							reset_traverse(context.heads);
							dfs_change_aabb(seg, oaa, obb, print_self);
							//std::cout << "xxxxx" << std::endl;
							check_consistency(context);

							oaa = seg->aa_; // refresh
							obb = seg->bb_;

							nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
							seg->bb_ = nbb;
							reset_traverse(context.heads);
							dfs_change_aabb(seg, oaa, obb, print_self);
							//std::cout << "xxxxx" << std::endl;
							check_consistency(context);
							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
							break;
						}
					}
					else // if seg collides, others must follow to prevent constraint broken
					{
						naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
						seg->aa_ = naa;
						reset_traverse(context.heads);
						dfs_change_aabb(seg, oaa, obb, print_self);
						//std::cout << "xxxxx" << std::endl;
						check_consistency(context);

						oaa = seg->aa_; // refresh
						obb = seg->bb_;

						nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
						seg->bb_ = nbb;
						reset_traverse(context.heads);

						dfs_change_aabb(seg, oaa, obb, print_self);
						//std::cout << "xxxxx" << std::endl;
						check_consistency(context);
						////std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
					}
				}


			}
		}


		bool allin = allin_range(context);


		double nenergy = compute_energy_chrono_reverse3(context);

		/* Very special */
		if (nenergy <= 1)
		{
			nenergy = 1.000000000001;
		}

		double eee1 = nenergy;
		double eee2 = oenergy;

		int divid = -1;
		for (int tt = 0; tt < 20; tt++)
		{
			if (eee1 > division[tt])
			{
				divid = tt;
			}
		}

		if (divid < 0 || divid >= 18)
		{
			std::cout << "Something wrong!" << std::endl;
		}

		divid++;
		eee1 = eee1 / division[divid] + 10;
		eee2 = eee2 / division[divid] + 10;

		eee1 = pow(eee1, 25);
		eee2 = pow(eee2, 25);

		double a1 = exp(-1 / hot * eee1);
		double a2 = exp(-1 / hot * eee2);

		double alpha = a1 / a2;

		alpha = alpha < 1 ? alpha : 1;

		double rd = (rand() % 1001) / 1000.0;

		if (itern % 1000 == 0)
		{
			

			std::cout << itern << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
			if (rd < alpha && allin)
				std::cout << "  Accept"<<std::endl;
			else
				std::cout << "  Refuse" << std::endl;

		}

		// 接受，更新状态
		if (rd < alpha && allin) // accept
		{
			oenergy = nenergy;

			//bool no_reverted_tube = true;

			//for (int i = 0; i < tube_num - 1; i++)
			//{
			//	for (int j = i + 1; j < tube_num; j++)
			//	{
			//		Segment* s1 = context.all_segs[i][0];
			//		Segment* s2 = context.all_segs[j][0];
			//		int ds1 = s1->a_ - s2->a_;
			//		int ds2 = s1->aa_ - s2->aa_;
			//		// 两个tube的first segment's position 逆转了，且逆转帧数超过 chrono_slack
			//		if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
			//		{
			//			no_reverted_tube = false;
			//			break;
			//		}
			//	}
			//	if (!no_reverted_tube)
			//	{
			//		break;
			//	}
			//}


			if (nenergy < best_energy)
			{
				best_energy = nenergy;

				for (size_t u = 0; u < context.all_segs.size(); u++)
				{
					for (size_t v = 0; v < context.all_segs[u].size(); v++)
					{
						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
					}
				}
			}
		}
		// 拒绝
		else
		{
			for (size_t j = 0; j < Acopy.size(); j++)
			{
				for (size_t k = 0; k < Acopy[j].size(); k++)
				{
					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
				}
			}
		}
	}

	// Delete copy
	for (size_t j = 0; j < Acopy.size(); j++)
	{
		//std::cout << j << std::endl;
		for (size_t k = 0; k < Acopy[j].size(); k++)
		{
			delete Acopy[j][k];
		}
	}

	std::cout << "Best energy:" << best_energy << std::endl;
	return iter_chrono;
}

double mcmc2(Context& context, int& iter_best, double& min_s, double& max_s, bool can_preserve_order)
{
	//std::vector<std::vector<std::tuple<cv::Vec6d, int, int>>> framebb; // used for collision energy 
	//std::vector<int> framebb_count;
	//framebb.resize(context.synoplen);
	//framebb_count.resize(context.synoplen,0);

	//for (size_t i = 0; i < framebb.size(); i++)
	//{
	//	framebb[i].resize(100);
	//}


	double enout; // energy outside
	double encolli;
	double enscale;
	std::vector<double> tube_energy;
	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };
	double sigma = context.sigma;

	int tube_num = int(context.all_segs.size());
	
	

	std::vector<std::vector<int>> stat;
	
	stat.resize(tube_num);

	for (int i = 0; i < tube_num; i++)
	{
		stat[i].resize(context.all_segs[i].size(), 0);
	}

	std::vector<std::vector<Segment*>> bestCopy;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::vector<Segment*> tmp;
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			tmp.push_back(new Segment(*context.all_segs[i][j]));
		}
		bestCopy.push_back(tmp);
	}
	
	double hot = 1.5e23;
	cv::RNG rng;

	/*tube_energy.clear();
	tube_energy.resize(context.all_segs.size(), 0);
	encolli = compute_energy_collision3(context, tube_energy, framebb, framebb_count);
	double best_energy = encolli + 1;
	double oenergy = best_energy;*/
	
	double enself;
	double oenergy = compute_energy(context, tube_energy, enout, encolli, enscale,enself);
	double best_energy = oenergy;
	double best_enco = encolli;//
	double best_enself = enself;
	

	std::vector<std::vector<Segment*>> Acopy;
	for (size_t j = 0; j < context.all_segs.size(); j++)
	{
		std::vector<Segment*> tmp;
		for (size_t k = 0; k < context.all_segs[j].size(); k++)
		{
			tmp.push_back(new Segment(*context.all_segs[j][k]));
		}
		Acopy.push_back(tmp);
	}

	
	for (int itern = 0; itern < context.mcmc2_iteration_num; itern++)
	{
		if (itern % 50000 == 0 && itern != 0)
		{
			sigma = sigma * context.sigma_decay;
			srand(time(NULL));
		}

		for (size_t j = 0; j < context.all_segs.size(); j++)
		{
			for (size_t k = 0; k < context.all_segs[j].size(); k++)
			{
				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
			}
		}

		//int tube_id = rand() % tube_num;
		//bool jumpout = false;
		int tube_id = -1;
		double max_ener = -1000;
		std::vector<int> tubes_with_collision;
		for (size_t i = 0; i < tube_energy.size(); i++)
		{
			if (tube_energy[i] > max_ener)
			{
				max_ener = tube_energy[i];
				tube_id = i;
			}
			if (tube_energy[i] > 0)
			{
				tubes_with_collision.push_back(i);
			}
		}

		if (tubes_with_collision.empty())
		{
			tube_id = 0;
		}
		else
		{
			int NNS = rand() % int(tubes_with_collision.size());
			tube_id = tubes_with_collision[NNS];
		}

		std::vector<int> numbers; 
		for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
		{
			numbers.push_back(i);
		}
		std::random_shuffle(numbers.begin(), numbers.end());

		int change_num;
		if (context.all_segs[tube_id].size() == 1)
		{
			change_num = 1;
		}
		else
		{
			change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
		}

		change_num = 1;

		for (int i = 0; i < change_num; i++)
		{
			int seg_id = numbers[i];

			stat[tube_id][seg_id] += 1;

			Segment* seg = context.all_segs[tube_id][seg_id];
			//std::cout << tube_id << "," << seg_id << std::endl;
			//if (seg->len_ < 100)	// skip segments of length 1
			//{
			//	continue;
			//}
			if (seg->len_ < context.seg_len_lbound)	// skip segments of length 1
			{
				continue;
			}
			int oaa = seg->aa_;
			int obb = seg->bb_;

			int naa, nbb;
			
			int C;
			if(context.remainspeed_component){
				C = 3 + rand() % 4;
			}else if (context.resizing_component) { //改变物体大小
				C = rand() % 7;
			}
			else
				C = rand() % 4;

			/*switch (C)
			{
			case 0:
				naa = random_aa(rng, seg, sigma, context.synoplen,min_s,max_s,context.seg_avg_len,context.bbox);
				seg->aa_ = naa;
				reset_traverse(context.heads);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 1:
				nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
				seg->bb_ = nbb;
				reset_traverse(context.all_segs);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 2:
				naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
				seg->aa_ = naa;
				reset_traverse(context.heads);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);

				oaa = seg->aa_;
				obb = seg->bb_;

				nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
				seg->bb_ = nbb;
				reset_traverse(context.all_segs);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 3:
				int segle;
				segle = seg->bb_ - seg->aa_;
				naa = seg->aa_ + int(rng.gaussian(1000));
				seg->aa_ = naa;
				seg->bb_ = naa + segle;
				reset_traverse(context.all_segs);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 4:
			case 5:
			case 6:
				double new_scale = random_scale(rng, seg);
				seg->scale_ = new_scale;
				reset_traverse(context.heads);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			}*/

			if (seg->colli_segs_.empty())
			{
				switch (C)
				{
				case 0:
					naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox,context);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 1:
					nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 2:
					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000));
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 3:
					naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					oaa = seg->aa_;
					obb = seg->bb_;

					nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 4:
				case 5:
				case 6:
					double new_scale = random_scale(rng, seg);
					seg->scale_ = new_scale;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				}
			}
			else
			{
				switch (C)
				{
				case 0:
					naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
					seg->aa_ = naa;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 1:
					nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 2:
					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000));
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 3:
					naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					oaa = seg->aa_;
					obb = seg->bb_;

					nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 4:
				case 5:
				case 6:
					double new_scale = random_scale(rng, seg);
					seg->scale_ = new_scale;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				}
			}
		}

		bool allin = true;
		double nenergy = 0.0;
		if (can_preserve_order)
		{
			allin = allin_range(context) && chronological_order_preserved(context, context.chrono_slack);
			nenergy = compute_energy(context, tube_energy, enout, encolli, enscale, enself);
		}
		else
		{
			allin = allin_range(context);
			nenergy = compute_energy(context, tube_energy, enout, encolli, enscale, enself) + compute_energy_chrono_reverse3(context);
		}
		//bool allin = allin_range(context) && chronological_order_preserved(context, context.chrono_slack);
		//bool allin = allin_range(context) && chronological_order_preserved(context, context.chrono_slack) && beg_chronological_order_preserved(context,context.beg_slack);

		//enout = compute_energy_outside(context);
		//tube_energy.clear();
		//tube_energy.resize(context.all_segs.size(), 0);
		//encolli = compute_energy_collision3(context, tube_energy, framebb, framebb_count);
		////encolli = 0;
		//double nenergy = encolli + 1;

		/* Very special */
		if (nenergy <= 1)
		{
			nenergy = 1.000000000001;
		}

		double eee1 = nenergy;
		double eee2 = oenergy;

		int divid = -1;
		for (int tt = 0; tt < 20; tt++)
		{
			if (eee1 > division[tt])
			{
				divid = tt;
			}
		}

		if (divid < 0 || divid >= 18)
		{
			std::cout << "Something wrong! divid: "<<divid <<","<< eee1 << std::endl;
			
		}

		divid++;
		eee1 = eee1 / division[divid] + 10;
		eee2 = eee2 / division[divid] + 10;

		eee1 = pow(eee1, 25);
		eee2 = pow(eee2, 25);

		double a1 = exp(-1 / hot * eee1);
		double a2 = exp(-1 / hot * eee2);

		double alpha = a1 / a2;

		alpha = alpha < 1 ? alpha : 1;

		double rd = (rand() % 1001) / 1000.0;


		if (itern % 1000 == 0)
		{
			std::cout << itern << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
			if (rd < alpha && allin)
				std::cout << "  Accept"<<std::endl;
			else
				std::cout << "  Refuse" << std::endl;
			std::cout << "Outside:" << enout << ",   Collision:" << encolli << ", Scale:" << enscale << std::endl;
			//std::cout << "Max scale:" << show_max_scale(context) << std::endl;
			std::cout << "Current Best Energy:" << best_energy << ", Best Collision:" << best_enco << ", Best Enself:" << best_enself << std::endl;
			//std::cout << std::endl;

		}

		/*如果需要满足某个范围*/
		//double max_thre = 1000000.0, min_thre = 950000.0; // 规定一个collision 范围，满足则直接提前终止mcmc2
		//if (best_enco < max_thre && best_enco > min_thre)
		//{
		//	std::cout << "break out,itern: " << itern << std::endl;
		//	break;
		//}

		if (rd < alpha && allin) // accept
		{
			oenergy = nenergy;

			/*for (size_t i = 0; i < tube_energy.size(); i++)
			{
				std::cout <<i<<": "<<tube_energy[i] << ",    ";
			}
			std::cout << std::endl;*/

			if (nenergy <= best_energy)
			{
				best_energy = nenergy;
				best_enco = encolli;
				best_enself = enself;

				for (size_t u = 0; u < context.all_segs.size(); u++)
				{
					for (size_t v = 0; v < context.all_segs[u].size(); v++)
					{
						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
					}
				}


				int tmp_energy = best_energy;
				if (tmp_energy == 1) {
					//jumpout = true;
					iter_best = itern;
					break;
				}
			}
		}
		else
		{
			for (size_t j = 0; j < Acopy.size(); j++)
			{
				for (size_t k = 0; k < Acopy[j].size(); k++)
				{
					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
				}
			}
		}
	}

	// 把 bestCopy 复制到 all_segs
	for (size_t u = 0; u < bestCopy.size(); u++)
	{
		for (size_t v = 0; v < bestCopy[u].size(); v++)
		{
			context.all_segs[u][v]->a_ = bestCopy[u][v]->a_;
			context.all_segs[u][v]->aa_ = bestCopy[u][v]->aa_;
			context.all_segs[u][v]->b_ = bestCopy[u][v]->b_;
			context.all_segs[u][v]->bb_ = bestCopy[u][v]->bb_;
		}
	}

	// Delete copy
	for (size_t j = 0; j < Acopy.size(); j++)
	{
		//std::cout << j << std::endl;
		for (size_t k = 0; k < Acopy[j].size(); k++)
		{
			delete Acopy[j][k];
		}
	}

	/*for (size_t i = 0; i < stat.size(); i++)
	{
		for (size_t j = 0; j < stat[i].size(); j++)
		{
			std::cout << i << "," << j << "," << stat[i][j] << std::endl;
			std::cout << context.all_segs[i][j]->aa_ << "," << context.all_segs[i][j]->bb_ << std::endl;
		}
	}*/
	
	std::cout << "Best energy:" << best_energy << ", Best Collision:" << best_enco << " ,Iter:" << iter_best << std::endl;

	std::ofstream ofile_colli(context.filepath + "/metrics.txt", std::ofstream::app);
	// best_enco = compute_energy_collision4(context, tube_energy, enself);
	ofile_colli << "Collision: " <<  best_enco <<"\n";

	return best_energy;
}

double mcmc_chro_colli(Context& context, int& iter_best, double& min_s, double& max_s,double limit_time)
{
	clock_t end_t, start_t = clock();

	double enout; // energy outside
	double encolli;
	double enscale;
	std::vector<double> tube_energy;
	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };
	double sigma = context.sigma;

	int tube_num = int(context.all_segs.size());

	std::vector<std::vector<int>> stat;

	stat.resize(tube_num);

	for (int i = 0; i < tube_num; i++)
	{
		stat[i].resize(context.all_segs[i].size(), 0);
	}

	std::vector<std::vector<Segment*>> bestCopy;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::vector<Segment*> tmp;
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			tmp.push_back(new Segment(*context.all_segs[i][j]));
		}
		bestCopy.push_back(tmp);
	}

	double hot = 1.5e23;
	cv::RNG rng;

	double enself;
	double oenergy = compute_energy(context, tube_energy, enout, encolli, enscale, enself) + compute_energy_chrono_reverse3(context); // 计算Ec和Et
	double best_energy = oenergy;
	double best_enco = encolli;//
	double best_enself = enself;
	double best_enout = enout;


	std::vector<std::vector<Segment*>> Acopy;
	for (size_t j = 0; j < context.all_segs.size(); j++)
	{
		std::vector<Segment*> tmp;
		for (size_t k = 0; k < context.all_segs[j].size(); k++)
		{
			tmp.push_back(new Segment(*context.all_segs[j][k]));
		}
		Acopy.push_back(tmp);
	}

	bool need_MCMC_chro = true;
	int iter_end = 0;
	for (int itern = 0; itern < context.mcmc2_iteration_num; itern++)
	{
		/*minimize Et*/
		/*if (need_MCMC_chro)
		{
			int opt = rand() % 1000;
			if (opt < 200) // the first option
			{
				std::vector<std::pair<int, int>> revert_pairs;
				std::vector<int> revert_degree;
				for (int i = 0; i < tube_num - 1; i++)
				{
					for (int j = i + 1; j < tube_num; j++)
					{
						Segment* s1 = context.all_segs[i][0];
						Segment* s2 = context.all_segs[j][0];
						int ds1 = s1->a_ - s2->a_;
						int ds2 = s1->aa_ - s2->aa_;
						// 两个tube的first segment's position 逆转了，且逆转帧数超过 chrono_slack
						if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
						{
							revert_pairs.push_back(std::pair<int, int>(i, j));
							revert_degree.push_back(abs(ds2));
						}
					}
				}
				if (revert_pairs.empty())
				{
					std::cout << "Find the solution. There is no reverted tubes." << std::endl;
					need_MCMC_chro = false;
					break;
				}
				int se_pair_id = rand() % revert_pairs.size();
				if (itern % 1000 == 0)
				{
					std::cout << "Number of reverted pairs:" << revert_pairs.size() << "............................" << std::endl;
					//std::cout << max_pair_revert << std::endl;
				}
				int tubeid1 = revert_pairs[se_pair_id].first;
				int tubeid2 = revert_pairs[se_pair_id].second;
				int t1oaa = context.all_segs[tubeid1][0]->aa_;
				int t1obb = context.all_segs[tubeid1][0]->bb_;
				int t1len = t1obb - t1oaa;
				int t2oaa = context.all_segs[tubeid2][0]->aa_;
				int t2obb = context.all_segs[tubeid2][0]->bb_;
				int t2len = t2obb - t2oaa;
				int opt2 = rand() % 2;
				if (opt2 == 0)
				{
					reset_traverse(context.heads);
					context.all_segs[tubeid1][0]->aa_ = t2oaa; // tube1 first segment's position 改为tube2的
					context.all_segs[tubeid1][0]->bb_ = t2oaa + t1len;
					dfs_change_aabb(context.all_segs[tubeid1][0], t1oaa, t1obb);
					check_consistency(context);
				}
				else
				{
					reset_traverse(context.heads);
					context.all_segs[tubeid2][0]->aa_ = t1oaa; // tube2 first segment's position 改为tube1的
					context.all_segs[tubeid2][0]->bb_ = t1oaa + t2len;
					dfs_change_aabb(context.all_segs[tubeid2][0], t2oaa, t2obb);
					check_consistency(context);
				}
			}
			else
			{
				std::vector<int> revert_tubeids;
				for (int i = 0; i < tube_num; i++)
				{
					for (int j = 0; j < tube_num; j++)
					{
						if (i == j)
							continue;
						Segment* s1 = context.all_segs[i][0];
						Segment* s2 = context.all_segs[j][0];
						int ds1 = s1->a_ - s2->a_;
						int ds2 = s1->aa_ - s2->aa_;
						if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
						{
							revert_tubeids.push_back(i);
							break;
						}
					}
				}
				if (revert_tubeids.empty())
				{
					std::cout << "Find the solution2" << std::endl;
					need_MCMC_chro = false;
					break;
				}
				int tube_id = revert_tubeids[rand() % revert_tubeids.size()];
				std::vector<int> numbers;
				for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
				{
					numbers.push_back(i);
				}
				std::random_shuffle(numbers.begin(), numbers.end());
				int change_num;
				if (context.all_segs[tube_id].size() == 1)
				{
					change_num = 1;
				}
				else
				{
					change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
				}
				for (int i = 0; i < change_num; i++)
				{
					bool print_self = false;
					int seg_id = numbers[i];
					Segment* seg = context.all_segs[tube_id][seg_id];
					if (seg->len_ < context.seg_len_lbound)	// skip segments of length lower than 15
					{
						continue;
					}
					int oaa = seg->aa_;
					int obb = seg->bb_;
					int naa, nbb;
					int C;
					if (context.remainspeed_component) {
						int segle;
						segle = seg->bb_ - seg->aa_;
						naa = seg->aa_ + int(rng.gaussian(1000)); // rng.gaussian(1000) 从均值为0，标准差为1000的高斯分布随机采样一个数
						seg->aa_ = naa;
						seg->bb_ = naa + segle;
						reset_traverse(context.all_segs);
						dfs_change_aabb(seg, oaa, obb, false);
						check_consistency(context);
					}
					else {
						C = rand() % 3;
						if (seg->colli_segs_.empty())
						{
							switch (C)
							{
							case 0:
								naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
								seg->aa_ = naa;
								//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
								reset_traverse(context.heads);
								dfs_change_aabb(seg, oaa, obb, print_self);
								check_consistency(context);
								break;
							case 1:
								nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
								seg->bb_ = nbb;
								//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
								reset_traverse(context.heads);
								dfs_change_aabb(seg, oaa, obb, print_self);
								
								check_consistency(context);
								break;
								
							case 2:
								naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
								seg->aa_ = naa;
								reset_traverse(context.heads);
								dfs_change_aabb(seg, oaa, obb, print_self);
								
								check_consistency(context);

								oaa = seg->aa_; // refresh
								obb = seg->bb_;

								nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
								seg->bb_ = nbb;
								reset_traverse(context.heads);
								dfs_change_aabb(seg, oaa, obb, print_self);
								
								check_consistency(context);
								//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
								break;
							}
						}
						else // if seg collides, others must follow to prevent constraint broken
						{
							naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
							seg->aa_ = naa;
							reset_traverse(context.heads);
							dfs_change_aabb(seg, oaa, obb, print_self);
							//std::cout << "xxxxx" << std::endl;
							check_consistency(context);

							oaa = seg->aa_; // refresh
							obb = seg->bb_;

							nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
							seg->bb_ = nbb;
							reset_traverse(context.heads);

							dfs_change_aabb(seg, oaa, obb, print_self);
							check_consistency(context);
						}
					}
				}
			}
		}*/
		
		/* minimize Ec*/
		if (itern % 50000 == 0 && itern != 0)
		{
			sigma = sigma * context.sigma_decay;
			srand(time(NULL));
		}

		for (size_t j = 0; j < context.all_segs.size(); j++)
		{
			for (size_t k = 0; k < context.all_segs[j].size(); k++)
			{
				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
			}
		}


		int tube_id = -1;
		double max_ener = -1000;
		std::vector<int> tubes_with_collision;
		for (size_t i = 0; i < tube_energy.size(); i++)
		{
			if (tube_energy[i] > max_ener)
			{
				max_ener = tube_energy[i];
				tube_id = i;
			}
			if (tube_energy[i] > 0)
			{
				tubes_with_collision.push_back(i);
			}
		}

		if (tubes_with_collision.empty())
		{
			tube_id = 0;
		}
		else
		{
			int NNS = rand() % int(tubes_with_collision.size());
			tube_id = tubes_with_collision[NNS];
		}

		std::vector<int> numbers;
		for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
		{
			numbers.push_back(i);
		}
		std::random_shuffle(numbers.begin(), numbers.end());

		int change_num;
		if (context.all_segs[tube_id].size() == 1)
		{
			change_num = 1;
		}
		else
		{
			change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
		}

		change_num = 1;

		for (int i = 0; i < change_num; i++)
		{
			int seg_id = numbers[i];

			stat[tube_id][seg_id] += 1;

			Segment* seg = context.all_segs[tube_id][seg_id];
			
			if (seg->len_ < context.seg_len_lbound)	// skip segments of length lower than lbound
			{
				continue;
			}
			int oaa = seg->aa_;
			int obb = seg->bb_;

			int naa, nbb;

			int C;
			if (context.remainspeed_component) {
				C = 3 + rand() % 4;
			}
			else if (context.resizing_component) { //改变物体大小
				C = rand() % 7;
			}
			else
				C = rand() % 4;

			if (seg->colli_segs_.empty())
			{
				switch (C)
				{
				case 0:
					naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 1:
					nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 2:
					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000));
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 3:
					naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					oaa = seg->aa_;
					obb = seg->bb_;

					nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 4:
				case 5:
				case 6:
					double new_scale = random_scale(rng, seg);
					seg->scale_ = new_scale;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				}
			}
			else
			{
				switch (C)
				{
				case 0:
					naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
					seg->aa_ = naa;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 1:
					nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 2:
					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000));
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 3:
					naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					oaa = seg->aa_;
					obb = seg->bb_;

					nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 4:
				case 5:
				case 6:
					double new_scale = random_scale(rng, seg);
					seg->scale_ = new_scale;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				}
			}
		}

		//bool allin = allin_range(context) && chronological_order_preserved(context, context.chrono_slack);
		//bool allin = allin_range(context);
		
		double nenergy = compute_energy(context, tube_energy, enout, encolli, enscale, enself) + compute_energy_chrono_reverse3(context); // 计算新的Ec和Et
		//enout = compute_energy_outside(context);
		//tube_energy.clear();
		//tube_energy.resize(context.all_segs.size(), 0);
		//encolli = compute_energy_collision3(context, tube_energy, framebb, framebb_count);
		////encolli = 0;
		//double nenergy = encolli + 1;

		/* Very special */
		if (nenergy <= 1)
		{
			nenergy = 1.000000000001;
		}

		double eee1 = nenergy;
		double eee2 = oenergy;

		int divid = -1;
		for (int tt = 0; tt < 20; tt++)
		{
			if (eee1 > division[tt])
			{
				divid = tt;
			}
		}

		if (divid < 0 || divid >= 18)
		{
			std::cout << "Something wrong! divid: " << divid << "," << eee1 << std::endl;

		}

		divid++;
		eee1 = eee1 / division[divid] + 10;
		eee2 = eee2 / division[divid] + 10;

		eee1 = pow(eee1, 25);
		eee2 = pow(eee2, 25);

		double a1 = exp(-1 / hot * eee1);
		double a2 = exp(-1 / hot * eee2);

		double alpha = a1 / a2;

		alpha = alpha < 1 ? alpha : 1;

		double rd = (rand() % 1001) / 1000.0;

		if (itern % 1000 == 0)
		{
			std::cout << itern << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
			//if (rd < alpha && allin)
			if (rd < alpha)
				std::cout << "  Accept" << std::endl;
			else
				std::cout << "  Refuse" << std::endl;
			std::cout << "Outside:" << enout << ",   Collision:" << encolli << ", Scale:" << enscale << std::endl;
			//std::cout << "Max scale:" << show_max_scale(context) << std::endl;
			std::cout << "Current Best Energy:" << best_energy << ", Best Collision:" << best_enco << ", Best Enself:" << best_enself << std::endl;
			std::cout << std::endl;
		}


		//if (rd < alpha && allin) // accept
		if (rd < alpha)
		{
			oenergy = nenergy;

			if (nenergy <= best_energy)
			{
				best_energy = nenergy;
				best_enco = encolli;
				best_enself = enself;
				best_enout = enout;

				for (size_t u = 0; u < context.all_segs.size(); u++)
				{
					for (size_t v = 0; v < context.all_segs[u].size(); v++)
					{
						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
					}
				}
			}
		}
		else
		{
			for (size_t j = 0; j < Acopy.size(); j++)
			{
				for (size_t k = 0; k < Acopy[j].size(); k++)
				{
					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
				}
			}
		}

		/*if (itern % 1000 == 0) 
		{*/
			end_t = clock();
			if (limit_time - 0.2 < (double)(end_t - start_t) / CLOCKS_PER_SEC && (double)(end_t - start_t) / CLOCKS_PER_SEC < limit_time + 0.2)
			{
				//std::cout << (double)(end_t - start_t) / CLOCKS_PER_SEC<<std::endl;
				iter_end = itern;
				break;
			}
		//}
	}

	// 把 bestCopy 复制到 all_segs
	for (size_t u = 0; u < bestCopy.size(); u++)
	{
		for (size_t v = 0; v < bestCopy[u].size(); v++)
		{
			context.all_segs[u][v]->a_ = bestCopy[u][v]->a_;
			context.all_segs[u][v]->aa_ = bestCopy[u][v]->aa_;
			context.all_segs[u][v]->b_ = bestCopy[u][v]->b_;
			context.all_segs[u][v]->bb_ = bestCopy[u][v]->bb_;
		}
	}

	// Delete copy
	for (size_t j = 0; j < Acopy.size(); j++)
	{
		for (size_t k = 0; k < Acopy[j].size(); k++)
		{
			delete Acopy[j][k];
		}
	}


	std::cout << "Best energy:" << best_energy << ", Best Collision:" << best_enco << " ,Iter:" << iter_best << std::endl;

	std::ofstream ofile_colli(context.filepath + "/metrics.txt", std::ofstream::app);
	// best_enco = compute_energy_collision4(context, tube_energy, enself);
	ofile_colli << "Collision: " << best_enco << "\n";
	ofile_colli << "Out: " << best_enout << "\n";

	return iter_end;
}

double mcmc3(Context& context, std::vector<std::vector<Segment*>>& bestCopy, int& iter_best, double& min_s, double& max_s)
{
	//std::vector<std::vector<std::tuple<cv::Vec6d, int, int>>> framebb; // used for collision energy 
	//std::vector<int> framebb_count;
	//framebb.resize(context.synoplen);
	//framebb_count.resize(context.synoplen,0);

	//for (size_t i = 0; i < framebb.size(); i++)
	//{
	//	framebb[i].resize(100);
	//}


	double enout; // energy outside
	double encolli;
	double enscale;
	std::vector<double> tube_energy;
	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };
	double sigma = context.sigma;

	int tube_num = int(context.all_segs.size());

	std::vector<std::vector<int>> stat;

	stat.resize(tube_num);

	for (int i = 0; i < tube_num; i++)
	{
		stat[i].resize(context.all_segs[i].size(), 0);
	}

	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		std::vector<Segment*> tmp;
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			tmp.push_back(new Segment(*context.all_segs[i][j]));
		}
		bestCopy.push_back(tmp);
	}

	double hot = 1.5e23;
	cv::RNG rng;

	/*tube_energy.clear();
	tube_energy.resize(context.all_segs.size(), 0);
	encolli = compute_energy_collision3(context, tube_energy, framebb, framebb_count);
	double best_energy = encolli + 1;
	double oenergy = best_energy;*/

	double enself;
	double oenergy = compute_energy_real(context, tube_energy, enout, encolli, enscale, enself);
	double best_energy = oenergy;
	double best_enco = encolli;//
	double best_enself = enself;




	std::vector<std::vector<Segment*>> Acopy;
	for (size_t j = 0; j < context.all_segs.size(); j++)
	{
		std::vector<Segment*> tmp;
		for (size_t k = 0; k < context.all_segs[j].size(); k++)
		{
			tmp.push_back(new Segment(*context.all_segs[j][k]));
		}
		Acopy.push_back(tmp);
	}

	int limit = iter_best;
	for (int itern = 0; itern < limit; itern++)
	{
		if (itern % 50000 == 0 && itern != 0)
		{
			sigma = sigma * context.sigma_decay;
			srand(time(NULL));
		}

		for (size_t j = 0; j < context.all_segs.size(); j++)
		{
			for (size_t k = 0; k < context.all_segs[j].size(); k++)
			{
				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
			}
		}

		//int tube_id = rand() % tube_num;
		//bool jumpout = false;
		int tube_id = -1;
		double max_ener = -1000;
		std::vector<int> tubes_with_collision;
		for (size_t i = 0; i < tube_energy.size(); i++)
		{
			if (tube_energy[i] > max_ener)
			{
				max_ener = tube_energy[i];
				tube_id = i;
			}
			if (tube_energy[i] > 0)
			{
				tubes_with_collision.push_back(i);
			}
		}

		if (tubes_with_collision.empty())
		{
			tube_id = 0;
		}
		else
		{
			int NNS = rand() % int(tubes_with_collision.size());
			tube_id = tubes_with_collision[NNS];
		}

		std::vector<int> numbers;
		for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
		{
			numbers.push_back(i);
		}
		std::random_shuffle(numbers.begin(), numbers.end());

		int change_num;
		if (context.all_segs[tube_id].size() == 1)
		{
			change_num = 1;
		}
		else
		{
			change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
		}

		change_num = 1;


		for (int i = 0; i < change_num; i++)
		{
			int seg_id = numbers[i];

			stat[tube_id][seg_id] += 1;

			Segment* seg = context.all_segs[tube_id][seg_id];
			//std::cout << tube_id << "," << seg_id << std::endl;
			//if (seg->len_ < 100)	// skip segments of length 1
			//{
			//	continue;
			//}
			if (seg->len_ < context.seg_len_lbound)	// skip segments of length seg_len_lbound
			{
				continue;
			}
			int oaa = seg->aa_;
			int obb = seg->bb_;

			int naa, nbb;

			int C;
			if (context.remainspeed_component) {
				C = 3 + rand() % 4;
			}
			else if (context.resizing_component) { //改变物体大小
				C = rand() % 7;
			}
			else
				C = rand() % 4;

			switch (C)
			{
			case 0:
				naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
				seg->aa_ = naa;
				reset_traverse(context.heads);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 1:
				nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
				seg->bb_ = nbb;
				reset_traverse(context.all_segs);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 2:
				naa = random_aa(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
				seg->aa_ = naa;
				reset_traverse(context.heads);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);

				oaa = seg->aa_;
				obb = seg->bb_;

				nbb = random_bb(rng, seg, sigma, context.synoplen, min_s, max_s, context.fps, context.bbox, context);
				seg->bb_ = nbb;
				reset_traverse(context.all_segs);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);


				break;
			case 3:
				int segle;
				segle = seg->bb_ - seg->aa_;
				naa = seg->aa_ + int(rng.gaussian(1000));
				seg->aa_ = naa;
				seg->bb_ = naa + segle;
				reset_traverse(context.all_segs);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			case 4:
			case 5:
			case 6:
				double new_scale = random_scale(rng, seg);
				seg->scale_ = new_scale;
				reset_traverse(context.heads);
				dfs_change_aabb(seg, oaa, obb, false);
				check_consistency(context);
				break;
			}

			/*if (seg->colli_segs_.empty())
			{
				switch (C)
				{
				case 0:
					naa = random_aa(rng, seg, sigma);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 1:
					nbb = random_bb(rng, seg, sigma);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 2:
					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000));
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 3:
					naa = random_aa(rng, seg, sigma);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					oaa = seg->aa_;
					obb = seg->bb_;

					nbb = random_bb(rng, seg, sigma);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 4:
				case 5:
				case 6:
					double new_scale = random_scale(rng, seg);
					seg->scale_ = new_scale;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				}
			}
			else
			{
				switch (C)
				{
				case 0:
					naa = random_aa_colli(rng, seg, sigma);
					seg->aa_ = naa;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 1:
					nbb = random_bb_colli(rng, seg, sigma);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 2:
					int segle;
					segle = seg->bb_ - seg->aa_;
					naa = seg->aa_ + int(rng.gaussian(1000));
					seg->aa_ = naa;
					seg->bb_ = naa + segle;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 3:
					naa = random_aa_colli(rng, seg, sigma);
					seg->aa_ = naa;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);

					oaa = seg->aa_;
					obb = seg->bb_;

					nbb = random_bb_colli(rng, seg, sigma);
					seg->bb_ = nbb;
					reset_traverse(context.all_segs);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				case 4:
				case 5:
				case 6:
					double new_scale = random_scale(rng, seg);
					seg->scale_ = new_scale;
					reset_traverse(context.heads);
					dfs_change_aabb(seg, oaa, obb, false);
					check_consistency(context);
					break;
				}
			}*/
		}

		bool allin = allin_range(context) && chronological_order_preserved(context, context.chrono_slack);
		double nenergy = compute_energy_real(context, tube_energy, enout, encolli, enscale, enself);
		//enout = compute_energy_outside(context);
		//tube_energy.clear();
		//tube_energy.resize(context.all_segs.size(), 0);
		//encolli = compute_energy_collision3(context, tube_energy, framebb, framebb_count);
		////encolli = 0;
		//double nenergy = encolli + 1;

		/* Very special */
		if (nenergy <= 1)
		{
			nenergy = 1.000000000001;
		}


		double eee1 = nenergy;
		double eee2 = oenergy;

		int divid = -1;
		for (int tt = 0; tt < 20; tt++)
		{
			if (eee1 > division[tt])
			{
				divid = tt;
			}
		}

		if (divid < 0 || divid >= 18)
		{
			std::cout << "Something wrong! divid: " << divid << "," << eee1 << std::endl;

		}

		divid++;
		eee1 = eee1 / division[divid] + 10;
		eee2 = eee2 / division[divid] + 10;

		eee1 = pow(eee1, 25);
		eee2 = pow(eee2, 25);

		double a1 = exp(-1 / hot * eee1);
		double a2 = exp(-1 / hot * eee2);

		double alpha = a1 / a2;

		alpha = alpha < 1 ? alpha : 1;

		double rd = (rand() % 1001) / 1000.0;

		if (itern % 1000 == 0)
		{
			std::cout << itern << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
			if (rd < alpha)
				std::cout << "  Accept" << std::endl;
			else
				std::cout << "  Refuse" << std::endl;
			std::cout << "Outside:" << enout << ",   Collision:" << encolli << ", Scale:" << enscale << std::endl;
			//std::cout << "Max scale:" << show_max_scale(context) << std::endl;
			std::cout << "Current Best Energy:" << best_energy << ", Best Collision:" << best_enco << ", Best Enself:" << best_enself << std::endl;

		}

		if (rd < alpha && allin) // accept
		{
			oenergy = nenergy;

			/*for (size_t i = 0; i < tube_energy.size(); i++)
			{
				std::cout <<i<<": "<<tube_energy[i] << ",    ";
			}
			std::cout << std::endl;*/

			if (nenergy < best_energy)
			{
				best_energy = nenergy;
				best_enco = encolli;
				best_enself = enself;

				for (size_t u = 0; u < context.all_segs.size(); u++)
				{
					for (size_t v = 0; v < context.all_segs[u].size(); v++)
					{
						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
					}
				}

				int tmp_energy = best_energy;
				//if (tmp_energy == 1) {
				//	//jumpout = true;
				//	iter_best = itern;
				//	break;
				//}
			}
		}
		else
		{
			for (size_t j = 0; j < Acopy.size(); j++)
			{
				for (size_t k = 0; k < Acopy[j].size(); k++)
				{
					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
				}
			}
		}
	}

	// Delete copy
	for (size_t j = 0; j < Acopy.size(); j++)
	{
		//std::cout << j << std::endl;
		for (size_t k = 0; k < Acopy[j].size(); k++)
		{
			delete Acopy[j][k];
		}
	}


	for (size_t i = 0; i < stat.size(); i++)
	{
		for (size_t j = 0; j < stat[i].size(); j++)
		{
			std::cout << i << "," << j << "," << stat[i][j] << std::endl;
			std::cout << context.all_segs[i][j]->aa_ << "," << context.all_segs[i][j]->bb_ << std::endl;
		}
	}



	std::cout << "Best energy:" << best_energy << ", Best Collision:" << best_enco << ",Iter:" << iter_best << std::endl;

	return best_energy;
}


//void mcmc_chrono(Context& context, std::vector<std::vector<Segment*>>& bestCopy,double& min_s,double& max_s)
//{
//	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };
//
//	int tube_num = int(context.all_segs.size());
//
//	for (size_t i = 0; i < context.all_segs.size(); i++)
//	{
//		std::vector<Segment*> tmp;
//		for (size_t j = 0; j < context.all_segs[i].size(); j++)
//		{
//			tmp.push_back(new Segment(*context.all_segs[i][j]));
//		}
//		bestCopy.push_back(tmp);
//	}
//
//	double hot = 1.5e23;
//	cv::RNG rng;
//
//	double oenergy = compute_energy_chrono_reverse3(context);
//	double best_energy = oenergy;
//
//
//	std::vector<std::vector<Segment*>> Acopy;
//	for (size_t j = 0; j < context.all_segs.size(); j++)
//	{
//		std::vector<Segment*> tmp;
//		for (size_t k = 0; k < context.all_segs[j].size(); k++)
//		{
//			tmp.push_back(new Segment(*context.all_segs[j][k]));
//		}
//		Acopy.push_back(tmp);
//	}
//
//
//	for (int itern = 0; itern < 1000000; itern++)
//	{
//		for (size_t j = 0; j < context.all_segs.size(); j++)
//		{
//			for (size_t k = 0; k < context.all_segs[j].size(); k++)
//			{
//				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
//				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
//			}
//		}
//
//
//		/* Two options: 
//		 * 1. Find a pair of tubes, and change their first segment' position
//		 * 2. Find a tube, and change its aabb
//		 */
//
//		
//
//
//		int opt = rand() % 1000;
//		if (opt < 200) // the first option
//		{
//			std::vector<std::pair<int, int>> revert_pairs;
//			std::vector<int> revert_degree;
//			for (int i = 0; i < tube_num - 1; i++)
//			{
//				for (int j = i + 1; j < tube_num; j++)
//				{
//					/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
//					{
//						continue;
//					}*/
//
//					Segment* s1 = context.all_segs[i][0];
//					Segment* s2 = context.all_segs[j][0];
//					int ds1 = s1->a_ - s2->a_;
//					int ds2 = s1->aa_ - s2->aa_;
//					if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
//					{
//						revert_pairs.push_back(std::pair<int, int>(i, j));
//						revert_degree.push_back(abs(ds2));
//					}
//				}
//			}
//
//			if (revert_pairs.empty())
//			{
//				std::cout << "Find the solution. There is no reverted tubes." << std::endl;
//				break;
//			}
//
//			/*int se_pair_id = -1;
//			int max_pair_revert = -100;
//			for (size_t i = 0; i < colli_pairs.size(); i++)
//			{
//				if (revert_degree[i] > max_pair_revert)
//				{
//					max_pair_revert = revert_degree[i];
//					se_pair_id = int(i);
//				}
//			}*/
//
//			int se_pair_id = rand() % revert_pairs.size();
//
//			if (itern % 1000 == 0)
//			{
//				std::cout << "Number of reverted pairs:" << revert_pairs.size() << "............................" << std::endl;
//				//std::cout << max_pair_revert << std::endl;
//			}
//
//			int tubeid1 = revert_pairs[se_pair_id].first;
//			int tubeid2 = revert_pairs[se_pair_id].second;
//
//
//			int t1oaa = context.all_segs[tubeid1][0]->aa_;
//			int t1obb = context.all_segs[tubeid1][0]->bb_;
//			int t1len = t1obb - t1oaa;
//
//			int t2oaa = context.all_segs[tubeid2][0]->aa_;
//			int t2obb = context.all_segs[tubeid2][0]->bb_;
//			int t2len = t2obb - t2oaa;
//
//			int opt2 = rand() % 2;
//
//			if (opt2 == 0)
//			{
//				reset_traverse(context.heads);
//				context.all_segs[tubeid1][0]->aa_ = t2oaa;
//				context.all_segs[tubeid1][0]->bb_ = t2oaa + t1len;
//				dfs_change_aabb(context.all_segs[tubeid1][0], t1oaa, t1obb);
//				check_consistency(context);
//			}
//			else
//			{
//				reset_traverse(context.heads);
//				context.all_segs[tubeid2][0]->aa_ = t1oaa;
//				context.all_segs[tubeid2][0]->bb_ = t1oaa + t2len;
//				dfs_change_aabb(context.all_segs[tubeid2][0], t2oaa, t2obb);
//				check_consistency(context);
//			}
//		}
//		else
//		{
//			std::vector<int> revert_tubeids;
//
//			for (int i = 0; i < tube_num; i++)
//			{
//				for (int j = 0; j < tube_num; j++)
//				{
//					if (i == j/* || context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1*/)
//						continue;
//					Segment* s1 = context.all_segs[i][0];
//					Segment* s2 = context.all_segs[j][0];
//					int ds1 = s1->a_ - s2->a_;
//					int ds2 = s1->aa_ - s2->aa_;
//					if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
//					{
//						revert_tubeids.push_back(i);
//						break;
//					}
//				}
//			}
//
//			if (revert_tubeids.empty())
//			{
//				std::cout << "Find the solution2" << std::endl;
//				break;
//			}
//
//			int tube_id = revert_tubeids[rand() % revert_tubeids.size()];
//
//			std::vector<int> numbers;
//			for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
//			{
//				numbers.push_back(i);
//			}
//			std::random_shuffle(numbers.begin(), numbers.end());
//
//			int change_num;
//			if (context.all_segs[tube_id].size() == 1)
//			{
//				change_num = 1;
//			}
//			else
//			{
//				change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
//			}
//
//
//			for (int i = 0; i < change_num; i++)
//			{
//
//				bool print_self = false;
//
//				/*if (itern == 43)
//				{
//					print_self = true;
//				}*/
//
//				int seg_id = numbers[i];
//
//				Segment* seg = context.all_segs[tube_id][seg_id];
//
//				if (seg->len_ < 10)	// skip segments of length 1
//				{
//					continue;
//				}
//
//				int oaa = seg->aa_;
//				int obb = seg->bb_;
//				
//
//
//				int naa, nbb;
//				int C;
//
//				if (context.remainspeed_component) {
//
//					int segle;
//					segle = seg->bb_ - seg->aa_;
//					naa = seg->aa_ + int(rng.gaussian(1000));
//					seg->aa_ = naa;
//					seg->bb_ = naa + segle;
//					reset_traverse(context.all_segs);
//					dfs_change_aabb(seg, oaa, obb, false);
//					check_consistency(context);
//
//					/*C= rand() % 2;
//
//					switch (C)
//					{
//					
//					case 0:
//						int segle;
//						segle = seg->bb_ - seg->aa_;
//						naa = seg->aa_ + int(rng.gaussian(1000));
//						seg->aa_ = naa;
//						seg->bb_ = naa + segle;
//						reset_traverse(context.all_segs);
//						dfs_change_aabb(seg, oaa, obb, false);
//						check_consistency(context);
//						break;
//					case 1:
//						double new_scale = random_scale(rng, seg);
//						seg->scale_ = new_scale;
//						reset_traverse(context.heads);
//						dfs_change_aabb(seg, oaa, obb, false);
//						check_consistency(context);
//						break;
//					}*/
//				}
//				else {
//					C = rand() % 3;
//
//					if (seg->colli_segs_.empty())
//					{
//						switch (C)
//						{
//						case 0:
//							naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s);
//							seg->aa_ = naa;
//							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//							break;
//						case 1:
//							nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s);
//							seg->bb_ = nbb;
//							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//							break;
//							//case 2:
//							//	/*if (context.collisionVec[tube_id] == 0)
//							//	{*/
//							//	naa = random_aa(rng, seg, true);
//							//	seg->aa_ = naa;
//							//	reset_traverse(context.heads);
//							//	dfs_change_aabb(seg, oaa, obb, print_self);
//							//	//std::cout << "xxxxx" << std::endl;
//							//	check_consistency(context);
//							//	nbb = random_bb(rng, seg, false);
//							//	seg->bb_ = nbb;
//							//	//}
//							//	//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							//	break;
//							//case 3:
//							//	/*if (context.collisionVec[tube_id] == 0)
//							//	{*/
//							//	nbb = random_bb(rng, seg, true);
//							//	seg->bb_ = nbb;
//							//	reset_traverse(context.heads);
//							//	dfs_change_aabb(seg, oaa, obb, print_self);
//							//	//std::cout << "xxxxx" << std::endl;
//							//	check_consistency(context);
//							//	naa = random_aa(rng, seg, false);
//							//	seg->aa_ = naa;
//							//	//}
//							//	//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							//	break;
//						case 2:
//							naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s);
//							seg->aa_ = naa;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//
//							oaa = seg->aa_; // refresh
//							obb = seg->bb_;
//
//							nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s);
//							seg->bb_ = nbb;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							break;
//						}
//					}
//					else // if seg collides, others must follow to prevent constraint broken
//					{
//						naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
//						seg->aa_ = naa;
//						reset_traverse(context.heads);
//						dfs_change_aabb(seg, oaa, obb, print_self);
//						//std::cout << "xxxxx" << std::endl;
//						check_consistency(context);
//
//						oaa = seg->aa_; // refresh
//						obb = seg->bb_;
//
//						nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
//						seg->bb_ = nbb;
//						reset_traverse(context.heads);
//
//						dfs_change_aabb(seg, oaa, obb, print_self);
//						//std::cout << "xxxxx" << std::endl;
//						check_consistency(context);
//						////std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//					}
//				}
//
//				
//			}
//		}
//		
//
//		bool allin = allin_range(context);
//		
//
//		double nenergy = compute_energy_chrono_reverse3(context);
//
//
//		/* Very special */
//		if (nenergy <= 1)
//		{
//			nenergy = 1.000000000001;
//		}
//
//		double eee1 = nenergy;
//		double eee2 = oenergy;
//
//		int divid = -1;
//		for (int tt = 0; tt < 20; tt++)
//		{
//			if (eee1 > division[tt])
//			{
//				divid = tt;
//			}
//		}
//
//		if (divid < 0 || divid >= 18)
//		{
//			std::cout << "Something wrong!" << std::endl;
//		}
//
//		divid++;
//		eee1 = eee1 / division[divid] + 10;
//		eee2 = eee2 / division[divid] + 10;
//
//		eee1 = pow(eee1, 25);
//		eee2 = pow(eee2, 25);
//
//		double a1 = exp(-1 / hot * eee1);
//		double a2 = exp(-1 / hot * eee2);
//
//		double alpha = a1 / a2;
//
//		alpha = alpha < 1 ? alpha : 1;
//
//		double rd = (rand() % 1001) / 1000.0;
//
//		if (itern % 10 == 0)
//		{
//			std::cout << itern << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
//			if (rd < alpha)
//				std::cout << "  Accept" << std::endl;
//			else
//				std::cout << "  Refuse" << std::endl;
//		}
//
//		if (rd < alpha && allin) // accept
//		{
//			oenergy = nenergy;
//
//			if (nenergy < best_energy)
//			{
//				best_energy = nenergy;
//
//				for (size_t u = 0; u < context.all_segs.size(); u++)
//				{
//					for (size_t v = 0; v < context.all_segs[u].size(); v++)
//					{
//						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
//						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
//						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
//						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
//					}
//				}
//			}
//		}
//		else
//		{
//			for (size_t j = 0; j < Acopy.size(); j++)
//			{
//				for (size_t k = 0; k < Acopy[j].size(); k++)
//				{
//					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
//					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
//				}
//			}
//		}
//	}
//
//	// Delete copy
//	for (size_t j = 0; j < Acopy.size(); j++)
//	{
//		//std::cout << j << std::endl;
//		for (size_t k = 0; k < Acopy[j].size(); k++)
//		{
//			delete Acopy[j][k];
//		}
//	}
//
//	std::cout << "Best energy:" << best_energy << std::endl;
//}




//int mcmc_beg_chrono(Context& context, std::vector<std::vector<Segment*>>& bestCopy, double& min_s, double& max_s)
//{
//	double division[20] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19 };
//
//	int tube_num = int(context.all_segs.size());
//
//	for (size_t i = 0; i < context.all_segs.size(); i++)
//	{
//		std::vector<Segment*> tmp;
//		for (size_t j = 0; j < context.all_segs[i].size(); j++)
//		{
//			tmp.push_back(new Segment(*context.all_segs[i][j]));
//		}
//		bestCopy.push_back(tmp);
//	}
//
//	double hot = 1.5e23;
//	cv::RNG rng;
//
//	double oenergy = compute_energy_chrono_reverse4(context);
//	double best_energy = oenergy;
//
//
//	std::vector<std::vector<Segment*>> Acopy;
//	for (size_t j = 0; j < context.all_segs.size(); j++)
//	{
//		std::vector<Segment*> tmp;
//		for (size_t k = 0; k < context.all_segs[j].size(); k++)
//		{
//			tmp.push_back(new Segment(*context.all_segs[j][k]));
//		}
//		Acopy.push_back(tmp);
//	}
//
//	int iter_chrono = 1000000;
//	for (int itern = 0; itern < 1000000; itern++)
//	{
//		for (size_t j = 0; j < context.all_segs.size(); j++)
//		{
//			for (size_t k = 0; k < context.all_segs[j].size(); k++)
//			{
//				Acopy[j][k]->aa_ = context.all_segs[j][k]->aa_;
//				Acopy[j][k]->bb_ = context.all_segs[j][k]->bb_;
//			}
//		}
//
//
//		/* Two options:
//		 * 1. Find a pair of tubes, and change their first segment's position
//		 * 2. Find a tube, and change its aabb
//		 */
//
//
//
//
//		int opt = rand() % 1000;
//		if (opt < 200) // the first option
//		{
//			std::vector<std::pair<int, int>> revert_pairs;
//			std::vector<int> revert_degree;
//			for (int i = 0; i < tube_num - 1; i++)
//			{
//				for (int j = i + 1; j < tube_num; j++)
//				{
//					/*if (context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1)
//					{
//						continue;
//					}*/
//
//					Segment* s1 = context.all_segs[i][0];
//					Segment* s2 = context.all_segs[j][0];
//					int ds1 = s1->a_ - s2->a_;
//					int ds2 = s1->aa_ - s2->aa_;
//					// 两个tube的first segment's position 逆转了，且逆转帧数超过 chrono_slack
//					if (abs(abs(ds1) - abs(ds2)) > context.beg_slack)
//					{
//						revert_pairs.push_back(std::pair<int, int>(i, j));
//						revert_degree.push_back(abs(ds2));
//					}
//				}
//			}
//
//			// 没有逆转程度超过chrono_slack的tube pair
//
//
//			//---------------
//
//			if (revert_pairs.empty())
//			{
//				std::cout << "Find the solution. There is no reverted tubes." << std::endl;
//				iter_chrono = itern + 1;
//				break;
//			}
//
//			/*if (fabs(oenergy) < 1e-15)
//			{
//				std::cout << "Find the solution. There is no reverted tubes." << std::endl;
//				iter_chrono = itern + 1;
//				break;
//			}*/
//
//
//
//			/*int se_pair_id = -1;
//			int max_pair_revert = -100;
//			for (size_t i = 0; i < colli_pairs.size(); i++)
//			{
//				if (revert_degree[i] > max_pair_revert)
//				{
//					max_pair_revert = revert_degree[i];
//					se_pair_id = int(i);
//				}
//			}*/
//
//			int se_pair_id = rand() % revert_pairs.size();
//
//			if (itern % 1000 == 0)
//			{
//				std::cout << "Number of reverted pairs:" << revert_pairs.size() << "............................" << std::endl;
//				//std::cout << max_pair_revert << std::endl;
//			}
//
//			int tubeid1 = revert_pairs[se_pair_id].first;
//			int tubeid2 = revert_pairs[se_pair_id].second;
//
//
//			int t1oaa = context.all_segs[tubeid1][0]->aa_;
//			int t1obb = context.all_segs[tubeid1][0]->bb_;
//			int t1len = t1obb - t1oaa;
//
//			int t2oaa = context.all_segs[tubeid2][0]->aa_;
//			int t2obb = context.all_segs[tubeid2][0]->bb_;
//			int t2len = t2obb - t2oaa;
//
//			int opt2 = rand() % 2;
//
//			if (opt2 == 0)
//			{
//				reset_traverse(context.heads);
//				context.all_segs[tubeid1][0]->aa_ = t2oaa; // tube1 first segment's position 改为tube2的
//				context.all_segs[tubeid1][0]->bb_ = t2oaa + t1len;
//				dfs_change_aabb(context.all_segs[tubeid1][0], t1oaa, t1obb);
//				check_consistency(context);
//			}
//			else
//			{
//				reset_traverse(context.heads);
//				context.all_segs[tubeid2][0]->aa_ = t1oaa; // tube2 first segment's position 改为tube1的
//				context.all_segs[tubeid2][0]->bb_ = t1oaa + t2len;
//				dfs_change_aabb(context.all_segs[tubeid2][0], t2oaa, t2obb);
//				check_consistency(context);
//			}
//		}
//		else
//		{
//			std::vector<int> revert_tubeids;
//
//			for (int i = 0; i < tube_num; i++)
//			{
//				for (int j = 0; j < tube_num; j++)
//				{
//					if (i == j/* || context.circulated_tube[i] == 1 || context.circulated_tube[j] == 1*/)
//						continue;
//					Segment* s1 = context.all_segs[i][0];
//					Segment* s2 = context.all_segs[j][0];
//					int ds1 = s1->a_ - s2->a_;
//					int ds2 = s1->aa_ - s2->aa_;
//					if (abs(abs(ds1) - abs(ds2)) > context.beg_slack)
//					{
//						revert_tubeids.push_back(i);
//						break;
//					}
//				}
//			}
//
//			if (revert_tubeids.empty())
//			{
//				std::cout << "Find the solution2" << std::endl;
//				iter_chrono = itern + 1;
//				break;
//			}
//
//			int tube_id = revert_tubeids[rand() % revert_tubeids.size()];
//
//			std::vector<int> numbers;
//			for (size_t i = 0; i < context.all_segs[tube_id].size(); i++)
//			{
//				numbers.push_back(i);
//			}
//			std::random_shuffle(numbers.begin(), numbers.end());
//
//			int change_num;
//			if (context.all_segs[tube_id].size() == 1)
//			{
//				change_num = 1;
//			}
//			else
//			{
//				change_num = rand() % (context.all_segs[tube_id].size() - 1) + 1;
//			}
//
//
//			for (int i = 0; i < change_num; i++)
//			{
//
//				bool print_self = false;
//
//				/*if (itern == 43)
//				{
//					print_self = true;
//				}*/
//
//				int seg_id = numbers[i];
//
//				Segment* seg = context.all_segs[tube_id][seg_id];
//
//				if (seg->len_ < context.seg_len_lbound)	// skip segments of length lower than 10
//				{
//					continue;
//				}
//
//				int oaa = seg->aa_;
//				int obb = seg->bb_;
//
//
//
//				int naa, nbb;
//				int C;
//
//				if (context.remainspeed_component) {
//
//					int segle;
//					segle = seg->bb_ - seg->aa_;
//					naa = seg->aa_ + int(rng.gaussian(1000)); // rng.gaussian(1000) 从均值为0，标准差为1000的高斯分布随机采样一个数
//					seg->aa_ = naa;
//					seg->bb_ = naa + segle;
//					reset_traverse(context.all_segs);
//					dfs_change_aabb(seg, oaa, obb, false);
//					check_consistency(context);
//
//					/*C= rand() % 2;
//
//					switch (C)
//					{
//
//					case 0:
//						int segle;
//						segle = seg->bb_ - seg->aa_;
//						naa = seg->aa_ + int(rng.gaussian(1000));
//						seg->aa_ = naa;
//						seg->bb_ = naa + segle;
//						reset_traverse(context.all_segs);
//						dfs_change_aabb(seg, oaa, obb, false);
//						check_consistency(context);
//						break;
//					case 1:
//						double new_scale = random_scale(rng, seg);
//						seg->scale_ = new_scale;
//						reset_traverse(context.heads);
//						dfs_change_aabb(seg, oaa, obb, false);
//						check_consistency(context);
//						break;
//					}*/
//				}
//				else {
//					C = rand() % 3;
//
//					if (seg->colli_segs_.empty())
//					{
//						switch (C)
//						{
//						case 0:
//							naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
//							seg->aa_ = naa;
//							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//							break;
//						case 1:
//							nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
//							seg->bb_ = nbb;
//							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//							break;
//							//case 2:
//							//	/*if (context.collisionVec[tube_id] == 0)
//							//	{*/
//							//	naa = random_aa(rng, seg, true);
//							//	seg->aa_ = naa;
//							//	reset_traverse(context.heads);
//							//	dfs_change_aabb(seg, oaa, obb, print_self);
//							//	//std::cout << "xxxxx" << std::endl;
//							//	check_consistency(context);
//							//	nbb = random_bb(rng, seg, false);
//							//	seg->bb_ = nbb;
//							//	//}
//							//	//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							//	break;
//							//case 3:
//							//	/*if (context.collisionVec[tube_id] == 0)
//							//	{*/
//							//	nbb = random_bb(rng, seg, true);
//							//	seg->bb_ = nbb;
//							//	reset_traverse(context.heads);
//							//	dfs_change_aabb(seg, oaa, obb, print_self);
//							//	//std::cout << "xxxxx" << std::endl;
//							//	check_consistency(context);
//							//	naa = random_aa(rng, seg, false);
//							//	seg->aa_ = naa;
//							//	//}
//							//	//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							//	break;
//						case 2:
//							naa = random_aa(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
//							seg->aa_ = naa;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//
//							oaa = seg->aa_; // refresh
//							obb = seg->bb_;
//
//							nbb = random_bb(rng, seg, context.sigma, context.synoplen, min_s, max_s, context.seg_avg_len, context.bbox);
//							seg->bb_ = nbb;
//							reset_traverse(context.heads);
//							dfs_change_aabb(seg, oaa, obb, print_self);
//							//std::cout << "xxxxx" << std::endl;
//							check_consistency(context);
//							//std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//							break;
//						}
//					}
//					else // if seg collides, others must follow to prevent constraint broken
//					{
//						naa = random_aa_colli(rng, seg, context.sigma, context.synoplen);
//						seg->aa_ = naa;
//						reset_traverse(context.heads);
//						dfs_change_aabb(seg, oaa, obb, print_self);
//						//std::cout << "xxxxx" << std::endl;
//						check_consistency(context);
//
//						oaa = seg->aa_; // refresh
//						obb = seg->bb_;
//
//						nbb = random_bb_colli(rng, seg, context.sigma, context.synoplen);
//						seg->bb_ = nbb;
//						reset_traverse(context.heads);
//
//						dfs_change_aabb(seg, oaa, obb, print_self);
//						//std::cout << "xxxxx" << std::endl;
//						check_consistency(context);
//						////std::cout << seg->aa_ << "," << seg->bb_ << std::endl;
//					}
//				}
//
//
//			}
//		}
//
//
//		bool allin = allin_range(context) && chronological_order_preserved(context, context.chrono_slack);
//
//
//		double nenergy = compute_energy_chrono_reverse4(context);
//
//		/* Very special */
//		if (nenergy <= 1)
//		{
//			nenergy = 1.000000000001;
//		}
//
//		double eee1 = nenergy;
//		double eee2 = oenergy;
//
//		int divid = -1;
//		for (int tt = 0; tt < 20; tt++)
//		{
//			if (eee1 > division[tt])
//			{
//				divid = tt;
//			}
//		}
//
//		if (divid < 0 || divid >= 18)
//		{
//			std::cout << "Something wrong!" << std::endl;
//		}
//
//		divid++;
//		eee1 = eee1 / division[divid] + 10;
//		eee2 = eee2 / division[divid] + 10;
//
//		eee1 = pow(eee1, 25);
//		eee2 = pow(eee2, 25);
//
//		double a1 = exp(-1 / hot * eee1);
//		double a2 = exp(-1 / hot * eee2);
//
//		double alpha = a1 / a2;
//
//		alpha = alpha < 1 ? alpha : 1;
//
//		double rd = (rand() % 1001) / 1000.0;
//
//		if (itern % 10 == 0)
//		{
//			std::cout << itern << "    old energy:" << oenergy << " new energy:" << nenergy << " a1:" << a1 << "	  a2:" << a2 << "   alpha:" << alpha;
//			if (rd < alpha && allin)
//				std::cout << "  Accept" << std::endl;
//			else
//				std::cout << "  Refuse" << std::endl;
//		}
//
//		// 接受，更新状态
//		if (rd < alpha && allin) // accept
//		{
//			oenergy = nenergy;
//
//			//bool no_reverted_tube = true;
//
//			//for (int i = 0; i < tube_num - 1; i++)
//			//{
//			//	for (int j = i + 1; j < tube_num; j++)
//			//	{
//			//		Segment* s1 = context.all_segs[i][0];
//			//		Segment* s2 = context.all_segs[j][0];
//			//		int ds1 = s1->a_ - s2->a_;
//			//		int ds2 = s1->aa_ - s2->aa_;
//			//		// 两个tube的first segment's position 逆转了，且逆转帧数超过 chrono_slack
//			//		if (ds1 * ds2 < 0 && abs(ds2) > context.chrono_slack)
//			//		{
//			//			no_reverted_tube = false;
//			//			break;
//			//		}
//			//	}
//			//	if (!no_reverted_tube)
//			//	{
//			//		break;
//			//	}
//			//}
//
//
//			if (nenergy < best_energy)
//			{
//				best_energy = nenergy;
//
//				for (size_t u = 0; u < context.all_segs.size(); u++)
//				{
//					for (size_t v = 0; v < context.all_segs[u].size(); v++)
//					{
//						bestCopy[u][v]->a_ = context.all_segs[u][v]->a_;
//						bestCopy[u][v]->aa_ = context.all_segs[u][v]->aa_;
//						bestCopy[u][v]->b_ = context.all_segs[u][v]->b_;
//						bestCopy[u][v]->bb_ = context.all_segs[u][v]->bb_;
//					}
//				}
//			}
//		}
//		// 拒绝
//		else
//		{
//			for (size_t j = 0; j < Acopy.size(); j++)
//			{
//				for (size_t k = 0; k < Acopy[j].size(); k++)
//				{
//					context.all_segs[j][k]->aa_ = Acopy[j][k]->aa_;
//					context.all_segs[j][k]->bb_ = Acopy[j][k]->bb_;
//				}
//			}
//		}
//	}
//
//	// Delete copy
//	for (size_t j = 0; j < Acopy.size(); j++)
//	{
//		//std::cout << j << std::endl;
//		for (size_t k = 0; k < Acopy[j].size(); k++)
//		{
//			delete Acopy[j][k];
//		}
//	}
//
//	std::cout << "Best energy:" << best_energy << std::endl;
//	return iter_chrono;
//}

LPCWSTR stringToLPCWSTR(std::string orig)
{
	size_t origsize = orig.length() + 1;
	//const size_t newsize = 100;
	size_t convertedChars = 0;
	wchar_t* wcstring = (wchar_t*)malloc(sizeof(wchar_t) * (orig.length() - 1));
	mbstowcs_s(&convertedChars, wcstring, origsize, orig.c_str(), _TRUNCATE);

	return wcstring;
}