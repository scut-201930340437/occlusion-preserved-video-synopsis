#include "energy.h"
#include "defines.h"

#define __DEBUG__


double compute_energy_outside(Context & context)
{
	double cost = 0;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		/*if (context.circulated_tube[i] == 1)
			continue;*/

		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			double tmp = 0;
			int aa = context.all_segs[i][j]->aa_ < 0 ? context.all_segs[i][j]->aa_ : 0;
			int bb = context.all_segs[i][j]->bb_ < 0 ? context.all_segs[i][j]->bb_ : 0;

			for (int aabb = aa; aabb < bb; aabb++)
			{
				int k = context.all_segs[i][j]->interp_frame_No(aabb);
				tmp += context.activity[i][k - context.all_segs[i][0]->a_];
			}

			aa = context.all_segs[i][j]->aa_ >= context.synoplen ? context.all_segs[i][j]->aa_ : context.synoplen;
			bb = context.all_segs[i][j]->bb_ >= context.synoplen ? context.all_segs[i][j]->bb_ : context.synoplen;

			for (int aabb = aa; aabb < bb; aabb++)
			{
				int k = context.all_segs[i][j]->interp_frame_No(aabb);
				tmp += context.activity[i][k - context.all_segs[i][0]->a_];
			}

			//tmp *= all_segs[i][j]->s_;
			cost += tmp;
		}
	}
	return cost;
}

double compute_energy_chorono_reverse2(std::vector<std::vector<Segment*>>& all_segs)
{
	double cost = 0;
	int tube_num = int(all_segs.size());
	
	for (int i = 0; i < tube_num - 1; i++)
	{
		for (int j = i + 1; j < tube_num; j++)
		{
			double ni = all_segs[i][0]->a_;
			double nj = all_segs[j][0]->b_;
			double fi = all_segs[i][0]->aa_;
			double fj = all_segs[j][0]->bb_;

			if ((ni - nj) * (fi - fj) < 0 && fabs(fi-fj) > 100)
			{
				cost += 1e10;
			}
		}
	}

	return cost;
}

double compute_energy_chrono_reverse3(Context& context)
{
	double ret = 0.0;
	for (int i = 0; i < int(context.all_segs.size())-1; i++)
	{
		for (size_t j = i + 1; j < context.all_segs.size(); j++)
		{
			Segment* s1 = context.all_segs[i][0];
			Segment* s2 = context.all_segs[j][0];

			int d1 = s1->a_ - s2->a_;
			int d2 = s1->aa_ - s2->aa_;

			if (d1 * d2 < 0 && abs(d2) > context.chrono_slack)
			{
				ret += fabs(d2) * 100.0;
			}
			//ret += std::max(abs(abs(d1) - abs(d2)) - context.beg_slack, 0) * 100.0;
		}
	}
	return ret;
}

//double compute_energy_chrono_reverse4(Context& context)
//{
//	double ret = 0.0;
//	for (int i = 0; i < int(context.all_segs.size()) - 1; i++)
//	{
//		for (size_t j = i + 1; j < context.all_segs.size(); j++)
//		{
//			Segment* s1 = context.all_segs[i][0];
//			Segment* s2 = context.all_segs[j][0];
//
//			int d1 = s1->a_ - s2->a_;
//			int d2 = s1->aa_ - s2->aa_;
//
//			/*if (d1 * d2 < 0 && abs(d2) > context.chrono_slack)
//			{
//				ret += fabs(d2) * 100.0;
//			}*/
//			ret += std::max(abs(abs(d1) - abs(d2)) - context.beg_slack, 0) * 100.0;
//		}
//	}
//
//	return ret;
//}


double compute_energy_chrono_reverse(std::vector<Segment*> & fronts)
{
	/* This part should be revised 
	 * The aim of this function is to roughly ensure appearing order of objects.
	 * Therefore, the cost can be defined directly according to each tubes, not fronts of roots
	*/
	double cost = 0;

	double xi = 150;
	double e = 1e-9;

	for (int i = 0; i < int(fronts.size()) - 1; i++)
	{
		for (size_t j = i + 1; j < fronts.size(); j++)
		{
			int aai = fronts[i]->aa_;
			int aaj = fronts[j]->aa_;

			int ai = fronts[i]->a_;
			int aj = fronts[j]->a_;

			double d1 = aai - aaj;
			double d2 = ai - aj;

			if (d1 * d2 < 0)
			{
				double omega = 1.0 / (pow((1.0 - mymin(fabs(d1), xi) / xi), 5) + e) - 1.0;
				cost += omega;
			}
		}
	}
	return cost;
}

double compute_energy_speed(std::vector<std::vector<Segment*>> & all_segs)
{
	//double cost = 0;
	//for (size_t i = 0; i < all_segs.size(); i++)
	//{
	//	double tmp = 0;
	//	for (size_t j = 0; j < all_segs[i].size(); j++)
	//	{
	//		if (all_segs[i][j]->s_ < 1.0 || all_segs[i][j]->s_>5)
	//		{
	//			tmp += 1e9;
	//		}
	//		/*double spee = all_segs[i][j]->s_;
	//		spee = exp(sig / (mymin(spee, 1.0 / spee)));
	//		tmp += spee;*/
	//	}
	//	tmp /= all_segs[i].size();
	//	cost += tmp;
	//}
	//return cost;
	return 0;
}

double compute_energy_collision2(std::vector<std::vector<Segment*>>& all_segs,
	std::vector<std::vector<cv::Vec6d>>& bbox, int synoplen, cv::Mat & collisionMat)
{
	double cost = 0;

	if (all_segs.size() < 2)
	{
		return 0;
	}

	for (size_t i = 0; i < all_segs.size() - 1; i++)
	{
		for (size_t j = i + 1; j < all_segs.size(); j++)
		{
			if (collisionMat.at<uchar>(i, j) == 1)
			{
				continue;
			}

			int iaa = all_segs[i][0]->aa_;
			int ibb = all_segs[i][all_segs[i].size() - 1]->bb_;

			int jaa = all_segs[j][0]->aa_;
			int jbb = all_segs[j][all_segs[j].size() - 1]->bb_;
			
			int minaa = mymin(iaa, jaa);
			int maxbb = mymax(ibb, jbb);

			std::vector<std::vector<cv::Vec6d>> framebb;
			framebb.resize(maxbb - minaa + 1);

			for (size_t u = 0; u < all_segs[i].size(); u++)
			{
				for (int k = all_segs[i][u]->aa_; k <= all_segs[i][u]->bb_; k++)
				{
					int frameNum = all_segs[i][u]->interp_frame_No(k);
					if (k >= 0 && k < synoplen)
					{
						framebb[k - minaa].push_back(bbox[i][frameNum - all_segs[i][0]->a_]);
					}
				}
			}
			for (size_t u = 0; u < all_segs[j].size(); u++)
			{
				for (int k = all_segs[j][u]->aa_; k <= all_segs[j][u]->bb_; k++)
				{
					int frameNum = all_segs[j][u]->interp_frame_No(k);
					if (k >= 0 && k < synoplen)
					{
						framebb[k - minaa].push_back(bbox[j][frameNum - all_segs[j][0]->a_]);
					}
				}
			}

			for (size_t k = 0; k < framebb.size(); k++)
			{
				if (framebb[k].size() == 2)
				{
					cv::Rect rt1(framebb[k][0][1], framebb[k][0][2], framebb[k][0][3], framebb[k][0][4]);
					cv::Rect rt2(framebb[k][1][1], framebb[k][1][2], framebb[k][1][3], framebb[k][1][4]);
					cv::Rect rt = rt1 & rt2;
					cost += rt.area();
				}
			}
		}
	}
	return cost;
}

double compute_energy_collision(Context & context, std::vector<double> & tube_energy)
{
	double cost = 0;

	std::vector<std::vector<std::tuple<cv::Vec6d,int, int>>> framebb;
	framebb.resize(context.synoplen);
	
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			for (int aabb = context.all_segs[i][j]->aa_; aabb <= context.all_segs[i][j]->bb_; aabb++)
			{
				int k = context.all_segs[i][j]->interp_frame_No(aabb);
				if (aabb >= 0 && aabb < context.synoplen)
				{
					framebb[aabb].push_back(std::tuple<cv::Vec6d, int, int>(context.bbox[i][k - context.all_segs[i][0]->a_], i, j));
				}
			}
		}
	}

	for (size_t i = 0; i < framebb.size(); i++)
	{
		for (int ii = 0; ii<int(framebb[i].size()) - 1; ii++)
		{
			for (size_t jj = ii + 1; jj < framebb[i].size(); jj++)
			{
				int w1, w2;
				int j1, j2;
				w1 = std::get<1>(framebb[i][ii]);
				w2 = std::get<1>(framebb[i][jj]);
				j1 = std::get<2>(framebb[i][ii]);
				j2 = std::get<2>(framebb[i][jj]);

				cv::Vec6d bbox1 = std::get<0>(framebb[i][ii]);
				cv::Vec6d bbox2 = std::get<0>(framebb[i][jj]);

				Segment* s1 = context.all_segs[w1][j1];
				Segment* s2 = context.all_segs[w2][j2];


				if (w1 == w2 || context.segCollisionMat.at<uchar>(s1->seg_index_, s2->seg_index_) == 1)
					continue;
				else
				{
					cv::Rect rt1(bbox1[1], bbox1[2], bbox1[3], bbox1[4]);
					cv::Rect rt2(bbox2[1], bbox2[2], bbox2[3], bbox2[4]);
					cv::Rect rt = rt1 & rt2;
					cost += rt.area(); 
					tube_energy[w1] += rt.area();
					tube_energy[w2] += rt.area();
				}
			}
		}
	}
	return cost;
}



double compute_energy_collision4(Context& context, std::vector<double>& tube_energy, double & cost_circuit_self)
{
	double cost = 0;

	for (int i = 0; i < int(context.all_segs.size())-1; i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			for (int ii = i+1; ii < int(context.all_segs.size()); ii++)
			{
				//if (i == ii) continue; // belonging to the same tube 
				for (size_t jj = 0; jj < context.all_segs[ii].size(); jj++)
				{

					Segment* s1 = context.all_segs[i][j];
					Segment* s2 = context.all_segs[ii][jj];
					

					if (context.segCollisionMat.at<uchar>(s1->seg_index_, s2->seg_index_) == 1) // colliding segments
						continue;

					int a = mymax(s1->aa_, s2->aa_);
					int b = mymin(s1->bb_, s2->bb_);
					
					if (a > b) continue;	// no collision artifact

					for (int aabb = a; aabb <= b; aabb+=5) //Ìø5Ö¡
					{
						if (aabb >= 0 && aabb < context.synoplen)
						{
							int k1 = s1->interp_frame_No(aabb);
							int k2 = s2->interp_frame_No(aabb);

							cv::Vec6d bbox1 = context.bbox[i][k1 - context.all_segs[i][0]->a_];
							cv::Vec6d bbox2 = context.bbox[ii][k2 - context.all_segs[ii][0]->a_];

							int dx = int(bbox1[3] * (s1->scale_ - 1) / 2.0);
							int dy = int(bbox1[4] * (s1->scale_ - 1) / 2.0);
							cv::Rect rt1(bbox1[1] - dx, bbox1[2] - dy, bbox1[3] + 2 * dx, bbox1[4] + 2 * dy);
							
							dx = int(bbox2[3] * (s2->scale_ - 1) / 2.0);
							dy = int(bbox2[4] * (s2->scale_ - 1) / 2.0);
							cv::Rect rt2(bbox2[1] - dx, bbox2[2] - dy, bbox2[3] + 2 * dx, bbox2[4] + 2 * dy);

							cv::Rect rt = rt1 & rt2;
							double aaa = rt.area();
							
							tube_energy[i] += aaa;
							tube_energy[ii] += aaa;

							cost = cost + aaa;

							bool isSameCir = false;
							for (std::set<int>::iterator iter=s1->circulation_id_.begin();iter != s1->circulation_id_.end();++iter)
							{
								if (s2->circulation_id_.count(*iter) > 0)
								{
									isSameCir = true;
									break;
								}
							}

							if (isSameCir)
							{
								cost_circuit_self = cost_circuit_self + aaa;
							}
						}
					}
				}
			}
		}
	}
	return cost;
}

double compute_energy_collision4_vis(Context& context)
{
	double cost = 0;

	for (int i = 0; i < int(context.all_segs.size()) - 1; i++)
	{
		if(!context.visualize_state[i])
		{
			continue;
		}
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			for (int ii = i + 1; ii < int(context.all_segs.size()); ii++)
			{
				if (!context.visualize_state[ii])
				{
					continue;
				}
				//if (i == ii) continue; // belonging to the same tube 
				for (size_t jj = 0; jj < context.all_segs[ii].size(); jj++)
				{

					Segment* s1 = context.all_segs[i][j];
					Segment* s2 = context.all_segs[ii][jj];


					if (context.segCollisionMat.at<uchar>(s1->seg_index_, s2->seg_index_) == 1) // colliding segments
						continue;

					int a = mymax(s1->aa_, s2->aa_);
					int b = mymin(s1->bb_, s2->bb_);

					if (a > b) continue;	// no collision artifact

					for (int aabb = a; aabb <= b; aabb += 5) //Ìø5Ö¡
					{
						if (aabb >= 0 && aabb < context.synoplen)
						{
							int k1 = s1->interp_frame_No(aabb);
							int k2 = s2->interp_frame_No(aabb);
							cv::Vec6d bbox1 = context.bbox[i][k1 - context.all_segs[i][0]->a_];
							cv::Vec6d bbox2 = context.bbox[ii][k2 - context.all_segs[ii][0]->a_];

							int dx = int(bbox1[3] * (s1->scale_ - 1) / 2.0);
							int dy = int(bbox1[4] * (s1->scale_ - 1) / 2.0);
							cv::Rect rt1(bbox1[1] - dx, bbox1[2] - dy, bbox1[3] + 2 * dx, bbox1[4] + 2 * dy);

							dx = int(bbox2[3] * (s2->scale_ - 1) / 2.0);
							dy = int(bbox2[4] * (s2->scale_ - 1) / 2.0);
							cv::Rect rt2(bbox2[1] - dx, bbox2[2] - dy, bbox2[3] + 2 * dx, bbox2[4] + 2 * dy);

							cv::Rect rt = rt1 & rt2;
							double aaa = rt.area();

							/*tube_energy[i] += aaa;
							tube_energy[ii] += aaa;*/

							cost = cost + aaa;

							bool isSameCir = false;
							for (std::set<int>::iterator iter = s1->circulation_id_.begin(); iter != s1->circulation_id_.end(); ++iter)
							{
								if (s2->circulation_id_.count(*iter) > 0)
								{
									isSameCir = true;
									break;
								}
							}

							/*if (isSameCir)
							{
								cost_circuit_self = cost_circuit_self + aaa;
							}*/
						}
					}
				}
			}
		}
	}
	return cost;
}

double compute_energy_collision5(Context& context, std::vector<double>& tube_energy, double& cost_circuit_self)
{
	double cost = 0;
	tube_energy.clear();
	tube_energy.resize(context.all_segs.size(), 0);


	for (int i = 0; i < int(context.all_segs.size()) - 1; i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{

			for (int ii = i + 1; ii < int(context.all_segs.size()); ii++)
			{
				//if (i == ii) continue; // belonging to the same tube 
				for (size_t jj = 0; jj < context.all_segs[ii].size(); jj++)
				{

					Segment* s1 = context.all_segs[i][j];
					Segment* s2 = context.all_segs[ii][jj];


					if (context.segCollisionMat.at<uchar>(s1->seg_index_, s2->seg_index_) == 1) // colliding segments
						continue;

					int a = mymax(s1->aa_, s2->aa_);
					int b = mymin(s1->bb_, s2->bb_);

					if (a > b) continue;	// no collision artifact

					for (int aabb = a; aabb <= b; aabb += 1) //Ìø1Ö¡
					{
						if (aabb >= 0 && aabb < context.synoplen)
						{
							int k1 = s1->interp_frame_No(aabb);
							int k2 = s2->interp_frame_No(aabb);
							cv::Vec6d bbox1 = context.bbox[i][k1 - context.all_segs[i][0]->a_];
							cv::Vec6d bbox2 = context.bbox[ii][k2 - context.all_segs[ii][0]->a_];

							int dx = int(bbox1[3] * (s1->scale_ - 1) / 2.0);
							int dy = int(bbox1[4] * (s1->scale_ - 1) / 2.0);
							cv::Rect rt1(bbox1[1] - dx, bbox1[2] - dy, bbox1[3] + 2 * dx, bbox1[4] + 2 * dy);

							dx = int(bbox2[3] * (s2->scale_ - 1) / 2.0);
							dy = int(bbox2[4] * (s2->scale_ - 1) / 2.0);
							cv::Rect rt2(bbox2[1] - dx, bbox2[2] - dy, bbox2[3] + 2 * dx, bbox2[4] + 2 * dy);

							cv::Rect rt = rt1 & rt2;
							double aaa = rt.area();

							tube_energy[i] += aaa;
							tube_energy[ii] += aaa;

							cost = cost + aaa;


							bool isSameCir = false;
							for (std::set<int>::iterator iter = s1->circulation_id_.begin(); iter != s1->circulation_id_.end(); ++iter)
							{
								if (s2->circulation_id_.count(*iter)>0)
								{
									isSameCir = true;
									break;
								}

							}


							if (isSameCir)
							{
								cost_circuit_self = cost_circuit_self + aaa;
							}
						}
					}
				}
			}
		}
	}

	return cost;
}

double compute_energy_collision3(Context& context, std::vector<double>& tube_energy, 
	std::vector<std::vector<std::tuple<cv::Vec6d, int, int>>>& framebb,
	std::vector<int> & frame_count)
{
	frame_count.clear();
	frame_count.resize(context.synoplen, 0);

	double cost = 0;
	
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			for (int aabb = context.all_segs[i][j]->aa_; aabb <= context.all_segs[i][j]->bb_; aabb++)
			{
				int k = context.all_segs[i][j]->interp_frame_No(aabb);
				if (aabb >= 0 && aabb < context.synoplen)
				{
					framebb[aabb][frame_count[aabb]] = (std::tuple<cv::Vec6d, int, int>(context.bbox[i][k - context.all_segs[i][0]->a_], i, j));
					frame_count[aabb] = frame_count[aabb] + 1;
				}
			}
		}
	}

	for (size_t i = 0; i < framebb.size(); i++)
	{
		for (int ii = 0; ii<frame_count[i]-1; ii++)
		{
			for (size_t jj = ii + 1; jj < frame_count[i]; jj++)
			{
				int w1, w2;
				int j1, j2;
				w1 = std::get<1>(framebb[i][ii]);
				w2 = std::get<1>(framebb[i][jj]);
				j1 = std::get<2>(framebb[i][ii]);
				j2 = std::get<2>(framebb[i][jj]);

				cv::Vec6d bbox1 = std::get<0>(framebb[i][ii]);
				cv::Vec6d bbox2 = std::get<0>(framebb[i][jj]);

				Segment* s1 = context.all_segs[w1][j1];
				Segment* s2 = context.all_segs[w2][j2];


				if (w1 == w2 || context.segCollisionMat.at<uchar>(s1->seg_index_, s2->seg_index_) == 1)
					continue;
				else
				{
					cv::Rect rt1(bbox1[1], bbox1[2], bbox1[3], bbox1[4]);
					cv::Rect rt2(bbox2[1], bbox2[2], bbox2[3], bbox2[4]);
					cv::Rect rt = rt1 & rt2;
					cost += rt.area();
					tube_energy[w1] += rt.area();
					tube_energy[w2] += rt.area();
				}
			}
		}
	}
	return cost;
}

double compute_energy_scale_size(Context& context)
{
	double cost = 0;
	const double base = 1;
	for (size_t i = 0; i < context.all_segs.size(); i++)
	{
		for (size_t j = 0; j < context.all_segs[i].size(); j++)
		{
			Segment* s = context.all_segs[i][j];
			double scale;

			double src_len = s->b_ - s->a_ + 1;
			double dst_len = s->bb_ - s->aa_ + 1;
			if (src_len < dst_len)
				scale = dst_len / src_len;
			else
				scale = src_len / dst_len;

			cost += base * scale * scale;
		}
	}
	return cost;
}

double compute_energy(Context & context, std::vector<double> & tube_energy, double & enout, double & encolli, double & enscale, double & enself)
{
	tube_energy.clear();
	tube_energy.resize(context.all_segs.size(), 0);
	enself = 0;
	
	double co = compute_energy_outside(context);

	
	//double ct = compute_energy_chrono_reverse(roots);
	//double ct = compute_energy_chorono_reverse2(context.all_segs);
	//double ct = 0;
	//double cs = compute_energy_speed(context.all_segs);
	
	//double cc2 = compute_energy_collision(context, tube_energy);
	double cc = compute_energy_collision4(context, tube_energy, enself);
	//cc = cc / 2.0;
	
	//std::cout << "**********************:" << cc << "," << cc2 << std::endl;
	//double cc = 0;


	//double cc = compute_energy_collision2(all_segs, bbox, synolen, collisionMat);
	//std::cout << "Collision cost: " << cc2 << "," << cc << "+++++++++++" << std::endl;
	
	double cs = compute_energy_scale_size(context);

	double wo = 1;
	/*double wt = 0;
	double ws = 0;*/
	double wc = 1;
	double ws = 0;

	//double total = wo * co + wt * ct + ws * cs + wc * cc + 1;
	double total = wo * co + wc * cc + ws * cs + 1;

	enout = co;
	encolli = cc;
	enscale = cs;

	//std::cout << "Outside:" << co << "                Collision:" << cc << std::endl;
#ifdef __DEBUG__
	//std::cout << "outside:" << co << "\tchrono-reverse:" << ct << "\tspeed:" << cs << "\tcollision:" << cc << "\ttotal:" << total << std::endl;
#endif // __DEBUG__

	return total;
}

double compute_energy_real(Context& context, std::vector<double>& tube_energy, double& enout, double& encolli, double& enscale, double& enself)
{
	tube_energy.clear();
	tube_energy.resize(context.all_segs.size(), 0);
	enself = 0;

	double co = compute_energy_outside(context);
	//double ct = compute_energy_chrono_reverse(roots);
	//double ct = compute_energy_chorono_reverse2(context.all_segs);
	//double ct = 0;
	//double cs = compute_energy_speed(context.all_segs);

	//double cc2 = compute_energy_collision(context, tube_energy);
	double cc = compute_energy_collision5(context, tube_energy, enself);
	//cc = cc / 2.0;

	//std::cout << "**********************:" << cc << "," << cc2 << std::endl;
	//double cc = 0;


	//double cc = compute_energy_collision2(all_segs, bbox, synolen, collisionMat);
	//std::cout << "Collision cost: " << cc2 << "," << cc << "+++++++++++" << std::endl;

	double cs = compute_energy_scale_size(context);

	double wo = 0;
	/*double wt = 0;
	double ws = 0;*/
	double wc = 1;
	double ws = 0;

	//double total = wo * co + wt * ct + ws * cs + wc * cc + 1;
	double total = wo * co + wc * cc + ws * cs + 1;

	enout = co;
	encolli = cc;
	enscale = cs;

	//std::cout << "Outside:" << co << "                Collision:" << cc << std::endl;
#ifdef __DEBUG__
	//std::cout << "outside:" << co << "\tchrono-reverse:" << ct << "\tspeed:" << cs << "\tcollision:" << cc << "\ttotal:" << total << std::endl;
#endif // __DEBUG__

	return total;
}