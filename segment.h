#pragma once
#include <vector>
#include <opencv2\opencv.hpp>

class Segment;

typedef struct EdgeNode
{
	Segment* pre_seg;
	Segment* seg;
	double weight;
	bool is_retained; // related to minimum spanning tree(最小生成树)
	bool is_traversed;
	bool is_break;
	int direction;	// edge along direction of tube or not for the computation of aa_ and bb_
	struct EdgeNode* next;
}EdgeNode;

void split_segment_using_a(Segment* seg, int pos, Segment*& seg1, Segment*& seg2);
void split_segment_using_b(Segment* seg, int pos, Segment*& seg1, Segment*& seg2);

class Segment
{
public:
	Segment(int a, int b, int n, int len, bool colli = false);
	Segment(int a, int b, int aa_, int bb_, int n, int len);
	Segment(const Segment& seg);
	~Segment();
	friend std::ostream& operator <<(std::ostream& out, const Segment& seg);

	int interp_frame_No(int aabb)
	{
		int k;
		if (aa_ == bb_)
		{
			k = a_;
		}
		else
		{
			k = (int)(((double)(aabb - aa_) / (bb_ - aa_)) * (b_ - a_) + a_);
		}

		return k;
	}

	void split_using_a(int pos, Segment* s1, Segment* s2);
	void split_using_b(int pos, Segment* s1, Segment* s2);
	void collect_right_neighbors(std::vector<Segment*>& ans);
	void collect_left_neighbors(std::vector<Segment*>& ans);

	//void MoveTo(int newPos)
	//{
	//	int oaa = this->aa_;
	//	int obb = this->bb_;
	//	int le = obb - oaa;
	//	
	//	if (direction() == 1)
	//	{
	//		this->aa_ = newPos;
	//		this->bb_ = this->aa_ + le;
	//	}
	//	else
	//	{
	//		this->bb_ = newPos;
	//		this->aa_ = this->bb_ - le;
	//	}
	//	
	//	for (size_t i = 0; i < flowout_segs_.size(); i++)
	//	{
	//		flowout_segs_[i]->change_accordingly(this, oaa, obb);
	//	}
	//}

	//void Change_speed(double newspeed)
	//{
	//	this->s_ = newspeed;
	//	int le = int(this->len_/newspeed);
	//	
	//	if (le < 1)
	//	{
	//		le = 1;
	//	}
	//	
	//	int oaa = this->aa_;
	//	int obb = this->bb_;

	//	if (direction() == 1)
	//	{
	//		this->bb_ = this->aa_ + le - 1;
	//	}
	//	else
	//	{
	//		this->aa_ = this->bb_ - le + 1;
	//	}
	//	
	//	
	//	for (size_t i = 0; i < flowout_segs_.size(); i++)
	//	{
	//		flowout_segs_[i]->change_accordingly(this, oaa, obb);
	//	}
	//}

	//int direction() // 0: backward, 1: forward
	//{
	//	if (flowin_seg_ != NULL)
	//	{
	//		if (flowin_seg_->a_ < this->a_)
	//		{
	//			return 1;
	//		}
	//		else
	//		{
	//			return 0;
	//		}
	//	}
	//	else
	//	{
	//		return 1;
	//	}
	//}

	//void change_accordingly(Segment* last_seg, int loaa, int lobb);

public:
	EdgeNode* firstedge;

	int a_;
	int b_;
	int aa_;
	int bb_;
	//double s_;	// speed
	int len_;
	int tube_id_; // begin from 0
	bool colli_;
	double scale_;

	std::vector<Segment*> colli_segs_;
	std::vector<bool> colli_break_;
	std::vector<cv::Point2i> colli_ab_;

	Segment* next_;
	Segment* last_;

	bool is_next_break_;
	bool is_last_break_;


	//std::vector<Segment*> flowout_segs_;
	//Segment* flowin_seg_;

	bool is_traversed;
	bool circulation_traversed;

	Segment* divide_seg1_;
	Segment* divide_seg2_;

	std::set<int> circulation_id_;

	int seg_index_;
};